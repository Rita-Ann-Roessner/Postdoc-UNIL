#!/usr/bin/env python3
"""
Author : roessner <ritaann.roessner@unil.ch>
Date   : 2026-03-05
Purpose: This script fetches and processes experimental structures of MHC I / MHC II : TCR complexes from https://tcr3d.ibbr.umd.edu. Basic analyses on the V/J and epitope statistics. 
"""

import argparse
from typing import NamedTuple, TextIO

import glob
import os
import shutil
from pathlib import Path
import pandas as pd
import numpy as np

from Bio.PDB import PDBList, PDBParser, NeighborSearch, PDBIO, Superimposer
from Bio.SeqUtils import seq1
from Bio.PDB import Structure, Model, Chain


class Args(NamedTuple):
    """ Command-line arguments """
    positional: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='TCR3d fetch and analyze',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('positional',
                        metavar='str',
                        help='Directory with CSV files from https://tcr3d.ibbr.umd.edu.')

    args = parser.parse_args()

    return Args(args.positional)
# --------------------------------------------------


# --------------------------------------------------
def merge_data(csv_files):
    """ Merge csv files into one dataframe. """

    lst = []
    for file in csv_files:
        df = pd.read_csv(file)[["PDB ID", "Epitope", "MHC Name", "TRAV gene", "TRBV gene"]]
        df = df.rename(columns={"PDB ID":"PDB", "TRAV gene":"TRAV", "TRBV gene":"TRBV"})
        df["MHC"] = os.path.basename(file).split('_')[0]
        lst.append(df)

    df = pd.concat(lst)
    df = df.sort_values(by="PDB")
    
    return df
# --------------------------------------------------


# --------------------------------------------------
def annotate_cdr1_cdr2(df):
    """ Annotate CDR1 and CDR2 based on V segment """
    
    genes = ['TRAV', 'TRBV']
    for gene in genes:
        df_gene = pd.read_csv(f'/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/HomoSapiens/{gene}.csv')
        df_gene = df_gene.rename(columns={df_gene.columns[0]:gene})
        df_gene[f'CDR1{gene[-2]}'] = df_gene['CDR1'].str.replace('-', '', regex=False)
        df_gene[f'CDR2{gene[-2]}'] = df_gene['CDR2'].str.replace('-', '', regex=False)
        df_gene = df_gene[[gene, f'CDR1{gene[-2]}', f'CDR2{gene[-2]}']]

        df = pd.merge(df, df_gene, on=gene)
        
    return df
# --------------------------------------------------


# --------------------------------------------------
def download_pdbs(df, out_dir="pdbs"):
    """
    Download a PDB structure from RCSB.
    Skips download if the file already exists.
    """
    for _, row in df.iterrows():
        pdb_id = row.PDB
        pdb_id = pdb_id.lower()
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # target file path
        pdb_file = out_dir / f"{pdb_id}.pdb"

        # skip if file already exists
        if pdb_file.exists():
            print(f"[INFO] {pdb_id} already downloaded, skipping.")
            return pdb_file

        try:
            pdbl = PDBList()
            raw_path = pdbl.retrieve_pdb_file(
                pdb_id,
                pdir=str(out_dir),
                file_format='pdb',
                overwrite=True
            )

            # rename downloaded file to standard name
            raw_path = Path(raw_path)
            raw_path.rename(pdb_file)

            print(f"[INFO] Downloaded {pdb_id}")
            return pdb_file
        
        except Exception as e:
            print(f"[ERROR] Could not download {pdb_id}: {e}")
            return None
# --------------------------------------------------


# --------------------------------------------------
def get_chain_sequences(pdb_file):
    """ 
    Create dictionary with all chain ids and the corresponding sequence. 
    """

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('tcr', pdb_file)
    
    chain_seqs = {}
    for chain in structure.get_chains():
        residues = [
            res for res in chain
            if res.id[0] == " "  # exclude hetero/water
        ]
        if not residues:
            continue
        
        seq = "".join(seq1(res.resname) for res in residues)
        chain_seqs[chain.id] = seq

    return chain_seqs
# --------------------------------------------------


# --------------------------------------------------
def find_epitope_chains(chain_seqs, epitope):
    return [
        cid for cid, seq in chain_seqs.items()
        if (seq.upper() in epitope.upper()) or (epitope.upper() in seq.upper())
    ]
# --------------------------------------------------


# --------------------------------------------------
def move_epitope_to_tmp_chain(structure, epitope_chain_id, epitope_seq):

    model = next(structure.get_models())
    old_chain = model[epitope_chain_id]

    if len(old_chain) < 50:
        return structure, epitope_chain_id
    else:
        tmp_chain = Chain.Chain('tmp')
        model.add(tmp_chain)

        residues = list(old_chain.get_residues())
        L = len(epitope_seq)

        for i in range(len(residues) - L + 1):
            window = residues[i:i+L]
            seq = "".join(seq1(res.get_resname()) for res in window)

            if seq == epitope_seq:
                for res in window:
                    old_chain.detach_child(res.id)
                    tmp_chain.add(res)
                break

        return structure, 'tmp'
# --------------------------------------------------


# --------------------------------------------------
def find_tcr_chains_by_cdr(chain_seqs, cdrs):
    hits = []
    for chain_id, seq in chain_seqs.items():
        n_hits = sum(cdr.lower() in seq.lower() for cdr in cdrs)
        if n_hits >= 2:
            hits.append(chain_id)
    return hits
# --------------------------------------------------


# --------------------------------------------------
def chains_in_contact_count(structure, ref_chain_id, candidate_chain_ids, cutoff=8.0):
    """
    Helper to get chains with maximal contacts to epitope.
    """
    model = structure[0]
    ref_chain = model[ref_chain_id]
    # unfold all C-alpha atoms in the reference chain
    ref_atoms = [atom for res in ref_chain for atom in res if atom.get_name() == 'CA']

    # all C-alpha atoms in the structure
    all_atoms = [atom for chain in model for res in chain for atom in res if atom.get_name() == 'CA']
    ns = NeighborSearch(all_atoms)

    counts = {}
    for cid in candidate_chain_ids:
        candidate_atoms = [atom for res in model[cid] for atom in res if atom.get_name() == 'CA']
        contact_atoms = set()
        for atom in candidate_atoms:
            neighbors = ns.search(atom.coord, cutoff)
            # count only neighbors in ref_atoms
            if any(n in ref_atoms for n in neighbors):
                contact_atoms.add(atom)
        counts[cid] = len(contact_atoms)

    return counts
# --------------------------------------------------


# --------------------------------------------------
def min_ca_distance(structure, ref_chain_id, candidate_chain_ids):
    """
    Helper to get chains with minimum CA distance to epitope.
    """
    model = structure[0]
    ref_chain = model[ref_chain_id]

    ref_atoms = [
        atom for res in ref_chain
        for atom in res if atom.get_name() == 'CA'
    ]

    distances = {}
    for cid in candidate_chain_ids:
        candidate_atoms = [
            atom for res in model[cid]
            for atom in res if atom.get_name() == 'CA'
        ]
        if not candidate_atoms:
            distances[cid] = np.inf
            continue

        min_dist = min(
            np.linalg.norm(a.coord - b.coord)
            for a in candidate_atoms
            for b in ref_atoms
        )
        distances[cid] = min_dist

    return distances
# --------------------------------------------------


# --------------------------------------------------
def is_beta2m(seq):
    return 85 <= len(seq) <= 120
# --------------------------------------------------


# --------------------------------------------------
def select_first_complex(structure, chain_seqs, row):
    """
    Get chain ids for first pMHC:TCR copy. 
    """

    # 1. pick epitope with lowest chain ID
    epitope_chains = find_epitope_chains(chain_seqs, row.Epitope)
    if not epitope_chains:
        return {}
    epitope_chain = sorted(epitope_chains)[0]
    
    # if epitope is same chain as TCR / MHC -> move to tmp
    structure, epitope_chain = move_epitope_to_tmp_chain(structure, epitope_chain, row.Epitope)
    if epitope_chain == 'tmp':
        epitope_chains = ['tmp']

    remaining = {cid: seq for cid, seq in chain_seqs.items() if cid != epitope_chain}

    # 2. pick TCR alpha/beta chains with most contacts (or smallest distance)
    tcr_alpha_candidates = find_tcr_chains_by_cdr(remaining, [row.CDR1A,row.CDR2A])
    tcr_beta_candidates  = find_tcr_chains_by_cdr(remaining, [row.CDR1B,row.CDR2B])
    
    tcr_alpha_counts = chains_in_contact_count(structure, epitope_chain, tcr_alpha_candidates)
    tcr_beta_counts  = chains_in_contact_count(structure, epitope_chain, tcr_beta_candidates)

    # TCR alpha
    if tcr_alpha_counts and max(tcr_alpha_counts.values()) > 0:
        tcr_alpha_chain = max(tcr_alpha_counts, key=tcr_alpha_counts.get)
    #elif tcr_alpha_candidates:
    #    alpha_distances = min_ca_distance(structure, epitope_chain, tcr_alpha_candidates)
    #    tcr_alpha_chain = min(alpha_distances, key=alpha_distances.get)
    else:
        tcr_alpha_chain = None

    # TCR beta
    if tcr_beta_counts and max(tcr_beta_counts.values()) > 0:
        tcr_beta_chain = max(tcr_beta_counts, key=tcr_beta_counts.get)
    #elif tcr_beta_candidates:
    #    beta_distances = min_ca_distance(structure, epitope_chain, tcr_beta_candidates)
    #    tcr_beta_chain = min(beta_distances, key=beta_distances.get)
    else:
        tcr_beta_chain = None

    used_ids = set(epitope_chains + tcr_alpha_candidates + tcr_beta_candidates)
    
    # 3. MHC chains: split remaining into candidates
    remaining = {cid: seq for cid, seq in remaining.items() if cid not in used_ids}

    # MHC-I: look for beta2m first
    beta2m_candidates = [cid for cid, seq in remaining.items() if is_beta2m(seq)]
    
    if beta2m_candidates:
        # MHC-I 
        mhci_alpha_candidates = [cid for cid in remaining if cid not in beta2m_candidates]
        
        # count contacts to epitope
        mhci_alpha_counts = chains_in_contact_count(structure, epitope_chain, mhci_alpha_candidates)
        mhci_alpha_chain = max(mhci_alpha_counts, key=mhci_alpha_counts.get)
        roles = {
            'epitope_chain': epitope_chain,
            'tcr_alpha_chain': tcr_alpha_chain,
            'tcr_beta_chain': tcr_beta_chain,
            'mhci_alpha_chain': mhci_alpha_chain,
        }
    else:
        # MHC-II
        mhcii_candidates = [cid for cid, seq in remaining.items()]

        if len(mhcii_candidates) < 2:
            return {}
        print(epitope_chain)
        mhcii_counts = chains_in_contact_count(
            structure, epitope_chain, mhcii_candidates
        )
        print(mhcii_counts)
        # get the two chains with the highest contact counts
        top_two = sorted(
            mhcii_counts.items(),
            key=lambda x: x[1],
            reverse=True)[:2]

        # order those two by chain ID
        chain_ids = sorted(cid for cid, _ in top_two)
        mhcii_alpha, mhcii_beta = chain_ids[0], chain_ids[1]


        roles = {
            'epitope_chain': epitope_chain,
            'tcr_alpha_chain': tcr_alpha_chain,
            'tcr_beta_chain': tcr_beta_chain,
            'mhcii_alpha_chain': mhcii_alpha,
            'mhcii_beta_chain': mhcii_beta
        }

    return roles
# --------------------------------------------------


# --------------------------------------------------
def clean_chains(pdb_file, chain_mapping, out_file):
    """
    Keep only specified chains, rename them, and remove water/ions/ligands.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    # create new structure
    new_structure = Structure.Structure(structure.id)
    new_model = Model.Model(0)
    new_structure.add(new_model)

    model = structure[0]

    for old_chain_id, new_chain_id in chain_mapping.items():
        if old_chain_id not in model:
            continue  # skip missing chains
        old_chain = model[old_chain_id]

        # create new chain
        new_chain = Chain.Chain(new_chain_id)

        for residue in old_chain:
            # skip water and hetero residues
            if residue.id[0] != ' ':
                continue
            new_chain.add(residue.copy())

        new_model.add(new_chain)

    # save cleaned PDB
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(out_file)
# --------------------------------------------------


# --------------------------------------------------
def clean_pdbs(df, indir='pdbs', outdir='pdbs_clean'):
    """
    Keep only specified chains, rename them, and remove water/ions/ligands.
    """

    failed_pdb_ids = []
    for _, row in df.iterrows():
        pdb_file = os.path.join(indir, f'{row.PDB}.pdb')

        if not os.path.exists(pdb_file):
            continue

        # get sequence for each chain
        chain_seqs = get_chain_sequences(pdb_file)

        # get chain ids for first pMHC:TCR copy 
        structure = PDBParser(QUIET=True).get_structure('structure', pdb_file)
        roles = select_first_complex(structure, chain_seqs, row)
        pdb_id = os.path.basename(pdb_file).split('.')[0]
        print(pdb_id, roles)
        # chain mappings old chain : new chain
        chain_mapping = {}
        if roles and all(pd.notna(v) for v in roles.values()):
            if "mhci_alpha_chain" in roles:
                if pd.notna(roles["mhci_alpha_chain"]):
                    chain_mapping[roles["mhci_alpha_chain"]] = 'A'
            if "mhcii_alpha_chain" in roles:
                if pd.notna(roles["mhcii_alpha_chain"]):
                    chain_mapping[roles["mhcii_alpha_chain"]] = 'A'
                if pd.notna(roles["mhcii_beta_chain"]):
                    chain_mapping[roles["mhcii_beta_chain"]] = 'B'
            if pd.notna(roles["epitope_chain"]):
                chain_mapping[roles["epitope_chain"]] = 'C'
            if pd.notna(roles["tcr_alpha_chain"]):
                chain_mapping[roles["tcr_alpha_chain"]] = 'D'
            if pd.notna(roles["tcr_beta_chain"]):
                chain_mapping[roles["tcr_beta_chain"]] = 'E'
        else:
            failed_pdb_ids.append(f'{pdb_id} : {roles}')
            continue
        
        # write specified chains to pdb file 
        clean_chains(pdb_file, chain_mapping, out_file=os.path.join(outdir, f'{row.PDB}.pdb'))
        with open('clean_failed.log', 'w') as f:
            f.write("\n".join(failed_pdb_ids))
# --------------------------------------------------


# --------------------------------------------------
def align_pdbs(df, indir, outdir):

    pdbs = glob.glob(f'{indir}/*.pdb')
    pdbs.sort()

    ref_file = pdbs[0]
    # Save reference to output directory
    ref_out_path = os.path.join(outdir, os.path.basename(ref_file))
    shutil.copy(ref_file, ref_out_path)

    parser = PDBParser(QUIET=True)
    io = PDBIO()

    # Load reference structure (first pdb)
    ref_structure = parser.get_structure("ref", ref_file)
    ref_chain = ref_structure[0]["A"]
    ref_atoms = list(chain for chain in ref_chain.get_atoms() if chain.get_name() == "CA")[:180]

    # Align all other PDBs to reference
    mobile_files = pdbs[1:]
    for pdb_path in mobile_files:
        # get MHC class
        pdb_id = os.path.basename(pdb_path).split('.')[0]
        MHC_class = df.loc[df['PDB'] == pdb_id, 'MHC'].values[0]

        structure = parser.get_structure(os.path.basename(pdb_path), pdb_path)
        if MHC_class == 'classI':
            chain = structure[0]["A"]
            moving_atoms = list(atom for atom in chain.get_atoms() if atom.get_name() == "CA")[:180]
        else:  # MHC class II
            chainA = structure[0]["A"]
            chainB = structure[0]["B"]

            # first 90 C-alpha atoms from each chain
            ca_A = [atom for atom in chainA.get_atoms() if atom.get_name() == "CA"][:90]
            ca_B = [atom for atom in chainB.get_atoms() if atom.get_name() == "CA"][:90]

            # combine them
            moving_atoms = ca_A + ca_B
        
        if len(moving_atoms) != len(ref_atoms):
            continue
            
        # Superimpose
        sup = Superimposer()
        sup.set_atoms(ref_atoms, moving_atoms)
        sup.apply(structure.get_atoms())

        # Save aligned PDB
        out_path = os.path.join(outdir, os.path.basename(pdb_path))
        io.set_structure(structure)
        io.save(out_path)
        print(f"Aligned {pdb_path} -> {out_path}")
# --------------------------------------------------


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()
    pos_arg = args.positional
    print(f'input data = "{pos_arg}"')

    csv_files = glob.glob(os.path.join(pos_arg, '*.csv'))
    df = merge_data(csv_files)
    df = df[df['MHC Name'].str.startswith('HLA-', na=False)] # remove mouse entries (the remaining mouse entries get removed by cdr annotation)

    # annotate CDR1 and CDR2 based on V segment
    df = annotate_cdr1_cdr2(df)
    #df = df[df['PDB'].isin(['1AO7', '3O6F', '3T0E', '4P4K', '6DFX', '8VD0'])]
    # download PDBs
    outdir = 'pdbs'
    #download_pdbs(df, outdir)

    # clean PDBs
    indir = 'pdbs'
    outdir = 'pdbs_clean'
    os.makedirs(outdir, exist_ok=True)
    clean_pdbs(df, indir, outdir)
    
    # align PDBs
    indir = 'pdbs_clean'
    outdir = 'pdbs_mhc_align'
    os.makedirs(outdir, exist_ok=True)
    align_pdbs(df, indir, outdir)

# --------------------------------------------------
if __name__ == '__main__':
    main()
