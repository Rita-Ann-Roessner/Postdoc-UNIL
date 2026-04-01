#!/usr/bin/env python3
"""
Author : roessner <ritaann.roessner@unil.ch>
Date   : 2026-03-26
Purpose: This script fetches and processes experimental structures of MHC I / MHC II : TCR complexes from https://tcr3d.ibbr.umd.edu. Annotate with MHC + peptide + V/J segments + CDR1-3 + resolution
"""

import argparse
from typing import NamedTuple, TextIO

import glob
import os
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio.PDB import PDBList, PDBParser, NeighborSearch, PDBIO, Superimposer, Select
from Bio.SeqUtils import seq1
from Bio.PDB import Structure, Model, Chain
from Bio.Align import PairwiseAligner

import re
import requests
from difflib import SequenceMatcher


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
        df = df.rename(columns={"PDB ID":"PDB", "TRAV gene":"TRAV", "TRBV gene":"TRBV", "Epitope":"Peptide", "MHC Name":"MHC"})
        df["MHC Class"] = os.path.basename(file).split('_')[0]
        lst.append(df)

    df = pd.concat(lst)
    df = df.sort_values(by="PDB")

    # annotate species
    df["Species"] = np.where(df["MHC"].str.startswith("HLA", na=False), "HomoSapiens", np.where(df["MHC"].str.startswith("H2", na=False), "MusMusculus", pd.NA))
    df = df.dropna(subset=["Species"])

    return df
# --------------------------------------------------


# --------------------------------------------------
def annotate_cdr1_cdr2(df):
    """ Annotate CDR1 and CDR2 based on V segment """
    
    species = ['HomoSapiens', 'MusMusculus']
    genes = ['TRAV', 'TRBV']

    lst = []
    for spec in species:
        df_spec = df[df['Species'] == spec]
        for gene in genes:
            df_gene = pd.read_csv(f'/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/{spec}/{gene}.csv')
            df_gene = df_gene.rename(columns={df_gene.columns[0]:gene})
            df_gene[f'CDR1{gene[-2]}'] = df_gene['CDR1'].str.replace('-', '', regex=False)
            df_gene[f'CDR2{gene[-2]}'] = df_gene['CDR2'].str.replace('-', '', regex=False)
            df_gene[f'{gene}seq'] = df_gene['full'].str.replace('-', '', regex=False)
            df_gene = df_gene[[gene, f'CDR1{gene[-2]}', f'CDR2{gene[-2]}', f'{gene}seq']]

            df_spec = pd.merge(df_spec, df_gene, on=gene)
        lst.append(df_spec)
    
    df = pd.concat(lst)
        
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
        
        except Exception as e:
            print(f"[ERROR] Could not download {pdb_id}: {e}")
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
def longest_common_substring(s1, s2):
    m = [[0]*(1+len(s2)) for _ in range(1+len(s1))]
    longest = 0
    
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            if s1[i-1] == s2[j-1]:
                m[i][j] = m[i-1][j-1] + 1
                longest = max(longest, m[i][j])
    
    return longest
# --------------------------------------------------


# --------------------------------------------------
def find_tcr_chains_by_seq(chain_seqs, seg_seq, min_len=10):
    """ Find TCR chain by identifying overlap with gene segement. """
    hits = []
    
    seg_seq = seg_seq.lower()

    for chain_id, seq in chain_seqs.items():
        seq = seq.lower()
        
        lcs = longest_common_substring(seq, seg_seq)
        
        if lcs >= min_len:
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
    
    # if finding tcr by cdr fails (e.g. due to mutation we try to find it based on the whole seqhence)
    if len(tcr_alpha_candidates) == 0:
        tcr_alpha_candidates = find_tcr_chains_by_seq(remaining, row.TRAVseq)
    if len(tcr_beta_candidates) == 0:
        tcr_beta_candidates = find_tcr_chains_by_seq(remaining, row.TRBVseq)
    
    tcr_alpha_counts = chains_in_contact_count(structure, epitope_chain, tcr_alpha_candidates)
    tcr_beta_counts  = chains_in_contact_count(structure, epitope_chain, tcr_beta_candidates)

    # TCR alpha
    if tcr_alpha_counts and max(tcr_alpha_counts.values()) > 0:
        tcr_alpha_chain = max(tcr_alpha_counts, key=tcr_alpha_counts.get)
    elif tcr_alpha_candidates:
        alpha_distances = min_ca_distance(structure, epitope_chain, tcr_alpha_candidates)
        tcr_alpha_chain = min(alpha_distances, key=alpha_distances.get)
    else:
        tcr_alpha_chain = np.nan

    # TCR beta
    if tcr_beta_counts and max(tcr_beta_counts.values()) > 0:
        tcr_beta_chain = max(tcr_beta_counts, key=tcr_beta_counts.get)
    elif tcr_beta_candidates:
        beta_distances = min_ca_distance(structure, epitope_chain, tcr_beta_candidates)
        tcr_beta_chain = min(beta_distances, key=beta_distances.get)
    else:
        tcr_beta_chain = np.nan

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
        
        mhcii_counts = chains_in_contact_count(
            structure, epitope_chain, mhcii_candidates
        )
  
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
def clean_chains(structure, chain_mapping, out_file):
    """
    Keep only specified chains, rename them, and remove water/ions/ligands.
    """
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
            roles['pdb_id'] = pdb_id
            print('failed', roles)
            failed_pdb_ids.append(roles)
            continue

        # write specified chains to pdb file 
        clean_chains(structure, chain_mapping, out_file=os.path.join(outdir, f'{row.PDB}.pdb'))
    
    df_failed_pdb_ids = pd.DataFrame(failed_pdb_ids)
    df_failed_pdb_ids.to_csv("clean_failed.csv", index=False)
# --------------------------------------------------


# --------------------------------------------------
class OrderedResidueSelect(Select):
    def __init__(self, structure, mhc_class):
        self.keep = set()
        self.mhc_class = mhc_class
        model = structure[0]

        if mhc_class == 'classI':
            residues_A = [r for r in model['A'] if r.id[0] == ' ']
            self.keep.update(residues_A[:180])

        else:
            residues_A = [r for r in model['A'] if r.id[0] == ' ']
            residues_B = [r for r in model['B'] if r.id[0] == ' ']
            self.keep.update(residues_A[:90])
            self.keep.update(residues_B[:90])

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id

        if self.mhc_class == 'classI':
            if chain_id == 'A':
                return residue in self.keep
            return True

        else:
            if chain_id in ['A', 'B']:
                return residue in self.keep
            return True
# --------------------------------------------------
    

# --------------------------------------------------
def align_pdbs(df, indir, outdir):
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    pdbs = glob.glob(f'{indir}/*.pdb')
    pdbs.sort()

    # Save reference to output directory
    ref_file = pdbs[0]
    ref_out_path = os.path.join(outdir, os.path.basename(ref_file))
    pdb_id = os.path.basename(ref_file).split('.')[0]
    MHC_class = df.loc[df['PDB'] == pdb_id, 'MHC Class'].values[0]
    structure = parser.get_structure(os.path.basename(ref_file), ref_file)
    io.set_structure(structure)
    io.save(ref_out_path, OrderedResidueSelect(structure, MHC_class))

    # Load reference structure (first pdb)
    ref_structure = parser.get_structure("ref", ref_file)
    ref_chain = ref_structure[0]["A"]
    ref_atoms = list(chain for chain in ref_chain.get_atoms() if chain.get_name() == "CA")[:180]

    # Align all other PDBs to reference
    mobile_files = pdbs[1:]
    for pdb_path in mobile_files:
        # get MHC class
        pdb_id = os.path.basename(pdb_path).split('.')[0]
        MHC_class = df.loc[df['PDB'] == pdb_id, 'MHC Class'].values[0]

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
        io.save(out_path, OrderedResidueSelect(structure, MHC_class))
        print(f"Aligned {pdb_path} -> {out_path}")
# --------------------------------------------------


# --------------------------------------------------
def best_matching_gene(ref, chain_seq):
    """ Find gene segment with best overlap between TCR chain sequence and gene segment sequence. """
    
    aligner = PairwiseAligner()
    aligner.mode = "local"
    
    # tune these if needed
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    best_name = None
    best_score = float("-inf")
    best_full = None

    chain_seq = chain_seq.upper()

    for _, row in ref.iterrows():
        full_seq = str(row["full"]).upper()
        score = aligner.score(chain_seq, full_seq)

        if score > best_score:
            best_score = score
            best_name = row["Name"]
            best_full = full_seq

    return best_name
# --------------------------------------------------


# --------------------------------------------------
def annotate_J(df, pdb_dir):
    """ Infer J segments based on best alignment between TCR chains and reference sequences. """
    
    # reference files
    species = ['HomoSapiens', 'MusMusculus']
    gene_map = {'TRAJ':'D', 'TRBJ':'E'}

    lst=[]
    for spec in species:
        for gene in gene_map:
            df_gene = pd.read_csv(f'/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/{spec}/{gene}.csv')
            df_gene = df_gene.rename(columns={df_gene.columns[0]:'Name'})
            df_gene = df_gene[['Name', 'full']]
            df_gene['Species'] = spec
            df_gene['Gene'] = gene
            lst.append(df_gene)
    ref = pd.concat(lst)

    lst = []
    for _, row in df.iterrows():
        pdb_file = os.path.join(pdb_dir, f'{row.PDB}.pdb')

        if not os.path.exists(pdb_file):
            continue

        # get sequence for each chain
        chain_seqs = get_chain_sequences(pdb_file)

        # get TRAJ / TRBV
        for gene, chain in gene_map.items():
            tmp = ref[(ref['Species'] == row.Species) & (ref['Gene'] == gene)]
            chain_seq = chain_seqs[chain]

            best_name = best_matching_gene(tmp, chain_seq)
            row[gene] = best_name
        lst.append(row)

    df = pd.DataFrame(lst)
    return df

# --------------------------------------------------


# --------------------------------------------------
def extract_cdr3(sequence, min_len=7, max_len=23):
    """
    Extract CDR3 using conserved C (start) and the F/W-GxG motif (end). Ensures min_len <= CDR3 length <= max_len.
    """
    # find all conserved C positions
    c_positions = [m.start() for m in re.finditer("C", sequence)]
    if not c_positions:
        return None

    # try Cs from last to first (most likely real CDR3 start)
    for c_pos in reversed(c_positions):
        subseq = sequence[c_pos:]

        # find all FGXG / WGXG motifs after this C
        matches = list(re.finditer(r"([FW])G.G", subseq))
        if not matches:
            continue

        # check motifs from LAST to FIRST
        for match in reversed(matches):
            end = c_pos + match.start(1) + 1
            cdr3 = sequence[c_pos:end]

            L = len(cdr3)
            if min_len <= L <= max_len:
                return cdr3

    return None
# --------------------------------------------------


# --------------------------------------------------
def annotate_cdr3(df, pdb_dir):
    cdr_map = {'cdr3_TRA':'D', 'cdr3_TRB':'E'}

    lst = []
    for _, row in df.iterrows():
        pdb_file = os.path.join(pdb_dir, f'{row.PDB}.pdb')
        
        if not os.path.exists(pdb_file):
            continue

        # get sequence for each chain
        chain_seqs = get_chain_sequences(pdb_file)
        
        # get TRAJ / TRBV
        for cdr, chain in cdr_map.items():
            chain_seq = chain_seqs[chain]
            row[cdr] = extract_cdr3(chain_seq)

        lst.append(row)

    df = pd.DataFrame(lst)
    return df

# --------------------------------------------------


# --------------------------------------------------
def longest_common_suffix_substring(query, ref):
    query = str(query).upper()
    ref = str(ref).upper()

    max_len = 0

    # try all suffixes of query
    for i in range(len(query)):
        suffix_q = query[i:]

        # try all suffixes of ref
        for j in range(len(ref)):
            suffix_r = ref[j:]

            # compare up to min length
            k = 0
            while k < min(len(suffix_q), len(suffix_r)) and suffix_q[k] == suffix_r[k]:
                k += 1

            if k > max_len:
                max_len = k

    return max_len
# --------------------------------------------------


# --------------------------------------------------
def best_matching_gene_by_cdr3_suffix(df, cdr3_seq, cdr3_col="CDR3", name_col="Name", min_match=3):
    """
    Return the gene name whose reference CDR3 end is the longest exact suffix
    of the given cdr3_seq.
    """
    cdr3_seq = str(cdr3_seq).strip().upper()

    tmp = df[[name_col, cdr3_col]].copy()
    tmp[cdr3_col] = tmp[cdr3_col].astype(str).str.strip().str.upper()

    tmp["match_len"] = tmp[cdr3_col].apply(
        lambda x: longest_common_suffix_substring(cdr3_seq, x)
    )

    # keep only meaningful matches
    tmp = tmp[tmp["match_len"] >= min_match]

    if tmp.empty:
        return None

    return tmp.sort_values("match_len", ascending=False).iloc[0][name_col]
# --------------------------------------------------


# --------------------------------------------------
def test(df, pdb_dir):
    """ Test if based on cdr3 extracted from structure we recover the same V/J segments. """

    # reference files
    species = ['HomoSapiens', 'MusMusculus']
    genes = ['TRAJ', 'TRBJ']

    lst=[]
    for spec in species:
        for gene in genes:
            df_gene = pd.read_csv(f'/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/{spec}/{gene}.csv')
            df_gene = df_gene.rename(columns={df_gene.columns[0]:'Name'})
            df_gene = df_gene[['Name', 'CDR3']]
            df_gene['Species'] = spec
            df_gene['Gene'] = gene
            lst.append(df_gene)
    ref = pd.concat(lst)

    chains = ['TRA', 'TRB']
    lst = []

    lst = []
    for _, row in df.iterrows():

        for chain in chains:
            cdr3 = row[f'cdr3_{chain}']
            # J-segment
            tmp = ref[(ref['Species'] == row.Species) & (ref['Gene'] == f'{chain}J')]
            row[f'{chain}J_test'] = best_matching_gene_by_cdr3_suffix(tmp, cdr3)
        
        lst.append(row)
    
    df = pd.DataFrame(lst)
    return df
# --------------------------------------------------


# --------------------------------------------------
def longest_common_substring_len(a, b):
    if not isinstance(a, str) or not isinstance(b, str):
        return 0
    match = SequenceMatcher(None, a, b, autojunk=False).find_longest_match(0, len(a), 0, len(b))
    return match.size
# --------------------------------------------------


# --------------------------------------------------
def match_ref_sequence(chain_seq, ref, min_frac=0.45):
    """
    Find best matching reference row based on longest shared contiguous block.
    min_frac is relative to the shorter of the two sequences.
    """
    scores = ref['Sequence'].apply(lambda x: longest_common_substring_len(chain_seq, x))
    ref_tmp = ref.copy()
    ref_tmp['score'] = scores

    if ref_tmp.empty or ref_tmp['score'].max() == 0:
        return None

    best = ref_tmp.loc[ref_tmp['score'].idxmax()]
    
    shorter_len = min(len(chain_seq), len(best['Sequence']))
    frac = best['score'] / shorter_len if shorter_len > 0 else 0

    if frac >= min_frac:
        return best
    return None
# --------------------------------------------------


# --------------------------------------------------
def clean_mhci_annotation(df, pdb_dir, ref):
    """ Clean MHC I annotation by comparing sequence with MHC motif atlas entries (ref). """

    lst = []
    for _, row in df.iterrows():
        if row['MHC Class'] == 'classI':
            pdb_file = os.path.join(pdb_dir, f'{row.PDB}.pdb')

            # get sequence for each chain
            chain_seqs = get_chain_sequences(pdb_file)
            
            # get MHCI sequence
            chain_seq = chain_seqs.get('A')

            best = match_ref_sequence(chain_seq, ref)
            
            if best is not None:
                row['MHC'] = best['Allele']
            else:
                row['MHC'] = 'NA'
        else:
            pass

        lst.append(row)
        print('pups', row.PDB)

    df = pd.DataFrame(lst)
    return df
# --------------------------------------------------


# --------------------------------------------------
def clean_mhcii_annotation(df, pdb_dir, ref):
    """Clean MHC II annotation by comparing chain sequence with ref, allowing end overhangs."""
    chains = ['A', 'B']
    lst = []

    for _, row in df.iterrows():
        row = row.copy()

        if row['MHC Class'] == 'classII':
            pdb_file = os.path.join(pdb_dir, f"{row.PDB}.pdb")
            chain_seqs = get_chain_sequences(pdb_file)

            mhcii = []
            for chain in chains:
                chain_seq = chain_seqs.get(chain)

                best = match_ref_sequence(chain_seq, ref)
                
                if best is not None:
                    mhcii.append(best['Allele'])
                else:
                    mhcii.append('NA')

            row['MHC'] = f'{mhcii[0]}-{mhcii[1]}'

        lst.append(row)
        print('keks', row.PDB)

    df = pd.DataFrame(lst)
    return df
# --------------------------------------------------


# --------------------------------------------------
def annotate_resolution(df):

    lst = []
    for _, row in df.iterrows():
        pdb_id = row.PDB.lower()
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

        r = requests.get(url)
        r.raise_for_status()
        data = r.json()

        res = data.get("rcsb_entry_info", {}).get("resolution_combined", None)

        if res:
            resolution = res[0]  # usually a list like [2.0]
        else:
            resolution = None  # e.g. NMR structures

        row['Resolution'] = np.round(resolution,1)
        lst.append(row)
    
    df = pd.DataFrame(lst)
    return df
# --------------------------------------------------


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """
    
    args = get_args()
    pos_arg = args.positional
    print(f'input data = "{pos_arg}"')

    csv_files = glob.glob(os.path.join(pos_arg, '*.csv'))
    df = merge_data(csv_files)

    # annotate CDR1 and CDR2 based on V segment; as well as sequence for the V segment
    df = annotate_cdr1_cdr2(df)
    df.to_csv(f'{pos_arg}.csv', index=False)
    """
    # download PDBs
    #outdir = 'pdbs'
    #download_pdbs(df, outdir)

    # clean PDBs
    indir = 'pdbs'
    outdir = 'pdbs_clean'
    os.makedirs(outdir, exist_ok=True)

    # exclude non-complex structures (TCR does not engage epitope)
    exclude_pdbs = ['6P64', '6UK2', '6UK4', '6ULN', '6ULR', '6UZ1', '6VM7', '6VM8', '6VM9', '6VM9', '6VMA', '6VMC', '7L1D', '7RRG', '7SU9', '9C3E', '9EJG', '9EJH', '9EJI', '1ZGL']
    df = df[~df['PDB'].isin(exclude_pdbs)]

    # exclude structures with epitopes that have non-canonical amino acids
    exclude_pdbs = ['3D39', '3D3V', '8SHI', '9NIG']
    df = df[~df['PDB'].isin(exclude_pdbs)]

    clean_pdbs(df, indir, outdir)

    # align PDBs
    indir = 'pdbs_clean'
    outdir = 'pdbs_mhc_align'
    os.makedirs(outdir, exist_ok=True)
    align_pdbs(df, indir, outdir)
"""  
    outdir = 'pdbs_mhc_align' # remove !!!
    pdbs_aligned = [f.split('.')[0] for f in os.listdir(outdir)]
    df = df[df['PDB'].isin(pdbs_aligned)]
    #df = df[df['PDB'].isin(['2YPL', '3KXF', '6Q3S'])]
    # infer J-segment
    df = df.drop(columns=['CDR1A', 'CDR2A', 'TRAVseq', 'CDR1B', 'CDR2B', 'TRBVseq'])
    df = annotate_J(df, outdir)

    # annotate CDR3s
    df = annotate_cdr3(df, outdir)
    
    # test
    df = test(df, outdir)
    fails = df[(df['TRAJ'] != df['TRAJ_test']) | (df['TRBJ'] != df['TRBJ_test'])]

    # clean MHC annotation
    ref_seqs = pd.read_csv('/Users/roessner/Documents/PostDoc/Data/TCR_data/mhc_alleles/data_classI_MHC_I_sequences_reduced.txt', sep='\t')
    df = clean_mhci_annotation(df, outdir, ref_seqs)
    
    ref_seqs = pd.read_csv('/Users/roessner/Documents/PostDoc/Data/TCR_data/mhc_alleles/data_classII_MHC_II_sequences_reduced.txt', sep='\t')
    df = clean_mhcii_annotation(df, outdir, ref_seqs)

    # annotate resolution
    df = annotate_resolution(df)
    df = df.drop(columns=['TRAJ_test', 'TRBJ_test'])
    df.to_csv(f'{pos_arg}_structures.csv', index=False)

# --------------------------------------------------
if __name__ == '__main__':
    main()
