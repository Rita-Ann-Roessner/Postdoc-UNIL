"""Build input.csv for fold.py from a TCR table with V/J gene names + CDR3.

Usage:
    python3 build_input.py <tcr_table.csv> [output.csv]

Input CSV must have columns: cdr3_TRA, cdr3_TRB, TRAV, TRAJ, TRBV, TRBJ, peptide, MHC_allele_a
Optional columns: id, species, MHC_allele_b
"""

import sys
import csv
import re
from pathlib import Path

MIXTCRVIZ_DATA = Path("/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne/MixTCRviz/data_raw/CDR123")
ALLELES_TAB = Path("/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne/preMSA/Alleles_tab.txt")


def load_alleles_tab():
    alleles = {}
    with open(ALLELES_TAB) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            alleles[row["Allele"]] = row["seq"]
    return alleles


def load_gene_table(species, chain, gene_type):
    """Load V or J gene reference table.

    Returns dict: gene_name -> {"full": str, "CDR3": str}
    """
    fname = f"TR{chain}{gene_type}.csv"
    path = MIXTCRVIZ_DATA / species / fname
    genes = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row[""]
            full = row["full"].replace("g", "").replace("-", "")
            cdr3 = row["CDR3"]
            genes[name] = {"full": full, "CDR3": cdr3}
    return genes


def get_seq_outside_cdr3(gene_type, gene_name, cdr3, gene_table):
    """Get the part of the V/J gene sequence that is outside the CDR3.

    For V: returns the sequence BEFORE the CDR3.
    For J: returns the sequence AFTER the CDR3.
    """
    if gene_name not in gene_table:
        return None

    gene_info = gene_table[gene_name]
    full_seq = gene_info["full"]
    gene_cdr3 = gene_info["CDR3"]

    if gene_type == "V":
        pattern = re.escape(gene_cdr3) + "$"
        result = re.sub(pattern, "", full_seq)
        if result == full_seq:
            for n in range(5, 0, -1):
                prefix = cdr3[:n]
                idx = full_seq.rfind(prefix)
                if idx >= 0 and len(full_seq) - idx <= 15:
                    result = full_seq[:idx]
                    break
        return result if result != full_seq else None
    else:
        pattern = "^" + re.escape(gene_cdr3)
        result = re.sub(pattern, "", full_seq)
        if result == full_seq:
            for n in range(5, 0, -1):
                suffix = cdr3[-n:]
                idx = full_seq.find(suffix)
                if idx >= 0 and idx + len(suffix) <= 15:
                    result = full_seq[idx + len(suffix):]
                    break
        return result if result != full_seq else None


def infer_species(mhc_allele):
    if mhc_allele.startswith("HLA") or re.match(r"^[ABC]\d", mhc_allele):
        return "HomoSapiens"
    elif mhc_allele.startswith("H2"):
        return "MusMusculus"
    return "HomoSapiens"


def infer_b2m(mhc_allele, species):
    if species == "MusMusculus":
        return "Mouse_B2M"
    return "B2M"


def normalize_mhc_name(name):
    if re.match(r"^[ABC]\d", name):
        return "HLA_" + name
    return name


def build_tcr_sequence(row, chain, gene_tables):
    """Build full TCR chain sequence: V_before_CDR3 + CDR3 + J_after_CDR3."""
    cdr3 = row[f"cdr3_TR{chain}"]
    v_gene = row[f"TR{chain}V"]
    j_gene = row[f"TR{chain}J"]

    v_table = gene_tables[f"TR{chain}V"]
    j_table = gene_tables[f"TR{chain}J"]

    v_seq = get_seq_outside_cdr3("V", v_gene, cdr3, v_table)
    j_seq = get_seq_outside_cdr3("J", j_gene, cdr3, j_table)

    if v_seq is None or j_seq is None:
        return None, v_gene if v_seq is None else None, j_gene if j_seq is None else None
    return v_seq + cdr3 + j_seq, None, None


def main():
    tcr_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else "input.csv"

    alleles = load_alleles_tab()

    with open(tcr_file) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print("No rows in input file.")
        sys.exit(1)

    first_mhc = normalize_mhc_name(rows[0].get("MHC_allele_a", rows[0].get("MHC", "")))
    species_default = infer_species(first_mhc)

    gene_tables = {}
    for sp in ("HomoSapiens", "MusMusculus"):
        sp_path = MIXTCRVIZ_DATA / sp
        if sp_path.exists():
            for chain in ("A", "B"):
                for gt in ("V", "J"):
                    key = f"TR{chain}{gt}"
                    table_path = sp_path / f"{key}.csv"
                    if table_path.exists():
                        gene_tables.setdefault(sp, {})[key] = load_gene_table(sp, chain, gt)

    output_rows = []
    n_bad = 0

    for i, row in enumerate(rows):
        tcr_id = row.get("id", f"tcr{i+1:04d}")

        mhc_a_name = normalize_mhc_name(row.get("MHC_allele_a", row.get("MHC", "")))
        species = row.get("species", species_default)
        if species not in gene_tables:
            print(f"  WARNING: skipping {tcr_id} — no gene data for species {species}")
            n_bad += 1
            continue

        sp_tables = gene_tables[species]

        tcra, bad_v_a, bad_j_a = build_tcr_sequence(row, "A", sp_tables)
        tcrb, bad_v_b, bad_j_b = build_tcr_sequence(row, "B", sp_tables)

        bad_genes = [g for g in [bad_v_a, bad_j_a, bad_v_b, bad_j_b] if g]
        if tcra is None or tcrb is None:
            print(f"  WARNING: skipping {tcr_id} — unknown gene(s): {', '.join(bad_genes)}")
            n_bad += 1
            continue

        if not re.match(r"^[A-Z]+$", tcra + tcrb):
            print(f"  WARNING: skipping {tcr_id} — non-standard amino acids in TCR sequence")
            n_bad += 1
            continue

        mhc_seq = alleles.get(mhc_a_name)
        if mhc_seq is None:
            print(f"  WARNING: skipping {tcr_id} — unknown MHC allele: {mhc_a_name}")
            n_bad += 1
            continue

        mhc_b_name = row.get("MHC_allele_b") or infer_b2m(mhc_a_name, species)
        if mhc_b_name.lower() == "none":
            b2m_seq = ""
        else:
            b2m_seq = alleles.get(mhc_b_name)
            if b2m_seq is None:
                print(f"  WARNING: skipping {tcr_id} — unknown B2M/MHC_b allele: {mhc_b_name}")
                n_bad += 1
                continue

        peptide = row.get("peptide", row.get("epitope", ""))
        if not peptide:
            print(f"  WARNING: skipping {tcr_id} — missing peptide sequence")
            n_bad += 1
            continue

        output_rows.append({
            "ID": tcr_id,
            "TCRA": tcra,
            "TCRB": tcrb,
            "MHC": mhc_seq,
            "B2M": b2m_seq,
            "PEPTIDE": peptide,
        })

    with open(out_file, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ID", "TCRA", "TCRB", "MHC", "B2M", "PEPTIDE"])
        w.writeheader()
        w.writerows(output_rows)

    print(f"Wrote {len(output_rows)} complexes to {out_file} ({n_bad} skipped)")


if __name__ == "__main__":
    main()
