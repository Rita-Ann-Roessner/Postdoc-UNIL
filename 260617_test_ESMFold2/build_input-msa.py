"""Build input.csv + per-complex MSA files for fold-msa.py.

Usage:
    python3 build_input-msa.py <tcr_table.csv> [output.csv] [msas_dir]

Input CSV must have columns: cdr3_TRA, cdr3_TRB, TRAV, TRAJ, TRBV, TRBJ, peptide, MHC_allele_a
Optional columns: id, species, MHC_allele_b

Outputs:
    - output.csv (default: input.csv) with columns ID,TCRA,TCRB,MHC,B2M,PEPTIDE
    - msas/{ID}/TCRA.a3m, TCRB.a3m, MHC.a3m, B2M.a3m per complex
"""

import sys
import csv
import json
import re
import os
from pathlib import Path

MIXTCRVIZ_DATA = Path("/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne/MixTCRviz/data_raw/CDR123")
ALLELES_TAB = Path("/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne/preMSA/Alleles_tab.txt")
PREMSA_ROOT = Path("/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne/preMSA/AF3_models")


# ── Allele / gene tables ────────────────────────────────────────────────────

def load_alleles_tab():
    alleles = {}
    with open(ALLELES_TAB) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            alleles[row["Allele"]] = row["seq"]
    return alleles


def load_gene_table(species, chain, gene_type):
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


# ── MSA helpers ──────────────────────────────────────────────────────────────

def gene_to_filename(gene_name):
    """Convert gene name to pre-computed MSA filename (matches R convention)."""
    name = gene_name.replace("*", "__").replace("/", ".").replace("-", ".")
    name = name.replace(".", "_")
    name = name.lower()
    name = re.sub(r"__.*$", "", name)
    return name + "_data.json"


def convert_header(header):
    """Convert AF3 taxonomy headers to ESMFold2 key=N format."""
    if header.startswith(">query"):
        return header
    m = re.search(r"TaxID=(\d+)", header) or re.search(r"OX=(\d+)", header)
    if m:
        return header + f" key={m.group(1)}"
    return header + " key=-1"


def parse_af3_msa_text(msa_text):
    """Parse AF3 FASTA-format MSA text into list of (header, sequence) tuples."""
    entries = []
    lines = msa_text.strip().split("\n")
    header = None
    seq_parts = []
    for line in lines:
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_parts)))
            header = line
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        entries.append((header, "".join(seq_parts)))
    return entries


def load_premsa_json(json_path):
    """Load AF3 pre-computed MSA JSON, return (sequence, unpaired_entries, paired_entries)."""
    with open(json_path) as f:
        data = json.load(f)
    prot = data["sequences"][0]["protein"]
    seq = prot["sequence"]
    unpaired = parse_af3_msa_text(prot.get("unpairedMsa", ""))
    paired = parse_af3_msa_text(prot.get("pairedMsa", ""))
    return seq, unpaired, paired


def write_a3m(entries, path):
    """Write list of (header, sequence) as A3M file with ESMFold2-compatible headers."""
    with open(path, "w") as f:
        for header, seq in entries:
            f.write(convert_header(header) + "\n")
            f.write(seq + "\n")


def stitch_tcr_msa(full_tcr_seq, v_gene, j_gene, cdr3, species, gene_tables):
    """Build a stitched TCR MSA from pre-computed V and J gene MSAs.

    Returns list of (header, sequence) entries, or None if MSAs not available.
    """
    chain = "A" if v_gene.startswith("TRAV") else "B"
    v_table = gene_tables[f"TR{chain}V"]
    j_table = gene_tables[f"TR{chain}J"]

    gene_info_v = v_table.get(v_gene)
    gene_info_j = j_table.get(j_gene)
    if gene_info_v is None or gene_info_j is None:
        return None

    v_ref_full = gene_info_v["full"]
    j_ref_full = gene_info_j["full"]

    v_outside = get_seq_outside_cdr3("V", v_gene, cdr3, v_table)
    j_outside = get_seq_outside_cdr3("J", j_gene, cdr3, j_table)
    if v_outside is None or j_outside is None:
        return None

    n_trim_v = len(v_ref_full) - len(v_outside)
    n_trim_j = len(j_ref_full) - len(j_outside)
    tcr_len = len(full_tcr_seq)

    # Load pre-computed MSAs
    v_subdir = "TRAV" if chain == "A" else "TRBV"
    j_subdir = "TRAJ" if chain == "A" else "TRBJ"

    v_file = PREMSA_ROOT / species / v_subdir / gene_to_filename(v_gene)
    j_file = PREMSA_ROOT / species / j_subdir / gene_to_filename(j_gene)

    v_entries = []
    j_entries = []

    if v_file.exists():
        v_msa_seq, v_unpaired, v_paired = load_premsa_json(v_file)
        v_msa_len = len(v_msa_seq)
        n_gap_v = tcr_len + n_trim_v - v_msa_len

        for header, seq in v_unpaired + v_paired:
            if header.startswith(">query"):
                continue
            trimmed = trim_alignment_end(seq, n_trim_v)
            v_entries.append((header, trimmed + "-" * n_gap_v))

    if j_file.exists():
        j_msa_seq, j_unpaired, j_paired = load_premsa_json(j_file)
        j_msa_len = len(j_msa_seq)
        n_gap_j = tcr_len + n_trim_j - j_msa_len

        for header, seq in j_unpaired + j_paired:
            if header.startswith(">query"):
                continue
            trimmed = trim_alignment_start(seq, n_trim_j)
            j_entries.append((header, "-" * n_gap_j + trimmed))

    if not v_entries and not j_entries:
        return None

    result = [(">query", full_tcr_seq)]
    result.extend(v_entries)
    result.extend(j_entries)
    return result


def trim_alignment_end(seq, n_trim):
    """Remove n_trim match-state positions from the end of an alignment sequence.

    Match-state = uppercase letter or gap '-'. Lowercase = insertions (also removed
    when adjacent to trimmed positions).
    """
    if n_trim <= 0:
        return seq
    count = 0
    cut_pos = len(seq)
    for i in range(len(seq) - 1, -1, -1):
        ch = seq[i]
        if ch.islower() or ch == ".":
            cut_pos = i
            continue
        count += 1
        cut_pos = i
        if count >= n_trim:
            break
    return seq[:cut_pos]


def trim_alignment_start(seq, n_trim):
    """Remove n_trim match-state positions from the start of an alignment sequence."""
    if n_trim <= 0:
        return seq
    count = 0
    cut_pos = 0
    for i in range(len(seq)):
        ch = seq[i]
        if ch.islower() or ch == ".":
            cut_pos = i + 1
            continue
        count += 1
        cut_pos = i + 1
        if count >= n_trim:
            break
    return seq[cut_pos:]


# Cache for pre-computed non-TCR MSAs (MHC alleles, B2M) — converted once
_msa_cache = {}


def get_direct_msa(label, species, allele_name):
    """Load and convert a pre-computed MSA for MHC or B2M (cached)."""
    cache_key = (label, species, allele_name)
    if cache_key in _msa_cache:
        return _msa_cache[cache_key]

    if label == "MHC":
        subdir = "MHC"
        fname = gene_to_filename(allele_name)
    elif label == "B2M":
        subdir = "OTHER"
        fname = "b2m_data.json"
        if allele_name == "Mouse_B2M":
            fname = "mouse_b2m_data.json"
    else:
        return None

    path = PREMSA_ROOT / species / subdir / fname
    if not path.exists():
        _msa_cache[cache_key] = None
        return None

    msa_seq, unpaired, paired = load_premsa_json(path)
    all_entries = []
    seen_headers = set()
    for header, seq in unpaired + paired:
        if header in seen_headers:
            continue
        seen_headers.add(header)
        all_entries.append((header, seq))

    _msa_cache[cache_key] = all_entries
    return all_entries


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    tcr_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else "input.csv"
    msas_dir = Path(sys.argv[3] if len(sys.argv) > 3 else "msas")

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
    n_msa_ok = 0
    n_msa_partial = 0

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

        # ── Build MSAs for this complex ──────────────────────────────────
        complex_msa_dir = msas_dir / tcr_id
        os.makedirs(complex_msa_dir, exist_ok=True)
        msa_count = 0

        cdr3_a = row[f"cdr3_TRA"]
        cdr3_b = row[f"cdr3_TRB"]
        v_gene_a = row["TRAV"]
        j_gene_a = row["TRAJ"]
        v_gene_b = row["TRBV"]
        j_gene_b = row["TRBJ"]

        tcra_msa = stitch_tcr_msa(tcra, v_gene_a, j_gene_a, cdr3_a, species, sp_tables)
        if tcra_msa:
            write_a3m(tcra_msa, complex_msa_dir / "TCRA.a3m")
            msa_count += 1

        tcrb_msa = stitch_tcr_msa(tcrb, v_gene_b, j_gene_b, cdr3_b, species, sp_tables)
        if tcrb_msa:
            write_a3m(tcrb_msa, complex_msa_dir / "TCRB.a3m")
            msa_count += 1

        mhc_msa = get_direct_msa("MHC", species, mhc_a_name)
        if mhc_msa:
            write_a3m(mhc_msa, complex_msa_dir / "MHC.a3m")
            msa_count += 1

        if b2m_seq:
            b2m_msa = get_direct_msa("B2M", species, mhc_b_name)
            if b2m_msa:
                write_a3m(b2m_msa, complex_msa_dir / "B2M.a3m")
                msa_count += 1

        if msa_count == (4 if b2m_seq else 3):
            n_msa_ok += 1
        else:
            n_msa_partial += 1

    with open(out_file, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ID", "TCRA", "TCRB", "MHC", "B2M", "PEPTIDE"])
        w.writeheader()
        w.writerows(output_rows)

    print(f"Wrote {len(output_rows)} complexes to {out_file} ({n_bad} skipped)")
    print(f"MSAs: {n_msa_ok} complete, {n_msa_partial} partial (missing some chains)")
    print(f"MSA files in {msas_dir}/")


if __name__ == "__main__":
    main()
