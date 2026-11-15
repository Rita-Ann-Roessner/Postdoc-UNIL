#!/usr/bin/env python3
"""
Assemble the master benchmark table for the PD-1 reprediction campaign.

Merges the per-species / per-method reprediction metric CSVs with the BLI
ground truth in PD1_duo_DesignSummary.csv, and augments them with *structural*
features recomputed from the predicted complex PDBs -- including cross-method
agreement features (do ESMFold2 and Protenix-v2 place the binder the same way?)
which only exist when the two methods are considered TOGETHER.

Outputs (into Analysis/data/):
  merged_human.csv   one row per design, all method scores, plot.py-compatible
  merged_mouse.csv
  master_long.csv    tidy long format: one row per (design, species, method)
  structural_raw.csv per-(design, species, method) structural features

Chain convention in every predicted complex: A = target (PD-1), B = binder.
B-factor column carries per-residue pLDDT on a 0-100 scale.
"""

import os
import sys
import itertools
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "Analysis", "data")
os.makedirs(OUT, exist_ok=True)

SUMMARY_CSV = os.path.join(ROOT, "PD1_duo_DesignSummary.csv")

# method key -> (pretty name, per-species metrics csv, MPNN dir)
METHOD_SOURCES = {
    "ESMFold2": {
        "human": os.path.join(ROOT, "ESMFold2/human/esmfold2_repredict_out"),
        "mouse": os.path.join(ROOT, "ESMFold2/mouse/esmfold2_repredict_out"),
        "csv": "reprediction_metrics.csv",
    },
    "Protenix-v2": {
        "human": os.path.join(ROOT, "Protenix-V2/human/protenix_repredict_out"),
        "mouse": os.path.join(ROOT, "Protenix-V2/mouse/protenix_repredict_out"),
        "csv": "reprediction_metrics.csv",
    },
}

TARGET_PDB = {"human": os.path.join(ROOT, "Target/hPD1.pdb"),
              "mouse": os.path.join(ROOT, "Target/mPD1.pdb")}

# The natural ligand (human PD-L1 ectodomain) is a positive control, not a
# design: 227 aa vs a 80-150 aa design median, and its complex is a real,
# PDB-deposited interface both predictors have effectively memorised. It is
# removed from every output table (not merely flagged), so it cannot leak into
# any downstream analysis; its own row is written to positive_control.csv for
# reference. Its predicted epitope is still used, before removal, as the
# reference PD-L1 footprint for the native_epitope_jaccard feature.
POSITIVE_CONTROL = "HuPDL1-Fc"

SPECIES = ["human", "mouse"]
METHODS = ["ESMFold2", "Protenix-v2"]

CONTACT_CUTOFF = 5.0      # heavy-atom cutoff defining an interface contact (A)
AA3 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}


# ---------------------------------------------------------------- PDB parsing
def parse_complex(path):
    """Return {'A': chain_dict, 'B': chain_dict} where chain_dict holds
    resids, per-residue mean pLDDT (B-factor), CA coords and heavy-atom coords.
    Only MODEL 1 is read (these files carry a single model)."""
    heavy = {"A": [], "B": []}          # (resid, x, y, z)
    ca = {"A": {}, "B": {}}             # resid -> xyz
    plddt = {"A": {}, "B": {}}          # resid -> list of B-factors
    seq = {"A": {}, "B": {}}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("ENDMDL"):
                break
            if not ln.startswith("ATOM"):
                continue
            ch = ln[21]
            if ch not in ("A", "B"):
                continue
            elem = ln[76:78].strip() or ln[12:16].strip()[0]
            if elem == "H":
                continue
            resid = int(ln[22:26])
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            heavy[ch].append((resid,) + xyz)
            plddt[ch].setdefault(resid, []).append(float(ln[60:66]))
            if ln[12:16].strip() == "CA":
                ca[ch][resid] = xyz
                seq[ch][resid] = AA3.get(ln[17:20].strip(), "X")
    out = {}
    for ch in ("A", "B"):
        arr = np.array(heavy[ch], dtype=float)
        out[ch] = {
            "heavy_resid": arr[:, 0].astype(int) if len(arr) else np.array([], int),
            "heavy_xyz": arr[:, 1:] if len(arr) else np.zeros((0, 3)),
            "ca": ca[ch],
            "plddt": {r: float(np.mean(v)) for r, v in plddt[ch].items()},
            "seq": seq[ch],
        }
    return out


def interface_features(cx):
    """Heavy-atom contact map across the A/B interface + interface pLDDT."""
    A, B = cx["A"], cx["B"]
    if len(A["heavy_xyz"]) == 0 or len(B["heavy_xyz"]) == 0:
        return None
    tree = cKDTree(A["heavy_xyz"])
    pairs = tree.query_ball_point(B["heavy_xyz"], r=CONTACT_CUTOFF)
    n_contacts = 0
    ep_t, ep_b = set(), set()
    for i, hits in enumerate(pairs):
        if not hits:
            continue
        n_contacts += len(hits)
        ep_b.add(int(B["heavy_resid"][i]))
        for j in hits:
            ep_t.add(int(A["heavy_resid"][j]))

    def mean_plddt(chain, resids):
        vals = [chain["plddt"][r] for r in resids if r in chain["plddt"]]
        return float(np.mean(vals)) if vals else np.nan

    # centre of mass of the epitope, used for coarse epitope-shift distances
    ep_t_xyz = np.array([A["ca"][r] for r in sorted(ep_t) if r in A["ca"]])
    return {
        "n_contacts": n_contacts,
        "iface_nres_target": len(ep_t),
        "iface_nres_binder": len(ep_b),
        "iface_pLDDT_binder": mean_plddt(B, ep_b),
        "iface_pLDDT_target": mean_plddt(A, ep_t),
        "binder_pLDDT": float(np.mean(list(B["plddt"].values()))) if B["plddt"] else np.nan,
        "target_pLDDT": float(np.mean(list(A["plddt"].values()))) if A["plddt"] else np.nan,
        "epitope": frozenset(ep_t),
        "paratope": frozenset(ep_b),
        "epitope_com": ep_t_xyz.mean(axis=0) if len(ep_t_xyz) else None,
    }


# -------------------------------------------------------------- superposition
def kabsch_rmsd(P, Q, apply_to=None):
    """Superpose P onto Q (both N x 3, paired). Returns (rmsd_on_P_Q,
    transformed `apply_to`)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    P0, Q0 = P - Pc, Q - Qc
    V, S, Wt = np.linalg.svd(P0.T @ Q0)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    P_rot = P0 @ R
    rmsd = float(np.sqrt(((P_rot - Q0) ** 2).sum() / len(P)))
    moved = None
    if apply_to is not None:
        moved = (apply_to - Pc) @ R + Qc
    return rmsd, moved


def jaccard(a, b):
    if not a and not b:
        return np.nan
    return len(a & b) / len(a | b)


# --------------------------------------------------- cross-species h<->m map
def species_residue_map(ref_complex):
    """Align the hPD1 and mPD1 target chains so epitopes can be compared across
    species. hPD1 and mPD1 are only ~65% identical and carry an indel, so a
    constant numbering offset would be wrong. The sequences are read from the
    *predicted complexes* (not Target/*.pdb) because the reprediction backends
    renumber chain A to 1..115 in both species, and epitope residue ids
    downstream are in that model numbering.

    ref_complex: {species: path to one representative predicted complex}"""
    from Bio import Align
    seqs = {}
    for sp in SPECIES:
        cx = parse_complex(ref_complex[sp])
        ch = cx["A"]
        resids = sorted(ch["seq"])
        seqs[sp] = (resids, "".join(ch["seq"][r] for r in resids))
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    try:
        from Bio.Align import substitution_matrices
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:
        aligner.match_score, aligner.mismatch_score = 1, -1
    aln = aligner.align(seqs["human"][1], seqs["mouse"][1])[0]
    h2m = {}
    for (hs, he), (ms, me) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(he - hs):
            h2m[seqs["human"][0][hs + k]] = seqs["mouse"][0][ms + k]
    ident = sum(seqs["human"][1][hs + k] == seqs["mouse"][1][ms + k]
                for (hs, he), (ms, me) in zip(aln.aligned[0], aln.aligned[1])
                for k in range(he - hs)) / max(1, len(h2m))
    print(f"  hPD1<->mPD1 alignment: {len(h2m)}/115 positions mapped, "
          f"{ident:.0%} identity")
    return h2m


# ------------------------------------------------------------------ main flow
def main():
    summary = pd.read_csv(SUMMARY_CSV)
    summary["name"] = summary["name"].astype(str).str.strip()
    names = summary["name"].tolist()

    # ---- 1. metric CSVs
    metrics = {}
    for meth, src in METHOD_SOURCES.items():
        for sp in SPECIES:
            p = os.path.join(src[sp], src["csv"])
            d = pd.read_csv(p)
            d["name"] = d["name"].astype(str).str.strip()
            metrics[(meth, sp)] = d
            print(f"  {meth:12s} {sp:5s}: {len(d)} rows from {os.path.relpath(p, ROOT)}")

    # ---- 2. structural features from the predicted complexes
    print("\nRecomputing structural features from predicted complexes ...")
    struct = {}      # (name, species, method) -> dict
    for meth, src in METHOD_SOURCES.items():
        for sp in SPECIES:
            mpnn = os.path.join(src[sp], "MPNN")
            n_ok = 0
            for nm in names:
                pdb = os.path.join(mpnn, f"{nm}_model1.pdb")
                if not os.path.isfile(pdb):
                    continue
                cx = parse_complex(pdb)
                feat = interface_features(cx)
                if feat is None:
                    continue
                feat["_cx"] = cx
                struct[(nm, sp, meth)] = feat
                n_ok += 1
            print(f"  {meth:12s} {sp:5s}: {n_ok}/{len(names)} complexes parsed")

    ref_complex = {
        sp: os.path.join(METHOD_SOURCES["ESMFold2"][sp], "MPNN", f"{names[0]}_model1.pdb")
        for sp in SPECIES}
    h2m = species_residue_map(ref_complex)

    # ---- 3. cross-method agreement (the "together" structural signal)
    print("\nComputing cross-method agreement ...")
    xmethod = {}
    for nm, sp in itertools.product(names, SPECIES):
        a = struct.get((nm, sp, "ESMFold2"))
        b = struct.get((nm, sp, "Protenix-v2"))
        if a is None or b is None:
            continue
        # superpose Protenix onto ESMFold2 using the shared target CAs, then
        # measure how far apart the two models place the binder.
        shared_t = sorted(set(a["_cx"]["A"]["ca"]) & set(b["_cx"]["A"]["ca"]))
        shared_b = sorted(set(a["_cx"]["B"]["ca"]) & set(b["_cx"]["B"]["ca"]))
        rec = {"xm_epitope_jaccard": jaccard(a["epitope"], b["epitope"]),
               "xm_paratope_jaccard": jaccard(a["paratope"], b["paratope"])}
        if len(shared_t) >= 3 and len(shared_b) >= 3:
            Pt = np.array([b["_cx"]["A"]["ca"][r] for r in shared_t])
            Qt = np.array([a["_cx"]["A"]["ca"][r] for r in shared_t])
            Pb = np.array([b["_cx"]["B"]["ca"][r] for r in shared_b])
            Qb = np.array([a["_cx"]["B"]["ca"][r] for r in shared_b])
            t_rmsd, Pb_moved = kabsch_rmsd(Pt, Qt, apply_to=Pb)
            rec["xm_target_rmsd"] = t_rmsd
            rec["xm_binder_rmsd"] = float(np.sqrt(((Pb_moved - Qb) ** 2).sum(1).mean()))
            # binder fold agreement alone (superposing on the binder itself),
            # separating "same fold, different pose" from "different fold"
            rec["xm_binder_fold_rmsd"] = kabsch_rmsd(Pb, Qb)[0]
        if a["epitope_com"] is not None and b["epitope_com"] is not None:
            # epitope-centre distance in each model's own target frame
            rec["xm_epitope_com_dist"] = float(np.linalg.norm(
                a["epitope_com"] - b["epitope_com"]))
        xmethod[(nm, sp)] = rec

    # ---- 4. cross-species epitope conservation, per method
    xspecies = {}
    for nm, meth in itertools.product(names, METHODS):
        a = struct.get((nm, "human", meth))
        b = struct.get((nm, "mouse", meth))
        if a is None or b is None:
            continue
        mapped = frozenset(h2m[r] for r in a["epitope"] if r in h2m)
        xspecies[(nm, meth)] = {"xs_epitope_jaccard": jaccard(mapped, b["epitope"])}

    # ---- 5. epitope overlap with the natural PD-L1 footprint
    native_ep = {}
    for sp, meth in itertools.product(SPECIES, METHODS):
        f = struct.get((POSITIVE_CONTROL, sp, meth))
        if f is not None:
            native_ep[(sp, meth)] = f["epitope"]

    # ---- 6. assemble
    STRUCT_COLS = ["n_contacts", "iface_nres_target", "iface_nres_binder",
                   "iface_pLDDT_binder", "iface_pLDDT_target",
                   "binder_pLDDT", "target_pLDDT"]
    METRIC_COLS = ["pLDDT", "pTM", "i_pTM", "pAE", "i_pAE"]

    long_rows, struct_rows = [], []
    for nm, sp, meth in itertools.product(names, SPECIES, METHODS):
        d = metrics[(meth, sp)]
        row = d.loc[d["name"] == nm]
        rec = {"name": nm, "species": sp, "method": meth}
        if len(row):
            r = row.iloc[0]
            rec["length"] = r["length"]
            rec["ok"] = bool(r["ok"])
            for c in METRIC_COLS:
                rec[c] = r[c] if c in d.columns else np.nan
        f = struct.get((nm, sp, meth))
        if f is not None:
            for c in STRUCT_COLS:
                rec[c] = f[c]
            ne = native_ep.get((sp, meth))
            rec["native_epitope_jaccard"] = jaccard(f["epitope"], ne) if ne else np.nan
            rec["epitope_resids"] = ";".join(str(x) for x in sorted(f["epitope"]))
        rec.update(xmethod.get((nm, sp), {}))
        rec.update(xspecies.get((nm, meth), {}))
        long_rows.append(rec)
        struct_rows.append(rec)

    long_df = pd.DataFrame(long_rows)
    long_df = long_df.merge(
        summary[["name", "VisterraName", "Classification",
                 "HuPD-1 Binding (1 uM)", "MuPD-1 Binding (1 uM)",
                 "Lysozyme Binding (1 uM)"]],
        on="name", how="left")
    is_ctrl = long_df["name"] == POSITIVE_CONTROL
    long_df[is_ctrl].to_csv(os.path.join(OUT, "positive_control.csv"), index=False)
    long_df = long_df[~is_ctrl].reset_index(drop=True)
    print(f"\nRemoved the {POSITIVE_CONTROL} positive control from the dataset "
          f"({int(is_ctrl.sum())} rows -> Analysis/data/positive_control.csv)")
    long_df.to_csv(os.path.join(OUT, "master_long.csv"), index=False)
    print(f"Wrote {os.path.join('Analysis/data', 'master_long.csv')}  "
          f"({len(long_df)} rows, {long_df['name'].nunique()} designs)")

    # ---- 7. wide, per-species tables (plot.py compatible column names)
    #  ESMFold2_ipTM / Protenix-v2_ipTM etc.
    for sp in SPECIES:
        wide = summary.loc[summary["name"] != POSITIVE_CONTROL,
                           ["name", "VisterraName", "Classification",
                            "HuPD-1 Binding (1 uM)", "MuPD-1 Binding (1 uM)",
                            "Lysozyme Binding (1 uM)"]].copy()
        sub = long_df[long_df["species"] == sp]
        wide = wide.merge(
            sub[sub["method"] == METHODS[0]][["name", "length"]], on="name", how="left")
        for meth in METHODS:
            m = sub[sub["method"] == meth].set_index("name")
            for c in METRIC_COLS + STRUCT_COLS + ["native_epitope_jaccard",
                                                  "xs_epitope_jaccard"]:
                if c not in m.columns:
                    continue
                alias = "ipTM" if c == "i_pTM" else c
                wide[f"{meth}_{alias}"] = wide["name"].map(m[c])
        # cross-method columns are method-agnostic: attach once
        xm = sub[sub["method"] == METHODS[0]].set_index("name")
        for c in [c for c in xm.columns if c.startswith("xm_")]:
            wide[c] = wide["name"].map(xm[c])
        p = os.path.join(OUT, f"merged_{sp}.csv")
        wide.to_csv(p, index=False)
        print(f"Wrote Analysis/data/merged_{sp}.csv  ({len(wide)} rows, "
              f"{len(wide.columns)} cols)")


if __name__ == "__main__":
    main()
