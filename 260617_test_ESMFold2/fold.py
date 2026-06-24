import sys
import os
import time
import csv

import torch
import numpy as np

from esm.models.esmfold2 import (
    ESMFold2InputBuilder,
    ProteinInput,
    StructurePredictionInput,
)
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model


def compute_iptm_pair_mean(result, chain_ids_unique):
    """Mean iPTM between TCR chains (indices 1,2) and pMHC chains (indices 0,4)."""
    if result.pair_chains_iptm is None:
        return None
    pc = result.pair_chains_iptm.numpy()
    i_tcr = [1, 2]  # TCRA, TCRB
    i_pmhc = [0, 4]  # MHC, PEPTIDE
    vals = []
    for i in i_tcr:
        for j in i_pmhc:
            if i < pc.shape[0] and j < pc.shape[1]:
                vals.append(pc[i, j])
                vals.append(pc[j, i])
    return float(np.mean(vals)) if vals else None


def fold_one(complex_id, seqs, model, builder):
    """Fold a single complex and return (result, mc, t_fold)."""
    spi = StructurePredictionInput(
        sequences=[
            ProteinInput(id="MHC", sequence=seqs["MHC"]),
            ProteinInput(id="TCRA", sequence=seqs["TCRA"]),
            ProteinInput(id="TCRB", sequence=seqs["TCRB"]),
            ProteinInput(id="B2M", sequence=seqs["B2M"]),
            ProteinInput(id="PEPTIDE", sequence=seqs["PEPTIDE"]),
        ]
    )

    t_start = time.time()
    result = builder.fold(
        model, spi, num_loops=3, num_sampling_steps=50,
        num_diffusion_samples=1, seed=0,
    )
    t_fold = time.time() - t_start
    return result, t_fold


def write_per_complex_csv(path, complex_id, result, t_fold):
    """Write the detailed per-complex metrics CSV."""
    mc = result.complex
    chain_lookup = mc.metadata.chain_lookup
    chain_ids_unique = sorted(set(mc.chain_id.tolist()))
    label_map = {0: "MHC", 1: "TCRA", 2: "TCRB", 3: "B2M", 4: "PEPTIDE"}

    with open(path, "w", newline="") as f:
        w = csv.writer(f)

        w.writerow(["metric", "value"])
        w.writerow(["id", complex_id])
        w.writerow(["fold_time_seconds", f"{t_fold:.2f}"])
        w.writerow(["ptm", f"{float(result.ptm):.4f}"])
        w.writerow(["iptm", f"{float(result.iptm):.4f}"])
        w.writerow(["plddt_mean", f"{float(result.plddt.mean()):.4f}"])
        w.writerow(["plddt_median", f"{float(result.plddt.median()):.4f}"])
        w.writerow(["plddt_min", f"{float(result.plddt.min()):.4f}"])
        w.writerow(["plddt_max", f"{float(result.plddt.max()):.4f}"])
        if result.pae is not None:
            w.writerow(["pae_mean", f"{result.pae.numpy().mean():.4f}"])
        iptm_pair = compute_iptm_pair_mean(result, chain_ids_unique)
        if iptm_pair is not None:
            w.writerow(["iptm_pair_mean", f"{iptm_pair:.4f}"])
        w.writerow(["total_tokens", len(mc.sequence)])
        w.writerow(["total_atoms", len(mc.atom_positions)])
        w.writerow([])

        w.writerow(["chain_id", "chain_label", "n_tokens", "plddt_mean", "plddt_min", "plddt_max"])
        for idx, cid in enumerate(chain_ids_unique):
            mask = mc.chain_id == cid
            plddt_chain = mc.plddt[mask]
            name = chain_lookup.get(cid, str(cid))
            label = label_map.get(idx, name)
            w.writerow([
                name, label, int(mask.sum()),
                f"{float(plddt_chain.mean()):.4f}",
                f"{float(plddt_chain.min()):.4f}",
                f"{float(plddt_chain.max()):.4f}",
            ])
        w.writerow([])

        if result.pair_chains_iptm is not None:
            pc_iptm = result.pair_chains_iptm.numpy()
            n = pc_iptm.shape[0]
            chain_names = {
                cid: chain_lookup.get(cid, str(cid)) for cid in chain_ids_unique
            }
            w.writerow(["pair_chains_iptm"])
            header_row = [""] + [
                label_map.get(i, chain_names.get(chain_ids_unique[i], str(i)))
                for i in range(n)
            ]
            w.writerow(header_row)
            for i in range(n):
                row_label = label_map.get(i, chain_names.get(chain_ids_unique[i], str(i)))
                row = [row_label] + [f"{pc_iptm[i, j]:.4f}" for j in range(n)]
                w.writerow(row)
        w.writerow([])

        w.writerow(["token_index", "chain_id", "chain_label", "residue_name", "plddt"])
        for tok_idx in range(len(mc.sequence)):
            cid = int(mc.chain_id[tok_idx])
            name = chain_lookup.get(cid, str(cid))
            chain_pos = list(chain_ids_unique).index(cid)
            label = label_map.get(chain_pos, name)
            w.writerow([
                tok_idx, name, label, mc.sequence[tok_idx],
                f"{mc.plddt[tok_idx]:.4f}",
            ])


# ── Main ─────────────────────────────────────────────────────────────────────
input_csv = sys.argv[1] if len(sys.argv) > 1 else "input.csv"

with open(input_csv) as f:
    reader = csv.DictReader(f)
    rows = list(reader)

print(f"Read {len(rows)} complexes from {input_csv}")

os.makedirs("models", exist_ok=True)
os.makedirs("metrics", exist_ok=True)

print("Loading model...")
model = ESMFold2Model.from_pretrained("biohub/ESMFold2-Fast").cuda().eval()
builder = ESMFold2InputBuilder()
print("Model loaded.")

summary_rows = []

for i, row in enumerate(rows):
    complex_id = row["ID"]
    print(f"\n[{i+1}/{len(rows)}] Folding {complex_id}...")

    seqs = {k: row[k] for k in ("TCRA", "TCRB", "MHC", "B2M", "PEPTIDE")}
    result, t_fold = fold_one(complex_id, seqs, model, builder)

    mc = result.complex
    chain_ids_unique = sorted(set(mc.chain_id.tolist()))

    iptm_pair = compute_iptm_pair_mean(result, chain_ids_unique)

    print(f"  pTM={float(result.ptm):.3f}  ipTM={float(result.iptm):.3f}"
          f"  pLDDT={float(result.plddt.mean()):.3f}"
          f"  iptm_pair_mean={iptm_pair:.3f}" if iptm_pair else ""
          f"  time={t_fold:.1f}s")

    cif_path = f"models/{complex_id}.cif"
    with open(cif_path, "w") as f:
        f.write(mc.to_mmcif())

    write_per_complex_csv(f"metrics/{complex_id}.csv", complex_id, result, t_fold)

    summary_rows.append({
        "ID": complex_id,
        "iptm": f"{float(result.iptm):.4f}",
        "ptm": f"{float(result.ptm):.4f}",
        "plddt_mean": f"{float(result.plddt.mean()):.4f}",
        "fold_time_seconds": f"{t_fold:.2f}",
        "iptm_pair_mean": f"{iptm_pair:.4f}" if iptm_pair is not None else "",
    })

    print(f"  Saved models/{complex_id}.cif and metrics/{complex_id}.csv")

out_csv = "output.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ID", "iptm", "ptm", "plddt_mean",
                                       "fold_time_seconds", "iptm_pair_mean"])
    w.writeheader()
    w.writerows(summary_rows)

print(f"\nDone. Summary written to {out_csv}")
