# Project: Binder/non-binder discrimination — cofolding confidence vs funnel metadynamics

## Goal

Compare how well deep-learning structure predictors discriminate binders from
non-binders, against an in-house funnel metadynamics (FMD) pipeline.

**Framing correction (important):** AF3, Protenix-v2, OpenDDE etc. are structure
predictors, not classifiers. What is actually being benchmarked is their
**confidence metrics** used as a discriminator. State this explicitly in any
writeup.

## Model lineup

| Model | Status | Notes |
|---|---|---|
| AF3 | in | Training cutoff 2021-09-30 |
| Protenix-v2 | in | 464M params, AF3 reproduction, +9–13 pts DockQ>0.23 over v1 on antibody benchmarks, has ranking mode |
| OpenDDE | in, off-label | Authors explicitly disclaim affinity ranking / virtual screening. Strength is Ab–Ag interface geometry |
| ~~ESMFold2~~ | **dropped** | Not a real model name. ESMFold is single-chain only, no multimer mode, no MSA pairing. Linker tricks give meaningless interface scores |
| Boltz-2 | suggested add | Has a trained affinity head — most interesting comparator |
| Chai-1 / AF2-Multimer | optional | Baselines |

### Scores to extract (best → worst for this task)
1. ipSAE, pDockQ2 — interface-local
2. iPAE (mean/min PAE across interface) — often beats ipTM
3. ipTM — global, poorly calibrated
4. iCS

**Do not use pTM.** Dominated by the larger chain; says nothing about interface.

Fix a **seed budget** and apply identically across all models. Report both top-1
(model's own ranking) and best-of-N. Only top-1 is achievable prospectively.

## Target class

Protein–protein, antibody–antigen, with a nanobody (VHH) focus.
Negatives: experimental, from display/panning campaigns.

## Dataset survey (completed)

| Resource | Nanobody content | Binder/non-binder labels |
|---|---|---|
| SKEMPI 2.0 | 3 complexes, 8 mutations (EGFR nanobodies 4KRL/4KRO/4KRP, all one 2013 SPR study) | No — ΔΔG only |
| AB-Bind | **Zero.** 1 monobody (3K2M, FN3 scaffold, 7 muts) | No — ΔΔG only |
| ANDD v2 xlsx (aggregator) | 30,119 rows → 24,797 unique VHH | **Labels stripped in merge.** Only 1,722 rows have Kd |
| **AVIDa-hIL6** | 573,891 pairs | **Yes, binary** |
| **AVIDa-SARS-CoV-2** | 77,003 pairs | **Yes, binary** |

**Lesson:** aggregators drop fields that don't fit their schema. Always go to the
primary release.

Other options not yet pursued: AbAgym (335k mutations, 67 Ab–Ag DMS experiments),
Mason trastuzumab–HER2 DMS, Lim et al. CTLA-4/PD-1, Chinery et al. (524,346
trastuzumab variants), AbDesign, FLAb, AbBiBench.

## Chosen datasets

### AVIDa-hIL6 (primary)
- HF: `COGNANO/AVIDa-hIL6` — DOI 10.57967/hf/4241 — **CC BY-NC 4.0**
- Mirror: https://www.cognano.co.jp/datasets/avida-hil6
- Code: https://github.com/cognano/AVIDa-hIL6
- Paper: arXiv 2306.03329 (NeurIPS 2023 D&B)
- 573,891 pairs = 20,980 binders + 552,911 non-binders (3.7% positive)
- 31 antigens (IL-6 WT + 30 artificial point mutants), ≥250 binders each
- 1 alpaca ("wizzy", male) — no individual-based split available
- Files: `AVIDa-hIL6.csv`, `antigen_sequences.csv` (105 MB total)

### AVIDa-SARS-CoV-2 (secondary)
- HF: `COGNANO/AVIDa-SARS-CoV-2` — **CC BY-NC 4.0**
- Code: https://github.com/cognano/AVIDa-SARS-CoV-2
- Paper: arXiv 2405.18749
- 77,003 pairs = 22,002 binders + 55,001 non-binders (28.6% positive — much better balance)
- 13 antigens: WT, D614G, Alpha, Alpha+K417N, Alpha+E484K, Beta, Delta, Kappa,
  Lambda, Omicron BA.1, PMS, S2-domain, OC43
- 2 alpacas: P (10,487 unique binders), C (3,651), only 60 shared
  → authors' benchmark trains on P, tests on C
- **427 VHHs labelled binder to some antigens, non-binder to others** — the best
  hard-negative set available for nanobodies
- Baseline to beat: AntiBERTa2-CSSP, F1 0.652, AUPRC 0.690 on the P→C split

## CRITICAL preprocessing

AVIDa VHH sequences are **not** bare domains. Every one is wrapped:

```
MKYLLPTAAAGLLLLAAQPAMA + QVQLQESGGG...VTVSS + HHHHHH
```

- N-term: 22-residue phagemid signal peptide
- C-term: His6 tag

**Strip both before any structure prediction or MD.** The signal peptide becomes a
disordered 22-residue tail that wrecks ipTM and funnel geometry. Real domain runs
`QVQLQESGGG` → `...VTVSS`.

### HF loader bug
`load_dataset("COGNANO/AVIDa-hIL6")` throws `DatasetGenerationCastError` — the two
CSVs have mismatched schemas and aren't declared as separate configs. Use:

```python
from huggingface_hub import hf_hub_download
import pandas as pd
pairs = hf_hub_download("COGNANO/AVIDa-hIL6", "AVIDa-hIL6.csv", repo_type="dataset")
ags   = hf_hub_download("COGNANO/AVIDa-hIL6", "antigen_sequences.csv", repo_type="dataset")
df = pd.read_csv(pairs).merge(pd.read_csv(ags), on="Ag_label", how="left")
```

## Label semantics

AVIDa labels are **proportion-shift calls, not affinities**. Binomial test on NGS
read counts before/after panning, p ≤ 0.05. Binder if proportion rose, non-binder
if it fell. ~97% of pairs were non-significant and discarded — these are the
confidently-separable extremes, not a random sample. No Kd anywhere.

Validation: 20 AVIDa-hIL6 VHHs confirmed by immunofluorescence + BLI; 9
AVIDa-SARS-CoV-2 VHHs confirmed binding spike variants in a separate paper.

## Benchmark design decisions

### Negative strata — report separately, never pooled
1. **Within-series** (library variants that lost binding) — hard; where FMD should win
2. **Cross-target / family-matched** (same VHH vs different antigen) — where
   interface confidence scores have a real shot
3. Random swapped — easy; include only as a sanity floor

### Swapped negatives: do NOT use for the FMD arm
FMD needs a bound starting pose. A swapped pair has none, so you'd have to dock
it — either with the models under test (circular) or with a third tool (unknown
quality on non-binders). The asymmetry correlates perfectly with the class label
and inflates FMD's AUC. Mutational negatives don't have this problem.

→ Use SAbDab-nano / swapped pairs for **cofolding only** (cheap, well-powered).
→ Run FMD only where every member has a legitimate structure.

### Leakage control
AF3 and Protenix cutoff 2021-09-30. FMD has no training set. Benchmarking on
pre-cutoff PDB entries hands the DL models a memory advantage. Use post-cutoff or
de novo systems; cluster at ≤40% identity on both VHH and antigen; report strata
separately. **Reviewers will go straight to this.**

### Size confound
ipTM family is sensitive to chain length and length ratio. Size-match negatives
to within ~20% of the antigen they replace. Run a permutation control: shuffle
labels within size bins, confirm AUC → 0.5.

### Metrics
- ROC-AUC **and** PR-AUC (PR matters at realistic class imbalance — hIL6 is 26:1)
- Enrichment factor at 1% and 5%
- DeLong's test, paired (all methods see the same systems)
- Bootstrap CIs **over campaigns/systems, not over variants** — variants within one
  campaign are not independent
- n≈30 gives AUC CI of roughly ±0.15. Powered for large differences only.
  Don't claim 0.78 > 0.72 at that n.

### Decision rule for FMD on non-binders
True non-binders have no well-defined bound minimum, so ΔG isn't defined.
**State the rule in advance**: depth of deepest minimum, or ΔG threshold with
unconverged runs called non-binder.

### Threshold
Fix the KD/label threshold before running anything — it silently sets both class
balance and how much free energy FMD must resolve.

## FMD technical risks (protein–protein, Ab–Ag)

FMD was designed for protein–ligand. Ab–Ag buries 700–900 Å² per side.
1. **Funnel geometry** — must contain the binder's rotational freedom near the
   bound state; large volume correction; sensitive to cone placement
2. **CVs** — projection distance alone insufficient. Fv can rotate with little COM
   change. Need orientational or path CVs. This choice will dominate ΔG
3. **Convergence** — profiles look converged long before they are. Block analysis,
   profile vs simulation time, ≥2 replicas per system with different seeds.
   Replicas disagreeing >1.5 kcal/mol → report, don't average away
4. **Plan B** if it won't converge: umbrella sampling / PMF along a pulling
   coordinate, or restraint-based ABFE

## Tiered compute plan

- Cofolding scores on the **full** dataset (cheap, well-powered DL-arm result)
- **Stratified subsample** for FMD: balanced across negative strata, spanning the
  range, including cases where DL scores confidently fail
- Head-to-head reported only on the subsample, full-set DL numbers as context
- **Select the subsample before running FMD and state the rule.** Picking FMD
  systems after seeing where cofolding fails gets the comparison dismissed

## SKEMPI panel (if conventional PPI systems are still wanted)

Key tension found in the data: systems with the most mutations (ovomucoid, BLIP,
barnase–barstar) have fM–pM wild types, so only 2–8% of mutants cross 1 µM.
Systems where mutations genuinely abolish binding have µM wild types and less data.
**Cannot get both from one system.**

| Complex | Target / binder | ~Size | Single muts | >1 µM | Role |
|---|---|---|---|---|---|
| 1BRS | Barnase / barstar | 110+89 | 49 | 0 | Calibration — most published PMF/metadynamics precedent. Run first |
| 4G0N | H-Ras / Raf-RBD | 166+80 | 48 | 18 | Best single discrimination system. Needs GppNHp + Mg²⁺ params |
| 1EMV + 2WPT | E9 DNase / Im9, Im2 | 134+86 | 46+45 | 0+8 | Cognate vs non-cognate, 16 fM vs 15 nM. Check Zn²⁺ in HNH motif |
| 1PPF | HLE / OMTKY3 | 218+56 | 203 | 16 | Statistics system. OMTKY3 is 56 res, disulfide-stapled, rigid. Mutations concentrated at P1 |
| 1JTG + 2G2U | TEM-1, SHV-1 / BLIP | 263+165 | 138+29 | 7+10 | Same binder, 2000-fold specificity. Largest — drop first if compute-bound |

Ras-family bonus: 1C1Y (Rap1a–Raf-RBD), 1K8R (Ras–Byr2), 1LFD (Ras–RalGDS) give
free cross-target negatives with compact RA/RBD binders.

Smaller swaps: 1FFW (CheY–CheA, ~190 res), 1XD3 (UCH-L3–ubiquitin).

**Avoid:** 1A4Y / 1Z7X (ribonuclease inhibitor, 460-res LRR horseshoe);
3BT1 / 1GC1 / 1A22 (ΔΔG range compressed, no discrimination signal).

Note: 1JTG, 1KTZ, 1FFW, 1AK4 appear in **both** SKEMPI and AB-Bind — dedupe if merging.

## Open questions

- Is the user's own display campaign VHH or conventional Fv/Fab?
- Provenance of the ANDD v2 xlsx — can the stripped AVIDa labels be recovered?
- Per-antigen binder counts in AVIDa-hIL6, and how many VHHs flip label across the
  31 IL-6 variants (the IL-6 analogue of the SARS-CoV-2 427-sequence set)
