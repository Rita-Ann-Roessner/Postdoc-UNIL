#!/usr/bin/env python3
"""Build the self-contained HTML report for the PD-1 reprediction benchmark."""
import base64, os, mimetypes
import pandas as pd

ROOT = "/Users/admin/Downloads/Reprediction_PD1_designs"
FIG = os.path.join(ROOT, "Analysis", "figures")
OUT = os.path.join(ROOT, "Analysis", "report.html")

def img(name):
    p = os.path.join(FIG, name)
    b = base64.b64encode(open(p, "rb").read()).decode()
    return f"data:image/png;base64,{b}"

auc = pd.read_csv(os.path.join(ROOT, "Analysis/data/auc_table.csv"))

def row(pred):
    out = {}
    for t in ["human", "mouse", "cross_reactive"]:
        r = auc[(auc.task == t) & (auc.predictor == pred)]
        out[t] = None if r.empty else (float(r.AUC.iloc[0]), float(r.AUC_lo.iloc[0]),
                                       float(r.AUC_hi.iloc[0]), r.estimate.iloc[0])
    return out

SCOREBOARD = [
    ("ESMFold2 interface pLDDT", "esm", "recomputed from PDB B-factors"),
    ("ESMFold2 ipTM", "esm", "pipeline metric"),
    ("ESMFold2 composite (CV)", "esm", "5 ESMFold2 metrics, out-of-fold"),
    ("Protenix-v2 ipTM", "prx", "pipeline metric"),
    ("Protenix-v2 interface pLDDT", "prx", "recomputed from PDB B-factors"),
    ("Protenix-v2 composite (CV)", "prx", "6 Protenix metrics, out-of-fold"),
    ("mean-z(ipTM)", "both", "unweighted z-average, no fitting"),
    ("max-z(ipTM)", "both", "either method may be high"),
    ("Both ipTM (CV)", "both", "logistic on both ipTMs, out-of-fold"),
    ("All metrics, both methods (CV)", "both", "11 features, out-of-fold"),
    ("cross-method pose RMSD", "cons", "do the two models agree on the pose?"),
    ("cross-method epitope overlap", "cons", "do they agree on the epitope?"),
    ("binder length (baseline)", "base", "confound check"),
]

def bar(v, group):
    """AUC rendered on a 0.50-1.00 scale; 0.50 is the chance line."""
    if v is None:
        return '<td class="num"><span class="na">--</span></td>'
    a, lo, hi, est = v
    w = max(0.0, min(1.0, (a - 0.5) / 0.5)) * 100
    llo = max(0.0, min(1.0, (lo - 0.5) / 0.5)) * 100
    lhi = max(0.0, min(1.0, (hi - 0.5) / 0.5)) * 100
    mark = "&#9670;" if est == "out-of-fold" else ""
    return (f'<td class="num"><div class="cell">'
            f'<span class="v">{a:.3f}{mark}</span>'
            f'<span class="ci">{lo:.2f}&ndash;{hi:.2f}</span>'
            f'<span class="track"><i class="ciband" style="left:{llo:.1f}%;'
            f'width:{max(lhi-llo,1):.1f}%"></i>'
            f'<i class="fill g-{group}" style="width:{w:.1f}%"></i></span>'
            f'</div></td>')

rows_html = []
best = {t: max(v[t][0] for _, _, _ in [(0,0,0)] for v in [row(p)] for p, _, _ in [] ) if False else None for t in []}
# best per task across the scoreboard, for the winner underline
tasks = ["human", "mouse", "cross_reactive"]
vals = {t: [] for t in tasks}
cache = {}
for pred, grp, note in SCOREBOARD:
    cache[pred] = row(pred)
    for t in tasks:
        if cache[pred][t]:
            vals[t].append(cache[pred][t][0])
bestv = {t: max(vals[t]) for t in tasks}

for pred, grp, note in SCOREBOARD:
    r = cache[pred]
    cells = ""
    for t in tasks:
        c = bar(r[t], grp)
        if r[t] and abs(r[t][0] - bestv[t]) < 1e-9:
            c = c.replace('<td class="num">', '<td class="num best">')
        cells += c
    rows_html.append(
        f'<tr><th scope="row"><span class="dot g-{grp}"></span>'
        f'<span class="pname">{pred}</span>'
        f'<span class="pnote">{note}</span></th>{cells}</tr>')
rows_html = "\n".join(rows_html)

FIGURES = [
    ("fig3_auc_forest.png", "Every score, every task",
     "AUC with 95% bootstrap CI for all 22 candidate predictors across the three tasks. "
     "Circles are apparent AUC; diamonds are out-of-fold estimates from 20&times;5-fold CV. "
     "The dotted line marks the best score in each panel. Nothing in the green or purple "
     "bands &mdash; the ways of using both methods together &mdash; clears the best blue bar."),
    ("fig1_score_distributions.png", "What the separation actually looks like",
     "Binder vs non-binder for the four most informative scores. The distributions overlap "
     "heavily: this is a triage signal, not a threshold. p-values are one-sided Mann-Whitney U."),
    ("fig2_roc.png", "ROC, alone and fused",
     "The shaded strip is the low-false-positive region that AUC<sub>0.1</sub> integrates &mdash; "
     "the part of the curve a real screening campaign lives in."),
    ("fig4_method_agreement.png", "The two methods rank different designs",
     "ESMFold2 vs Protenix-v2 ipTM, one point per design. Only 5 of the 16 hPD-1 binders land "
     "in the top-right quadrant where both methods are above their median. The raw material for "
     "an ensemble exists; this dataset is too small to exploit it."),
    ("fig6_budget_matched.png", "Spending a fixed budget",
     "For a fixed number of designs you can afford to express and test, which selection "
     "strategy finds the most real binders? Differences between one method, fused scores and a "
     "split budget are one or two binders &mdash; noise. The AND gate is reliably worst."),
    ("fig5_screening_yield.png", "Yield as a function of how much you order",
     "Recall (top) and enrichment over random (bottom) as the ordered fraction of the library "
     "grows. Enrichment settles around 2&ndash;3&times; and decays as the budget widens."),
]
figs_html = "\n".join(
    f'<figure class="plate">\n  <img src="{img(f)}" alt="{t}">\n'
    f'  <figcaption><span class="figtitle">{t}</span>{c}</figcaption>\n</figure>'
    for f, t, c in FIGURES)

FINDINGS = [
    ("Do the reprediction scores predict binding at all?", "yes, weakly",
     "ok",
     "Every ipTM-type score beats chance (p = 0.0006&ndash;0.06). But enrichment at a realistic "
     "budget is only 2&ndash;3&times; over random: ESMFold2's top 10 hPD-1 picks contain 3 of the "
     "16 true binders; its top 25 contain 9. Useful for triage, not for a go/no-go call."),
    ("Is one method better than the other?", "yes, on human PD-1",
     "ok",
     "ESMFold2 wins hPD-1 by &Delta;AUC 0.118 (DeLong p = 0.060) and cross-reactivity by 0.158 "
     "(p = 0.022). On mPD-1 they are indistinguishable (&Delta;AUC 0.003, p = 0.97). Protenix-v2's "
     "ipTM is compressed into a narrow 0.84&ndash;0.93 band, and that is what costs it discrimination."),
    ("Is there a better score than ipTM already in the outputs?", "yes &mdash; interface pLDDT",
     "ok",
     "Averaging per-residue pLDDT over interface residues only is the single best hPD-1 predictor "
     "in the benchmark (AUC 0.787 vs 0.758 for ipTM). It is not in the metrics table; it comes out "
     "of the B-factor column of PDBs you already have."),
    ("Does using both methods together beat using the better one alone?", "no, by any route",
     "null",
     "Score fusion (mean-z, min-z, max-z, rank-mean) lands inside the CI of the better single "
     "method. Learned fusion never beats ESMFold2 alone out-of-fold, and throwing all 11 metrics "
     "in makes it worse (hPD-1 0.710 vs 0.757, p = 0.007) &mdash; 11 features on 16 positives. "
     "Structural consensus &mdash; whether the two models agree on the pose (median 2.6 &Aring; RMSD) "
     "or the epitope (median Jaccard 0.83) &mdash; carries no signal at all (AUC 0.46&ndash;0.61). "
     "Two models agreeing is not evidence that a design binds."),
]
findings_html = "\n".join(
    f'<section class="finding">\n  <div class="qhead">'
    f'<span class="eyebrow">question</span>'
    f'<h3>{q}</h3>'
    f'<span class="chip c-{cls}">{verdict}</span></div>\n'
    f'  <p>{body}</p>\n</section>' for q, verdict, cls, body in FINDINGS)

html = f"""<title>PD-1 Reprediction Benchmark</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,500;0,6..72,600;1,6..72,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500;600&display=swap">
<style>
:root {{
  --paper:#F2F5F8; --surface:#FFFFFF; --surface-2:#E9EEF3;
  --ink:#141A21; --ink-2:#3C4854; --muted:#5F6B78; --rule:#D5DDE5;
  --esm:#2F5E96; --prx:#B85C34; --both:#2C6E4B; --cons:#6A5A96; --base:#7C858F;
  --ok:#2C6E4B; --null:#7C858F;
  --chip-ok-bg:#DCE9E1; --chip-null-bg:#E4E7EA;
  --shadow:0 1px 2px rgba(20,26,33,.06), 0 8px 24px -16px rgba(20,26,33,.18);
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper:#0D1218; --surface:#141B23; --surface-2:#1B242E;
    --ink:#E3E9EF; --ink-2:#B6C1CC; --muted:#8E99A5; --rule:#26313D;
    --esm:#7CA6D9; --prx:#E09168; --both:#61B287; --cons:#A296CE; --base:#8B949E;
    --ok:#61B287; --null:#8B949E;
    --chip-ok-bg:#1B2D24; --chip-null-bg:#1E252C;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -18px rgba(0,0,0,.8);
  }}
}}
:root[data-theme="dark"] {{
  --paper:#0D1218; --surface:#141B23; --surface-2:#1B242E;
  --ink:#E3E9EF; --ink-2:#B6C1CC; --muted:#8E99A5; --rule:#26313D;
  --esm:#7CA6D9; --prx:#E09168; --both:#61B287; --cons:#A296CE; --base:#8B949E;
  --ok:#61B287; --null:#8B949E;
  --chip-ok-bg:#1B2D24; --chip-null-bg:#1E252C;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -18px rgba(0,0,0,.8);
}}

* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:"IBM Plex Sans", ui-sans-serif, system-ui, sans-serif;
  font-size:16.5px; line-height:1.62; -webkit-font-smoothing:antialiased;
}}
.wrap {{ max-width:1120px; margin:0 auto; padding:0 clamp(18px,4vw,44px) 96px; }}
.prose {{ max-width:70ch; }}
h1,h2,h3 {{ font-family:Newsreader, Georgia, serif; font-weight:600; text-wrap:balance; margin:0; }}
h1 {{ font-size:clamp(2.1rem,4.6vw,3.3rem); line-height:1.08; letter-spacing:-.015em; }}
h2 {{ font-size:clamp(1.45rem,2.4vw,1.9rem); line-height:1.18; letter-spacing:-.01em; }}
h3 {{ font-size:1.2rem; line-height:1.28; font-weight:500; }}
p {{ margin:0; }}
a {{ color:var(--esm); }}
.eyebrow {{
  font-family:"IBM Plex Mono", ui-monospace, monospace; font-size:.665rem;
  letter-spacing:.14em; text-transform:uppercase; color:var(--muted); font-weight:500;
}}
code, .mono {{ font-family:"IBM Plex Mono", ui-monospace, monospace; font-size:.86em; }}

/* ---- masthead ---- */
header.top {{ padding:clamp(48px,7vw,84px) 0 0; }}
header.top .kicker {{ display:flex; gap:14px; align-items:baseline; flex-wrap:wrap; margin-bottom:20px; }}
header.top .rule {{ height:2px; background:var(--ink); width:52px; }}
header.top .dek {{
  font-family:Newsreader, Georgia, serif; font-size:clamp(1.05rem,1.9vw,1.3rem);
  color:var(--ink-2); max-width:60ch; margin-top:18px; line-height:1.45;
}}
.meta {{
  display:grid; grid-template-columns:repeat(auto-fit,minmax(158px,1fr)); gap:1px;
  background:var(--rule); border:1px solid var(--rule); border-radius:3px;
  margin:34px 0 0; overflow:hidden;
}}
.meta div {{ background:var(--surface); padding:12px 14px; }}
.meta dt {{ font-family:"IBM Plex Mono",monospace; font-size:.64rem; letter-spacing:.12em;
  text-transform:uppercase; color:var(--muted); margin:0 0 4px; }}
.meta dd {{ margin:0; font-size:.94rem; font-weight:500; font-variant-numeric:tabular-nums; }}

section.block {{ margin-top:clamp(52px,7vw,84px); }}
.blockhead {{ display:flex; align-items:baseline; gap:14px; margin-bottom:8px; }}
.blockhead .rule {{ flex:1; height:1px; background:var(--rule); }}
.lede {{ color:var(--ink-2); margin-top:14px; max-width:70ch; }}

/* ---- findings ---- */
.findings {{ display:flex; flex-direction:column; gap:1px; background:var(--rule);
  border:1px solid var(--rule); border-radius:3px; margin-top:26px; overflow:hidden; }}
.finding {{ background:var(--surface); padding:22px clamp(16px,2.4vw,26px); }}
.qhead {{ display:grid; grid-template-columns:1fr auto; gap:6px 16px;
  align-items:start; margin-bottom:10px; }}
.qhead .eyebrow {{ grid-column:1; }}
.qhead h3 {{ grid-column:1; }}
.chip {{ grid-column:2; grid-row:1 / span 2; align-self:center;
  font-family:"IBM Plex Mono",monospace; font-size:.7rem; letter-spacing:.06em;
  padding:5px 10px; border-radius:2px; white-space:nowrap; font-weight:500; }}
.c-ok {{ background:var(--chip-ok-bg); color:var(--ok); }}
.c-null {{ background:var(--chip-null-bg); color:var(--null); }}
.finding p {{ color:var(--ink-2); max-width:78ch; }}

/* ---- scoreboard ---- */
.tablewrap {{ overflow-x:auto; margin-top:26px; border:1px solid var(--rule);
  border-radius:3px; background:var(--surface); box-shadow:var(--shadow); }}
table.score {{ border-collapse:collapse; width:100%; min-width:760px; }}
table.score thead th {{
  font-family:"IBM Plex Mono",monospace; font-size:.66rem; letter-spacing:.1em;
  text-transform:uppercase; color:var(--muted); font-weight:500; text-align:left;
  padding:14px 16px 10px; border-bottom:1px solid var(--rule); vertical-align:bottom;
}}
table.score thead th.tnum {{ text-align:right; width:20%; }}
table.score thead th small {{ display:block; letter-spacing:.04em; text-transform:none;
  font-size:.68rem; color:var(--muted); font-weight:400; margin-top:3px; }}
table.score th[scope=row] {{ text-align:left; font-weight:400; padding:11px 16px;
  border-bottom:1px solid var(--rule); }}
.pname {{ display:block; font-size:.94rem; color:var(--ink); }}
.pnote {{ display:block; font-size:.76rem; color:var(--muted); }}
.dot {{ display:inline-block; width:7px; height:7px; border-radius:50%;
  margin-right:8px; vertical-align:1px; }}
.g-esm {{ background:var(--esm); }} .g-prx {{ background:var(--prx); }}
.g-both {{ background:var(--both); }} .g-cons {{ background:var(--cons); }}
.g-base {{ background:var(--base); }}
td.num {{ padding:11px 16px; border-bottom:1px solid var(--rule);
  border-left:1px solid var(--rule); vertical-align:middle; }}
.cell {{ display:grid; grid-template-columns:auto auto; gap:0 8px;
  justify-content:end; align-items:baseline; }}
.cell .v {{ font-family:"IBM Plex Mono",monospace; font-variant-numeric:tabular-nums;
  font-size:.98rem; font-weight:500; }}
.cell .ci {{ font-family:"IBM Plex Mono",monospace; font-size:.72rem; color:var(--muted);
  font-variant-numeric:tabular-nums; }}
.cell .track {{ grid-column:1 / -1; position:relative; height:5px; margin-top:6px;
  background:var(--surface-2); border-radius:1px; overflow:hidden; }}
.cell .track i {{ position:absolute; inset-block:0; display:block; }}
.cell .fill {{ left:0; border-radius:1px; opacity:.9; }}
.cell .ciband {{ background:var(--rule); }}
td.best .cell .v {{ text-decoration:underline; text-decoration-thickness:2px;
  text-underline-offset:3px; }}
.na {{ color:var(--muted); }}
.legend {{ display:flex; flex-wrap:wrap; gap:8px 22px; margin-top:14px;
  font-size:.78rem; color:var(--muted); }}
.legend span b {{ color:var(--ink-2); font-weight:500; }}

/* ---- figures ---- */
.plate {{ margin:34px 0 0; padding:0; background:var(--surface);
  border:1px solid var(--rule); border-radius:3px; overflow:hidden;
  box-shadow:var(--shadow); }}
.plate img {{ display:block; width:100%; height:auto; background:#fff; }}
.plate figcaption {{ padding:16px clamp(16px,2.4vw,24px) 20px; font-size:.9rem;
  color:var(--ink-2); border-top:1px solid var(--rule); }}
.figtitle {{ display:block; font-family:Newsreader,Georgia,serif; font-size:1.12rem;
  font-weight:600; color:var(--ink); margin-bottom:5px; }}

/* ---- prose bits ---- */
ul.tight {{ margin:16px 0 0; padding-left:1.15em; }}
ul.tight li {{ margin-bottom:9px; color:var(--ink-2); }}
ul.tight li strong {{ color:var(--ink); font-weight:600; }}
.callout {{ margin-top:26px; padding:20px clamp(16px,2.4vw,24px);
  background:var(--surface); border:1px solid var(--rule);
  border-radius:3px; }}
.callout .eyebrow {{ display:block; margin-bottom:8px; }}
.callout p + p {{ margin-top:12px; }}
pre {{ margin:18px 0 0; padding:16px 18px; background:var(--surface);
  border:1px solid var(--rule); border-radius:3px; overflow-x:auto;
  font-family:"IBM Plex Mono",monospace; font-size:.8rem; line-height:1.65;
  color:var(--ink-2); }}
table.files {{ border-collapse:collapse; width:100%; min-width:600px; margin-top:0; }}
table.files td {{ padding:9px 16px; border-bottom:1px solid var(--rule);
  font-size:.88rem; vertical-align:top; }}
table.files td:first-child {{ font-family:"IBM Plex Mono",monospace; font-size:.79rem;
  white-space:nowrap; color:var(--esm); }}
table.files td:last-child {{ color:var(--ink-2); }}
footer {{ margin-top:clamp(56px,7vw,88px); padding-top:22px;
  border-top:1px solid var(--rule); color:var(--muted); font-size:.82rem; }}
:focus-visible {{ outline:2px solid var(--esm); outline-offset:2px; }}
@media (prefers-reduced-motion:reduce) {{ * {{ animation:none !important; transition:none !important; }} }}
</style>

<div class="wrap">

<header class="top">
  <div class="kicker">
    <span class="rule"></span>
    <span class="eyebrow">PD1_duo campaign &middot; reprediction benchmark</span>
  </div>
  <h1>Does reprediction tell you which designs bind?</h1>
  <p class="dek">ESMFold2 and Protenix-v2 scored against BLI across 95 de novo designs, alone
  and in every combination worth trying. Both carry real signal. Neither is decisive, and
  putting them together does not help.</p>
  <dl class="meta">
    <div><dt>designs</dt><dd>95</dd></div>
    <div><dt>complexes</dt><dd>380</dd></div>
    <div><dt>ground truth</dt><dd>BLI @ 1 &micro;M</dd></div>
    <div><dt>hPD-1 binders</dt><dd>16 / 94</dd></div>
    <div><dt>mPD-1 binders</dt><dd>12 / 93</dd></div>
    <div><dt>cross-reactive</dt><dd>9 / 95</dd></div>
  </dl>
</header>

<section class="block">
  <div class="blockhead"><h2>The four questions</h2><span class="rule"></span></div>
  <div class="findings">
{findings_html}
  </div>
</section>

<section class="block">
  <div class="blockhead"><h2>Scoreboard</h2><span class="rule"></span></div>
  <p class="lede">AUC with a stratified-bootstrap 95% CI. Bars run from 0.50 (chance) to 1.00;
  the pale segment behind each bar is the confidence interval. The best score in each column is
  underlined. &#9670; marks an out-of-fold estimate from 20&times;5-fold cross-validation
  rather than an apparent one.</p>
  <div class="tablewrap">
    <table class="score">
      <thead><tr>
        <th>score</th>
        <th class="tnum">hPD-1<small>16 binders / 94</small></th>
        <th class="tnum">mPD-1<small>12 binders / 93</small></th>
        <th class="tnum">cross-reactive<small>9 / 95, worst-of-both</small></th>
      </tr></thead>
      <tbody>
{rows_html}
      </tbody>
    </table>
  </div>
  <div class="legend">
    <span><span class="dot g-esm"></span><b>ESMFold2</b> alone</span>
    <span><span class="dot g-prx"></span><b>Protenix-v2</b> alone</span>
    <span><span class="dot g-both"></span><b>both</b>, fused or learned</span>
    <span><span class="dot g-cons"></span><b>both</b>, structural consensus</span>
    <span><span class="dot g-base"></span><b>baseline</b></span>
  </div>
</section>

<section class="block">
  <div class="blockhead"><h2>Why &ldquo;no benefit&rdquo; is not &ldquo;proven equivalent&rdquo;</h2><span class="rule"></span></div>
  <div class="prose">
    <p class="lede">The two methods are only moderately correlated &mdash; Spearman
    &rho; = 0.56 on human ipTM, 0.66 on mouse &mdash; and they flag <em>different</em> designs.
    At a top-25% cutoff on hPD-1, ESMFold2 catches 9 of the 16 true binders, Protenix-v2 catches
    6, and only 5 are caught by both. The raw material for a useful ensemble is there.</p>
  </div>
  <div class="callout">
    <span class="eyebrow">the limit of this dataset</span>
    <p>With 9&ndash;16 positives, the minimum &Delta;AUC detectable at 80% power is
    <span class="mono">0.10&ndash;0.21</span> depending on the task. <strong>Every fusion effect
    measured here is smaller than that.</strong> So the finding is &ldquo;no detectable
    benefit&rdquo;, not &ldquo;no benefit&rdquo;.</p>
    <p>What would settle it: more positives. More designs with BLI data, a harder negative set,
    or quantitative affinity instead of a binary call at a single concentration.</p>
  </div>
</section>

<section class="block">
  <div class="blockhead"><h2>Figures</h2><span class="rule"></span></div>
{figs_html}
</section>

<section class="block">
  <div class="blockhead"><h2>Sanity checks</h2><span class="rule"></span></div>
  <ul class="tight">
    <li><strong>Binder length is not driving anything.</strong> AUC 0.44&ndash;0.50 on all
    three tasks, so the signal is not a length artefact.</li>
    <li><strong>The natural-ligand control is out of the dataset, not just out of the
    scoring.</strong> HuPDL1-Fc is labelled <span class="mono">Cross-reactive</span> in the BLI
    sheet, but it is human PD-L1: 227 aa and a deposited complex.
    <span class="mono">build_dataset.py</span> drops it from every table and parks its row in
    <span class="mono">positive_control.csv</span>, so it cannot reach any analysis. Including it
    had inflated the original ESMFold2 figure from 0.758 to 0.77.</li>
    <li><strong>Protenix-v2 failed outright on that control</strong> &mdash; from the parked row:
    ipTM 0.29 / 0.30, pLDDT 0.55 / 0.63, and 50 &Aring; away from the ESMFold2 pose &mdash; while
    ESMFold2 handled it normally (ipTM 0.86 / 0.84). Out of the benchmark, but it was the only
    complex here whose answer is independently known.</li>
    <li><strong>Independent reimplementation agrees.</strong> The campaign's own
    <span class="mono">plot.py</span>, pointed at the merged tables, reproduces
    0.758 / 0.749 / 0.640 / 0.746 and cross-reactive 0.819 / 0.661 exactly.</li>
  </ul>
</section>

<section class="block">
  <div class="blockhead"><h2>What to do with this</h2><span class="rule"></span></div>
  <ul class="tight">
    <li><strong>Rank on ESMFold2.</strong> ipTM or interface pLDDT &mdash; statistically
    interchangeable, with interface pLDDT ahead on hPD-1.</li>
    <li><strong>Don't spend a Protenix-v2 run to re-rank the same designs.</strong> On this
    evidence it adds nothing, and on hPD-1 it dilutes.</li>
    <li><strong>Don't filter on cross-method agreement.</strong> Pose and epitope consensus are
    at chance.</li>
    <li><strong>Expect 2&ndash;3&times; enrichment.</strong> Ordering the top quarter of a library
    recovers a bit over half the binders. There is no threshold here that gives a clean
    pass/fail.</li>
    <li><strong>AF3 and OpenDDE plug straight in.</strong> Both directories are still empty; the
    pipeline picks them up the moment a <span class="mono">reprediction_metrics.csv</span> lands,
    with no code change.</li>
  </ul>
</section>

<section class="block">
  <div class="blockhead"><h2>Code and tables</h2><span class="rule"></span></div>
  <div class="tablewrap">
    <table class="files">
      <tr><td>Analysis/build_dataset.py</td><td>merges the four metric CSVs and recomputes
        structural features &mdash; interface pLDDT, contacts, epitopes, cross-method pose RMSD
        &mdash; from all 384 predicted complex PDBs</td></tr>
      <tr><td>Analysis/predictive_power.py</td><td>AUC, partial AUC, precision/enrichment at K,
        DeLong tests, cross-validated models, power</td></tr>
      <tr><td>Analysis/plot_predictive_power.py</td><td>figures 1&ndash;6</td></tr>
      <tr><td>Analysis/plot.py</td><td>the campaign's own per-method figure, patched for both
        methods and for seaborn 0.11</td></tr>
      <tr><td>Analysis/data/master_long.csv</td><td>tidy: one row per design &times; species
        &times; method &mdash; 380 rows, 95 designs</td></tr>
      <tr><td>Analysis/data/positive_control.csv</td><td>the removed HuPDL1-Fc rows, parked for
        reference</td></tr>
      <tr><td>Analysis/data/auc_table.csv</td><td>every predictor &times; task with CI, partial
        AUC, precision@K, enrichment</td></tr>
      <tr><td>Analysis/data/head_to_head.csv</td><td>DeLong tests on paired &Delta;AUC</td></tr>
      <tr><td>Analysis/data/budget_matched.csv</td><td>binders found per strategy at fixed
        budget</td></tr>
      <tr><td>Analysis/data/binder_ranks.csv</td><td>every true binder's percentile under each
        method</td></tr>
      <tr><td>Analysis/REPORT.md</td><td>this write-up, in the repo</td></tr>
    </table>
  </div>
<pre>python3 Analysis/build_dataset.py
python3 Analysis/plot_predictive_power.py      # analysis + figures 1-6
python3 Analysis/plot.py --human-csv Analysis/data/merged_human.csv \\
    --mouse-csv Analysis/data/merged_mouse.csv \\
    --classification-csv PD1_duo_DesignSummary.csv --mode both \\
    --out Analysis/figures/method_benchmark.png</pre>
</section>

<footer>
  95 designs, 380 predicted complexes, 3 tasks, 22 candidate scores. Confidence intervals are
  stratified bootstrap (4000 replicates); cross-validated estimates are 20&times;5-fold.
  One loose end: <span class="mono">plot.py</span> labels the species axis
  <span class="mono">PDL1</span> &mdash; the target is PD-1.
</footer>

</div>
"""
open(OUT, "w").write(html)
print("wrote", OUT, f"{os.path.getsize(OUT)/1e6:.2f} MB")
