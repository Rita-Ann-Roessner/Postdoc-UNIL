"""Campaign score tables, so a figure cell can stand on its own.

`load_pdl1()` and `load_pd1()` rebuild exactly the `camp` and `pd1c` frames the
Figure 2 cells construct inline: one row per design, target and method, with

    Design  target  Experimental  Type  method  score  y  rank

`score` is oriented higher-is-better for every method (MD carries -deltaG), so
the bottleneck of an orthologue pair is always the minimum. Replicas and seeds
are averaged before anything else, as the protocol specifies.

This duplicates the loading logic in the Figure 2 cells. If one changes, change
the other, or point those cells here too.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PDL1 = Path('/Users/roessner/Documents/PostDoc/Data/261107_funnel_pdl1_reprediction')
PD1  = Path('/Users/roessner/Documents/PostDoc/Data/261111_reprediction_PD1_designs_LN')

PREDICTORS = ["AF3", "ESMFold2", "Protenix-v2", "OpenDDE"]
METHODS    = ["MD"] + PREDICTORS
TARGETS    = ["human", "mouse"]
COLS       = ['Design', 'target', 'Experimental', 'Type', 'method', 'score']


def _finish(df):
    df = df.copy()
    df['y']    = (df.Type == 'binder').astype(int)
    df['rank'] = df.groupby(['method', 'target'])['score'].rank(pct=True)
    return df


# --- PD-L1 ---------------------------------------------------------------
PDL1_DIR = {"AF3": "AF3_repredict", "ESMFold2": "ESMFold2_repredict",
            "Protenix-v2": "Protenix-v2_repredict", "OpenDDE": "OpenDDE_repredict"}
PDL1_SUB = {"human": "hPDL1", "mouse": "mPDL1"}


def load_pdl1():
    """The PD-L1 campaign: 27 designs against both orthologues."""
    anno = pd.read_csv(PDL1 / 'experimental_labels.csv')

    def binary(df):
        """binder = bound this target ('both' or this species); drops 'unknown'."""
        df = df.merge(anno, how='left', on='Design')
        df = df[df.Experimental != 'unknown'].copy()
        df['Type'] = np.select(
            [df.Experimental.eq('none'),
             df.Experimental.eq('both'),
             df.Experimental.eq(df.target)],
            ['decoy', 'binder', 'binder'], default='decoy')
        return df

    md = pd.read_pickle(PDL1 / 'deltaG_labled.pkl')
    md = md[md.Experimental != 'unknown']
    md = md.groupby(['Design', 'target', 'Experimental'], as_index=False).agg(
        deltaG=('deltaG', 'mean'))                    # protocol: mean over replicas
    md['Type'] = np.select(
        [md.Experimental.eq('none'),
         md.Experimental.eq('both'),
         md.Experimental.eq(md.target)],
        ['decoy', 'binder', 'binder'], default='decoy')
    md = md.assign(method='MD', score=-md.deltaG)

    def iptm(name):
        d = pd.concat(
            [pd.read_csv(PDL1 / PDL1_DIR[name] / PDL1_SUB[t] / 'output.csv')[['ID', 'iptm']]
               .groupby('ID', as_index=False).agg(iptm=('iptm', 'mean'))   # mean over seeds
               .rename(columns={'ID': 'Design'}).assign(target=t)
             for t in TARGETS], ignore_index=True)
        return binary(d).assign(method=name, score=lambda x: x.iptm)

    return _finish(pd.concat([md[COLS]] + [iptm(m)[COLS] for m in PREDICTORS],
                             ignore_index=True))


# --- PD-1 ----------------------------------------------------------------
# Label convention from the PD-1 Monitor notebook: the natural-ligand control is
# not a design, sticky designs count as non-binding, and a borderline call is
# dropped only from the target it is borderline for.
PD1_DIR     = {m: m for m in PREDICTORS}
PD1_SUB     = {"human": "hPD1", "mouse": "mPD1"}
PD1_CONTROL = 'Positive control (natural ligand)'
PD1_EXP = {'Cross-reactive': 'both', 'Human-specific': 'human',
           'Mouse-specific': 'mouse', 'Non-binding': 'none',
           'Likely nonspecific/sticky': 'none', 'Non-binding (no MuPD-1 data)': 'none',
           'Borderline (weak human trend)': 'unknown',
           'Borderline (weak mouse trend)': 'unknown'}
PD1_BINDER = {'human': {'Cross-reactive', 'Human-specific'},
              'mouse': {'Cross-reactive', 'Mouse-specific'}}
PD1_DROP   = {'human': {'Borderline (weak human trend)'},
              'mouse': {'Borderline (weak mouse trend)', 'Non-binding (no MuPD-1 data)'}}


def load_pd1():
    """The PD-1 campaign, restricted to the 39 designs that carry a funnel dG."""
    summ = (pd.read_csv(PD1 / 'PD1_duo_DesignSummary.csv', sep=';', decimal=',')
              .rename(columns={'EPFL Name': 'Design'}))
    summ = summ[summ.Classification != PD1_CONTROL].copy()
    summ['Experimental'] = summ.Classification.map(PD1_EXP)
    anno = summ[['Design', 'Experimental']]
    long = pd.concat([
        pd.DataFrame({'Design': s.Design, 'target': t,
                      'Type': np.where(s.Classification.isin(PD1_BINDER[t]),
                                       'binder', 'decoy')})
        for t in TARGETS for s in [summ[~summ.Classification.isin(PD1_DROP[t])]]
    ], ignore_index=True)

    def binary(df):
        return df.merge(long, on=['Design', 'target'], how='inner') \
                 .merge(anno, on='Design', how='left')

    def md_one(path, target):
        d = pd.read_csv(path, index_col=0)
        d['Design'] = d.binder.str.replace(r'_m$', '', regex=True)   # mouse names carry _m
        return d.assign(target=target)[['Design', 'target', 'deltaG']]

    md = pd.concat([md_one(PD1 / 'PD1_Human_deltaG.csv', 'human'),
                    md_one(PD1 / 'PD1_Mouse_deltaG.csv', 'mouse')], ignore_index=True)
    designs = set(md.Design)                          # MD covers 39 of the 95 designs
    md = md.groupby(['Design', 'target'], as_index=False).agg(
        deltaG=('deltaG', 'mean'))                    # protocol: mean over replicas
    md = binary(md).assign(method='MD')
    md['score'] = -md.deltaG

    def iptm(name):
        """3 seeds per design, averaged, as for PD-L1. Only output.csv is read."""
        d = pd.concat(
            [pd.read_csv(PD1 / PD1_DIR[name] / PD1_SUB[t] / 'output.csv')[['ID', 'iptm']]
               .groupby('ID', as_index=False).agg(iptm=('iptm', 'mean'))
               .rename(columns={'ID': 'Design'}).assign(target=t)
             for t in TARGETS], ignore_index=True)
        d = d[d.Design.isin(designs)]                 # same design set as MD
        return binary(d).assign(method=name, score=lambda x: x.iptm)

    return _finish(pd.concat([md[COLS]] + [iptm(m)[COLS] for m in PREDICTORS],
                             ignore_index=True))


def load_campaigns():
    """[(label, frame), ...] in figure order."""
    return [("PD-L1", load_pdl1()), ("PD-1", load_pd1())]
