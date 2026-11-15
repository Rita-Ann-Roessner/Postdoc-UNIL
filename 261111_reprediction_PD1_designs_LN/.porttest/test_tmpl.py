import json, sys, traceback
from alphafold3 import structure
from alphafold3.data import templates
from alphafold3.constants import mmcif_names

BUNDLE = "/iopsstor/scratch/cscs/rroessner/261111_reprediction_PD1_designs_LN/Target"
REV = """#
loop_
_pdbx_audit_revision_history.ordinal
_pdbx_audit_revision_history.data_content_type
_pdbx_audit_revision_history.major_revision
_pdbx_audit_revision_history.minor_revision
_pdbx_audit_revision_history.revision_date
1 'Structure model' 1 0 1970-01-01
"""

def try_parse(label, mmcif, qidx, tidx, qlen):
    try:
        s = structure.from_mmcif(mmcif, fix_mse_residues=True, fix_arginines=True,
                                 include_bonds=False, include_water=False, include_other=True)
        f = templates.get_polymer_features(
            chain=s, chain_poly_type=mmcif_names.PROTEIN_CHAIN,
            query_sequence_length=qlen,
            query_to_hit_mapping=dict(zip(qidx, tidx)))
        seq = f["template_sequence"]
        print(f"  {label}: OK  release_date={s.release_date}  template_sequence[{len(seq)}]={seq[:60]}...")
        return seq
    except Exception as e:
        print(f"  {label}: FAIL  {type(e).__name__}: {e}")
        return None

for t in ("hPD1", "mPD1"):
    print(f"=== {t} ===")
    tj = json.load(open(f"{BUNDLE}/{t}/templates/templates.json"))[0]
    seq = open(f"{BUNDLE}/{t}/{t}.fasta").read().split("\n")[1].strip()
    print(f"  target seq len={len(seq)}")
    a = try_parse("raw mmcif          ", tj["mmcif"], tj["queryIndices"], tj["templateIndices"], len(seq))
    patched = tj["mmcif"].rstrip("\n") + "\n" + REV
    b = try_parse("mmcif + rev history", patched, tj["queryIndices"], tj["templateIndices"], len(seq))
    got = b or a
    if got:
        core = (got.decode() if isinstance(got, bytes) else got).replace("-", "")
        print(f"  aligned == target sequence: {core == seq}")
