#!/usr/bin/env python3
"""Extract one chain from a deposited mmCIF into an AF3-ready template mmCIF.

This is how 9HAC_chainA.cif was produced from the deposited 9HAC.cif, for
AF3_DuoDesign/build_input-msa.py's --template-cif (see submit_bbf14_template.sh).

Usage (needs gemmi, which the system python3 does not have -- the protenix
conda env does):
    /work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne/conda_envs/protenix/bin/python3 \
        make_template_cif.py 9HAC.cif A 9HAC_chainA.cif

Why the deposited file can't just be handed to AF3 as-is, and why the two
date items below are re-attached by hand:

  - alphafold3/data/templates.py's get_polymer_features() raises "The
    structure must be filtered to a single chain" unless the template mmCIF
    holds exactly one polymer chain, so the other polymer chains (for 9HAC:
    the bound binder, chains B/D) and every ligand/water have to be dropped
    first. AF3's own data-pipeline-derived templates are single-chain too.

  - The same function raises "The structure must have a release date."
    when the structure has none, and
    alphafold3/structure/mmcif.py's get_release_date() reads exactly one
    item for that: min(_pdbx_audit_revision_history.revision_date). Building
    a fresh gemmi.Structure (rather than editing the parsed one in place)
    is what keeps the output clean and self-consistent -- single entity,
    regenerated _entity/_entity_poly/_chem_comp/_struct_asym -- but it also
    drops every header category, so the dates must be copied over
    explicitly or AF3 rejects the template at featurisation time.
    _pdbx_database_status.recvd_initial_deposition_date is carried along
    too: AF3 doesn't read it, but AF3's own mmCIF writer emits both items
    together (structure/structure.py), so this matches the shape of a
    template AF3 would have produced itself.

Note that the residue numbering (label_seq_id) of the source file is kept
as-is, so --template-seqid-start / --template-length stay meaningful against
the deposited numbering.
"""

import sys

import gemmi

DATE_ITEMS = (
    "_pdbx_audit_revision_history.revision_date",
    "_pdbx_database_status.recvd_initial_deposition_date",
)


def source_dates(cif_path):
    """Earliest value of each date item in the source mmCIF.

    Both are loops/pairs in a deposited file (revision_date has one row per
    revision); AF3 itself takes the oldest revision date, so do the same.
    """
    block = gemmi.cif.read(cif_path).sole_block()
    dates = {}
    for item in DATE_ITEMS:
        values = [v for v in block.find_values(item) if v not in (".", "?")]
        if values:
            dates[item] = min(v.strip("'\"") for v in values)
    return dates


def extract_chain(in_path, auth_chain_id, out_path):
    src = gemmi.read_structure(in_path)
    src_model = src[0]
    if auth_chain_id not in (ch.name for ch in src_model):
        sys.exit(f"ERROR: no chain '{auth_chain_id}' in {in_path} "
                 f"(has: {', '.join(ch.name for ch in src_model)})")
    src_chain = src_model[auth_chain_id]

    # The polymer subchain shares the chain's label_asym_id; ligands/waters
    # sitting under the same auth chain get their own subchain ids, and are
    # dropped here.
    polymer_subchain = src_chain.subchains()[0].subchain_id()

    st = gemmi.Structure()
    st.name = f"{src.name}_{auth_chain_id}"
    st.cell = src.cell
    st.spacegroup_hm = src.spacegroup_hm

    model = gemmi.Model("1")
    chain = gemmi.Chain(auth_chain_id)
    for res in src_chain:
        if res.subchain == polymer_subchain:
            chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()

    doc = st.make_mmcif_document()
    block = doc.sole_block()
    for item, value in source_dates(in_path).items():
        block.set_pair(item, value)
    doc.write_file(out_path)

    return len(chain), chain[0].seqid.num, chain[-1].seqid.num


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <input.cif> <auth_chain_id> <output.cif>")
    in_path, auth_chain_id, out_path = sys.argv[1:]

    dates = source_dates(in_path)
    if DATE_ITEMS[0] not in dates:
        sys.exit(f"ERROR: {in_path} has no {DATE_ITEMS[0]} -- AF3 would reject "
                 "the resulting template with 'The structure must have a release date.'")

    n_res, first, last = extract_chain(in_path, auth_chain_id, out_path)
    print(f"Wrote {out_path}: chain {auth_chain_id}, {n_res} residues "
          f"(label_seq_id {first}-{last})")
    for item, value in dates.items():
        print(f"  {item} = {value}")


if __name__ == "__main__":
    main()
