import pandas as pd
import os, time, traceback
from colabdesign import clear_mem
from colabdesign.mpnn import mk_mpnn_model

def mpnn_gen_sequence(trajectory_pdb, binder_chain, trajectory_interface_residues, advanced_settings):
    print("🧹 Clearing GPU memory...")
    clear_mem()

    print("⚙️ Initializing MPNN model...")
    mpnn_model = mk_mpnn_model(
        backbone_noise=advanced_settings["backbone_noise"],
        model_name=advanced_settings["model_path"],
        weights=advanced_settings["mpnn_weights"]
    )

    design_chains = 'A,' + binder_chain
    if advanced_settings["mpnn_fix_interface"]:
        fixed_positions = 'A,' + trajectory_interface_residues
        fixed_positions = fixed_positions.rstrip(",")
        print(f"🔒 Fixing interface residues: {trajectory_interface_residues}")
    else:
        fixed_positions = 'A'
        print("🔓 No interface residues fixed")

    print(f"📦 Preparing MPNN inputs for chains: {design_chains}, fixed: {fixed_positions}")
    mpnn_model.prep_inputs(
        pdb_filename=trajectory_pdb,
        chain=design_chains,
        fix_pos=fixed_positions,
    )

    print("🎲 Sampling sequences from MPNN...")
    mpnn_sequences = mpnn_model.sample(
        temperature=advanced_settings["sampling_temp"],
        num=1,
        batch=advanced_settings["num_seqs"]
    )

    print("✅ MPNN sampling complete.")
    return mpnn_sequences

# === Constants ===
input_csv_template_1 = "/Users/roessner/Documents/PostDoc/Data/260405_duo_design_CXCR4/Trajectory/trajectory_CXCR4_human_no_clashes.csv"
input_csv_template_2 = "/Users/roessner/Documents/PostDoc/Data/260405_duo_design_CXCR4/Trajectory/trajectory_CXCR4_mouse_no_clashes.csv"
template_name = "CXCR4_human" # target used for mpnn
pdb_folder = "/Users/roessner/Documents/PostDoc/Data/260405_duo_design_CXCR4/Trajectory/no_membrane_clashes"
output_csv = "/Users/roessner/Documents/PostDoc/Data/260405_duo_design_CXCR4/mpnn_seqs_new.csv"
binder_chain = "B"

# === Load design metadata ===
print(f"📥 Reading input CSVs: {input_csv_template_1}, {input_csv_template_2}")
df1 = pd.read_csv(input_csv_template_1)
df2 = pd.read_csv(input_csv_template_2)

# Keep only common designs
common_designs = set(df1["Design"]) & set(df2["Design"])

df1 = df1[df1["Design"].isin(common_designs)].copy()
df2 = df2[df2["Design"].isin(common_designs)].copy()

# Sort to ensure identical order
df1 = df1.sort_values("Design").reset_index(drop=True)
df2 = df2.sort_values("Design").reset_index(drop=True)

# === Merge interface residues ===
def merge_interface_residues(res1, res2):
    set1 = set(res1.split(",")) if pd.notna(res1) else set()
    set2 = set(res2.split(",")) if pd.notna(res2) else set()

    merged = sorted(set1 | set2)
    return ",".join(merged)

df = df1.copy()
df["DuoInterfaceResidues"] = [
    merge_interface_residues(r1, r2)
    for r1, r2 in zip(df1["InterfaceResidues"], df2["InterfaceResidues"])
]
print(f"📊 Found {len(df)} designs.")

output_rows = []

# === MPNN config ===
advanced_settings = {
    "mpnn_fix_interface": True,
    "omit_AAs": "",
    "force_reject_AA": True,
    "num_seqs": 20,
    "max_mpnn_sequences": 2,
    "sampling_temp": 0.1,
    "backbone_noise": 0.00,
    "model_path": "v_48_020",
    "mpnn_weights": "soluble",
}

# === Process each design ===
for idx, row in df.iterrows():
    design_name = row["Design"]
    trajectory_pdb = os.path.join(pdb_folder, f"{design_name}{template_name}.pdb")
    trajectory_interface_residues = row["DuoInterfaceResidues"]
    length = int(row["Length"])

    print(f"\n🚧 Processing design {idx+1}/{len(df)}: {design_name} (length={length})")

    if not os.path.isfile(trajectory_pdb):
        print(f"❌ PDB not found: {trajectory_pdb}")
        continue

    try:
        start_time = time.time()
        print(f"📂 Loading PDB: {trajectory_pdb}")

        # Run MPNN
        mpnn_trajectories = mpnn_gen_sequence(
            trajectory_pdb,
            binder_chain,
            trajectory_interface_residues,
            advanced_settings
        )

        print(f"📈 Extracting sequences of final length {length}...")
        mpnn_sequences = [{
            'seq': mpnn_trajectories['seq'][n][-length:],
            'score': mpnn_trajectories['score'][n],
            'seqid': mpnn_trajectories['seqid'][n]
        } for n in range(advanced_settings["num_seqs"])]

        print("🧮 Sorting sequences by MPNN score...")
        mpnn_sequences.sort(key=lambda x: x['score'])

        for i, seq_data in enumerate(mpnn_sequences):
            sequence_name = f"{design_name}_mpnn{i+1}"
            print(f"📌 {sequence_name}: Score = {seq_data['score']:.4f}")
            output_rows.append({
                "Name": sequence_name,
                "Design": design_name,
                "Sequence": seq_data["seq"],
                "Score": seq_data["score"],
                "SeqID": seq_data["seqid"],
                "FixedInterface": trajectory_interface_residues,
                "Length": length
            })

        elapsed = time.time() - start_time
        print(f"✅ Finished {design_name} in {elapsed:.2f} seconds.")

    except Exception as e:
        print(f"⚠️ Error generating sequences for {design_name}: {e}")
        traceback.print_exc()
        continue

# === Save output ===
print(f"\n💾 Writing output to: {output_csv}")
output_df = pd.DataFrame(output_rows)
output_df.to_csv(output_csv, index=False)
print(f"✅ Done. Total sequences written: {len(output_df)}")