# get_molecule_tab.py
import streamlit as st
import pandas as pd
import tempfile
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw
import pubchempy as pcp  # Optional: used for naming if needed
import time
import numpy as np

class GetMoleculeTab:
    def __init__(self):
        RDLogger.DisableLog('rdApp.*')
        
    def convert_smiles_qcmol(self, smile: str):
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) == -1:
            raise ValueError("Failed to generate a 3D conformer for SMILES: " + smile)
        AllChem.UFFOptimizeMolecule(mol)
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        conf = mol.GetConformer()
        coords = list(conf.GetPositions().flatten())
        bohr_unit = 1.8897259886
        coords_bohr = [x * bohr_unit for x in coords]
        import qcelemental as qcel
        qc_mol = qcel.models.Molecule(
            symbols=symbols,
            #geometry=coords,
            geometry=coords_bohr,
            molecular_charge=Chem.GetFormalCharge(mol),
            extras={"smiles": smile}
        )
        return qc_mol

    def run_csv(self, key_suffix="default"):
        st.title("Add Molecule via CSV")
        st.markdown("""
        ### 📄 CSV Upload Guidelines
        - Upload a valid `.csv` file containing one or more molecules.
        - The file **must include**:
        - A `SMILES` column (case-insensitive) with valid SMILES strings.
        - A `Name` column (case-insensitive) with molecule names.
        - Additional columns (e.g., molecular weight, double bond counts) are optional and supported.

        **Example CSV Format:**
                    """)
        st.image("gif/csv.png", caption="Example CSV Format for Upload")
        uploaded_file = st.file_uploader("Upload CSV file", type=["csv"], key=f"csv_upload_{key_suffix}")
        if uploaded_file is not None:
            try:
                csv_df = pd.read_csv(uploaded_file)
            except Exception as e:
                st.error(f"Error reading CSV file: {e}")
                return
            

            name_column = None
            smiles_column = None
            for col in csv_df.columns:
                if col.lower() =="name":
                    name_column = col
                if col.lower() == "smiles":
                    smiles_column = col

            if name_column:
                st.success(f"Found name column: '{name_column}'")
            else:
                st.warning("No 'name' column found. Molecule names will not be used. " \
                "Please add a name column if you have custom names for molecules")
            if smiles_column:
                st.success(f"Found smiles column: '{smiles_column}'")
            else:
                st.error("No 'SMILES' column found, CSV must contain a column named 'SMILES' (case-insensitive).")
            st.subheader("CSV Data")
            st.dataframe(csv_df)
            if smiles_column:
                user_choice = st.radio("Choose an option", 
                                        ["Process one SMILE", "Process Multiple SMILES", "Process all SMILES"])
                if user_choice == "Process one SMILE":
                    if name_column:
                        options = {
                            f"{row[name_column]} (Index {idx})": (row[smiles_column], row[name_column])
                            for idx, row in csv_df.iterrows()
                        }
                    else:
                        options = {
                            f"SMILES {idx}": (row[smiles_column], None)
                            for idx, row in csv_df.iterrows()
                        }
                    
                    selected_label = st.selectbox("Select entry:", list(options.keys()))
                    selected_smiles, selected_name = options[selected_label]
                    #selected_index = list(options.keys()).index(selected_label)

                    selected_name = selected_label.split(" (Index")[0] if name_column else None
                    # Create RDKit mol
                    mol = Chem.MolFromSmiles(selected_smiles)
                    if mol:
                        img = Draw.MolToImage(mol, size=(300, 300))
                        st.image(img, caption=selected_name)
                        
                        # Save mol and SMILES for preview
                        #st.session_state.sdf_mol_pairs = [(selected_smiles, mol)]

                        if st.button("Process Selected Molecule", key="csv_process_selected"):
                            st.session_state.molecules = [(selected_name, selected_smiles)]
                            st.success("Molecule added successfully! Proceed to the Request Calculation tab")
                    else:
                        st.error("Invalid SMILES string — could not generate molecule.")
                     # Extract name from label if available
                    st.write(f"**Selected SMILES:** `{selected_smiles}`")
                    if selected_name:
                        st.write(f"**Name:** `{selected_name}`")
                    #st.session_state.smile_code1 = selected_smiles

            # --------------------------------------------------------------------------------------------------------
                elif user_choice == "Process Multiple SMILES":
                    # Build label → (smiles, name, index) mapping
                    if name_column:
                        options = {
                            f"{row[name_column]} (Index {idx})": (row[smiles_column], row[name_column], idx)
                            for idx, row in csv_df.iterrows()
                        }
                    else:
                        options = {
                            f"SMILES {idx}": (row[smiles_column], None, idx)
                            for idx, row in csv_df.iterrows()
                        }

                    selected_labels = st.multiselect("Select multiple entries:", list(options.keys()))

                    # Extract selected tuples
                    selected_data = [(name, smiles, idx) for label in selected_labels
                                    for smiles, name, idx in [options[label]]]

                    st.write(f"**Selected {len(selected_data)} molecules:**")

                    # Show names, SMILES, and expandable images
                    for name, smiles, idx in selected_data:
                        st.write(f"- **{name or 'Unnamed'}**: `{smiles}`")
                        with st.expander(f"Click to view {name or 'Unnamed'}"):
                            mol = Chem.MolFromSmiles(smiles)
                            if mol:
                                img = Draw.MolToImage(mol, size=(300, 300))
                                st.image(img, caption=f"{name or 'Unnamed'} ({smiles})")
                            else:
                                st.error("Invalid SMILES string.")

                    # Save molecule data to session
                    if st.button("Process Selected SMILES", key="csv_process_multiple"):
                        st.session_state.molecules = [(name, smiles) for name, smiles, _ in selected_data]
                        st.success("Molecules added successfully! Proceed to the Request Calculation tab")
                
            # -------------------------------------------------------------------------------------------------------
                elif user_choice == "Process all SMILES":
                    all_smiles = csv_df[smiles_column].tolist()
                    all_indices = list(range(len(csv_df)))  # since no "Index" column is assumed in CSV

                    if name_column:
                        all_names = csv_df[name_column].tolist()
                    else:
                        all_names = [None] * len(all_smiles)

                    st.write(f"**Selected {len(all_smiles)} molecules**")

                    # Save name + SMILES in molecules list (no image preview)
                    if st.button("Process All SMILES", key="csv_process_all"):
                        st.session_state.molecules = list(zip(all_names, all_smiles))
                        st.success("All molecules added successfully! Proceed to the Request Calculation tab")

            else:
                st.error("CSV file must contain SMILES and Name columns.")
    
    def convert_mol_qcmol(self, smiles, mol):
        # Extract element symbols and 3D coordinates
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        conf = mol.GetConformer()
        coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        bohr_unit = 1.8897259886
        coords_bohr = [x * bohr_unit for x in coords]
        import qcelemental as qcel
        qc_mol = qcel.models.Molecule(
            symbols=symbols,
            #geometry=coords,
            geometry=coords_bohr,
            molecular_charge=Chem.GetFormalCharge(mol),
            extras={"smiles": smiles}
        )
        return qc_mol
    
    def run_sdf(self, key_suffix="default"):
        st.title("Add Molecule via SDF")
           # -- NEW: SDF Guidelines Section --
        st.markdown("""
        ### 📄 SDF Upload Guidelines
        - Upload a valid `.sdf` file containing one or more molecules.
        - **Custom molecule names** should be included in the `<name>` field of the SDF properties.
        - Molecules without a `<name>` entry will automatically be assigned a generic name (e.g., `Molecule_0`, `Molecule_1`, etc.).
        - Ensure that the SDF file has the correct atom and bond block definitions.

        Example of an SDF molecule entry:

        ```
        Fluorometholone Acetate
            RDKit          2D

        30 33  0  0  0  0  0  0  0  0999 V2000
        -3.0000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...
        M  END
        >  <name>
        Fluorometholone Acetate

        $$$$
        ```

        **Notes:**
        - Properties like `MW`, `logP`, and `Canonical SMILES` can also be included and will be automatically extracted.
        """)
        # -- END of Guidelines Section --
        uploaded_sdffile = st.file_uploader("Upload SDF file", type=["sdf"], key=f"sdf_upload_{key_suffix}")
        if uploaded_sdffile is not None:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as temp_file:
                temp_file.write(uploaded_sdffile.getvalue())
                temp_filepath = temp_file.name
            suppl = Chem.SDMolSupplier(temp_filepath, removeHs=False)
            mols = [mol for mol in suppl if mol is not None]
            if not mols:
                st.error("No valid molecules found in the SDF file.")
            else:
                st.success(f"Loaded {len(mols)} molecules from the SDF file.")
                data = []
                for i, mol in enumerate(mols):
                    smiles = Chem.MolToSmiles(mol) # retrieve smiles
                    props = mol.GetPropsAsDict()
                    lower_props = {k.lower(): v for k, v in props.items()} #normalize all keys to lower case
                     # Get name from "name" (case-insensitive), _Name, or fallback
                    name = lower_props.get("name", mol.GetProp("_Name") if mol.HasProp("_Name") else f"Molecule_{i}")
                    props["SMILES"] = smiles
                    props["Index"] = i
                    props["NAME"] = name
                    data.append(props)
                sdf_df = pd.DataFrame(data)
                st.subheader("Extracted Molecule Data")
                st.dataframe(sdf_df)
                if "SMILES" in sdf_df.columns:
                    user_choice = st.radio("Choose an option", 
                                           ["Process one molecule", "Process multiple molecules", "Process all molecules"])
                    

                    if user_choice == "Process one molecule":
                        options = {
                        f"{row['NAME']} (Index {row['Index']})": row['Index']
                        for _, row in sdf_df.iterrows()
                        }
                        selected_label = st.selectbox("Select Molecule:", options.keys())
                        selected_index = options[selected_label]
                        selected_row = sdf_df.loc[sdf_df["Index"] == selected_index].iloc[0]
                        selected_smiles = selected_row["SMILES"]
                        selected_name = selected_row["NAME"]
                        st.write(f"**Selected SMILES:** '{selected_smiles}'")
                        img = Draw.MolToImage(mols[selected_index], size=(300, 300))
                        st.image(img, caption="Selected Molecule")
                        #for sdf changes
                        st.session_state.sdf_mol_pairs = [(selected_smiles, mols[selected_index])] # not necessary
                        
                        if st.button("Process Selected Molecule", key="sdf_process_selected"):
                            st.session_state.sdfmolecules = [(selected_name, selected_smiles, mols[selected_index])]
                            st.success("Molecule added successfully! Proceed to the Request Calculation tab")

                    elif user_choice == "Process multiple molecules":
                        # Create options with labels showing both NAME and Index
                        options = {
                            f"{row['NAME']} (Index {row['Index']})": row['Index']
                            for _, row in sdf_df.iterrows()
                        }

                        selected_labels = st.multiselect(
                            "Select multiple molecules:",
                            options.keys()
                        )

                        selected_indices = [options[label] for label in selected_labels]

                        if selected_indices:
                            selected_rows = sdf_df[sdf_df["Index"].isin(selected_indices)]
                            selected_smiles_list = selected_rows["SMILES"].tolist()
                            selected_name_list = selected_rows["NAME"].tolist()

                            st.write(f"**Selected {len(selected_smiles_list)} molecules:**")
                            #for name, smiles in zip(selected_name_list, selected_smiles_list):
                                #st.write(f"- **{name}**: `{smiles}`")

                            for idx, name, smiles in zip(selected_indices, selected_name_list, selected_smiles_list):
                                st.write(f"- **{name}**: `{smiles}`")
                                with st.expander(f"Click to view{name}"):
                                    mol = mols[idx]
                                    img = Draw.MolToImage(mol, size=(300, 300))
                                    st.image(img, caption=f"{name} ({smiles})")

                            st.session_state.sdf_mol_pairs = [
                                (smiles, mols[i]) for smiles, i in zip(selected_smiles_list, selected_indices)
                            ]

                            if st.button("Process Selected Molecules", key="sdf_process_multiple"):
                                st.session_state.sdfmolecules = [
                                    (name, smiles, mols[i]) for name, smiles, i in zip(selected_name_list, selected_smiles_list, selected_indices)
                                ]
                                st.success("Molecules added successfully! Proceed to the Request Calculation tab")
                    
                    elif user_choice == "Process all molecules":
                        all_indices = sdf_df["Index"].tolist()

                        st.write(f"**Selected {len(all_indices)} molecules**")

                        # Save all molecules into sdf_mol_pairs individually (smiles, mol)
                        all_smiles = sdf_df["SMILES"].tolist()

                        st.session_state.sdf_mol_pairs = [
                            (smiles, mols[idx]) for smiles, idx in zip(all_smiles, all_indices)
                        ]

                        if st.button("Process All Molecules", key="sdf_process_all"):
                            all_names = sdf_df["NAME"].tolist()
                            st.session_state.sdfmolecules = [
                                (name, smiles, mols[idx]) for name, smiles, idx in zip(all_names, all_smiles, all_indices)
                            ]
                            st.success("All molecules added successfully! Proceed to the Request Calculation tab")
                else:
                    st.error("SDF file must contain valid molecular structures.")

    def run(self, key_suffix="default"):
        st.title("Add Molecule(s)")
        # Create inner tabs for CSV and SDF functionality
        tab_csv, tab_sdf = st.tabs(["Add by CSV file", "Add by SDF file"])
        with tab_csv:
            self.run_csv(key_suffix=key_suffix)
        with tab_sdf:
            self.run_sdf(key_suffix=key_suffix)
