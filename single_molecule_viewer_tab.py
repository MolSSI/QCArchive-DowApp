import streamlit as st
import pandas as pd
import py3Dmol 
from calculations_display_tab import CalculationsDisplayTab
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import GetPeriodicTable

from archived import ArchiveManager

import re

from rdkit.Chem import rdmolfiles

import os

class SingleMoleculeViewerTab(CalculationsDisplayTab):
    def __init__(self, dataset_name, dataset_type):
        super().__init__(dataset_name, dataset_type)
        self.entry_selection_key = "smv_selected_entry"  # Use unique session key

    def show_3d_viewer(self, molecule):
        """Render a 3D molecule viewer using XYZ format."""
        try:
            xyz = molecule.to_string("xyz")  # This works in your case

            viewer = py3Dmol.view(width=600, height=400)
            viewer.addModel(xyz, "xyz")  # Format is still XYZ
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()

            # Use _make_html instead of render()
            st.components.v1.html(viewer._make_html(), height=500)
        except Exception as e:
            st.warning(f"Unable to render 3D view: {e}")

    def show_2d_image(self, smiles=None):
        """Render a 2D image of the molecule using SMILES if available."""
        # st.write(type(smiles))
        try:
            if smiles and smiles != None:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(600, 500))
                    st.image(img)
                else:
                    st.warning("Could not parse SMILES into a valid molecule.")
            else:
                st.info("No 2D image available — no SMILES string found for this entry.")
        except Exception as e:
            st.warning(f"Error generating 2D structure: {e}")

    def show_molecule_info(self, molecule, smiles=""):
        """Display molecule information in smaller text."""
        st.markdown("#### Molecule Info")  # Smaller heading

        formula = molecule.get_molecular_formula()
        weight = molecule.molecular_weight()
        symbols = molecule.symbols
        charge = molecule.molecular_charge
        multiplicity = molecule.molecular_multiplicity
        num_atoms = len(symbols)
        heavy_atoms = sum(1 for s in symbols if s != "H")

        col1, col2, col3 = st.columns(3)
        col1.caption(f"**Formula:** {formula}")
        col2.caption(f"**Weight:** {weight:.2f} g/mol")
        col3.caption(f"**Atoms:** {num_atoms}")

        col4, col5, col6 = st.columns(3)
        col4.caption(f"**Heavy Atoms:** {heavy_atoms}")
        col5.caption(f"**Charge:** {charge}")
        col6.caption(f"**Multiplicity:** {multiplicity}")

        if smiles:
            st.caption(f"**SMILES:** `{smiles}`")

    def show_molecule(self):
        records_df = self.load_records()
        if records_df is None:
            return

        st.markdown("## Select a Molecule")
        unique_entries = records_df["Entry Name"].dropna().unique()

        selected_entry = st.selectbox(
            "Select entry name",
            sorted(unique_entries),
            key=self.entry_selection_key
        )

        filtered_df = records_df[records_df["Entry Name"] == selected_entry]

        if filtered_df.empty:
            st.warning("No records found for the selected entry.")
            return

        try:
            entry = self.dataset.get_entry(selected_entry)
            if st.session_state.dataset_info['dataset_type'] == "singlepoint":
                molecule = entry.molecule
            elif st.session_state.dataset_info['dataset_type'] == "optimization":
                molecule = entry.final_molecule

            raw_smiles = filtered_df["SMILES"].values[0] if "SMILES" in filtered_df.columns else None
            smiles = str(raw_smiles) if isinstance(raw_smiles, str) or pd.notna(raw_smiles) else ""

            # === Viewer Tabs ===
            tab2d, tab3d, tab_xyz = st.tabs(["2D Viewer", "3D Viewer", "XYZ Coordinates"])

            with tab2d:
                self.show_2d_image(smiles=smiles)

            with tab3d:
                try:
                    self.show_3d_viewer(molecule)
                except Exception as e:
                    st.warning(f"Could not render 3D structure: {e}")

            with tab_xyz:
                symbols = molecule.symbols
                coords = molecule.geometry.reshape(-1, 3)

                xyz_df = pd.DataFrame({
                    "Element": symbols,
                    "X": coords[:, 0],
                    "Y": coords[:, 1],
                    "Z": coords[:, 2]
                })
                st.dataframe(xyz_df, hide_index=True)

                xyz_string = molecule.to_string("xyz")
                st.download_button(
                    label="Download XYZ File",
                    data=xyz_string,
                    file_name=f"{selected_entry}.xyz",
                    mime="text/plain"
                )

            self.show_molecule_info(molecule, smiles=smiles)

            # self.add_entry_to_other_dataset(selected_entry, molecule)

        except Exception as e:
            st.warning(f"Could not retrieve molecular info: {e}")

        # Record table
        st.markdown("#### Records")
        st.dataframe(filtered_df, hide_index=True)


    def add_entry_to_other_dataset(self, selected_entry, molecule):
        """Add the current entry to another dataset, excluding archived datasets and prioritizing recent ones."""
        st.markdown("#### Add This Entry to Another Dataset")

        metadata_file = f"./users/{st.session_state.username}/App.metadata/dataset_list.csv"
        if not os.path.exists(metadata_file):
            st.info("No dataset metadata found.")
            return

        df = pd.read_csv(metadata_file)

        # Load archive manager
        archive_manager = ArchiveManager()
        dataset_list = df.to_dict(orient="records")
        active_datasets = archive_manager.filter_active_datasets(dataset_list)

        # Load recent list
        recent_path = os.path.join(f"./users/{st.session_state.username}/App.metadata", "recent_datasets.txt")
        recent = []
        if os.path.exists(recent_path):
            with open(recent_path, "r") as f:
                recent = f.read().splitlines()

        # Sort active datasets by recency
        def sort_key(ds):
            return recent.index(ds["dataset_name"]) if ds["dataset_name"] in recent else float('inf')

        active_datasets.sort(key=sort_key)

        # Build dropdown options
        dataset_options = [
            f"{row['dataset_name']} (Type: {row['dataset_type']}, Records: {row['record_count']})"
            for row in active_datasets
        ]

        selected_dataset_label = st.selectbox(
            "Select target dataset:",
            dataset_options,
            key="target_dataset_select"
        )

        if selected_dataset_label:
            selected_index = dataset_options.index(selected_dataset_label)
            selected_dataset_info = active_datasets[selected_index]
            dataset_type = selected_dataset_info["dataset_type"]
            dataset_name = selected_dataset_info["dataset_name"]

            if st.button("Add Entry to Selected Dataset"):
                try:
                    target_dataset = st.session_state.client.get_dataset(dataset_type, dataset_name)
                    entry = self.dataset.get_entry(selected_entry)
                    target_dataset.add_entry(selected_entry, molecule)
                    st.success(f"Successfully added entry `{selected_entry}` to dataset `{dataset_name}`.")
                except Exception as e:
                    st.error(f"Failed to add entry: {e}")


    def run(self):
        # Navigation button to go back to dataset selection
        records_file, last_updated = self.get_latest_records_file()

        if records_file:
            st.write(f"**Last Updated:** {last_updated.strftime('%Y-%m-%d %H:%M:%S')}")  # Display last update time
            self.compare_local_remote(records_file)
        else:
            st.write("No records found locally. Please retrieve them from QC Archive:")

        
        if st.button("Retrieve Records from QCArchive", key=2):
            self.save_records()
            self.save_filters_metadata()  # New function to store filter metadata
            st.rerun()
        
        self.show_molecule()  # Display molecule page



class OptimizationMoleculeViewerTab(CalculationsDisplayTab):
    def __init__(self, dataset_name, dataset_type):
        super().__init__(dataset_name, dataset_type)
        self.entry_selection_key = "optmv_entry"
        self.spec_key = "optmv_spec"

    def show_3d_viewer(self, molecule):
        """Render a 3D molecule viewer using XYZ format."""
        try:
            xyz = molecule.to_string("xyz")
            viewer = py3Dmol.view(width=600, height=400)
            viewer.addModel(xyz, "xyz")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()
            st.components.v1.html(viewer._make_html(), height=500)
        except Exception as e:
            st.warning(f"Unable to render 3D view: {e}")

    def show_2d_viewer(self, smiles):
        """Render a 2D viewer from SMILES."""
        try:
            if smiles and smiles != "None":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img, caption="2D Structure from SMILES")
                else:
                    st.info("Could not parse SMILES into molecule.")
            else:
                st.info("No SMILES available for 2D rendering.")
        except Exception as e:
            st.warning(f"Error generating 2D view: {e}")

    def show_molecule_info(self, molecule, smiles=""):
        """Display molecule information in smaller text."""
        st.markdown("#### Molecule Info")  # Smaller heading

        formula = molecule.get_molecular_formula()
        weight = molecule.molecular_weight()
        symbols = molecule.symbols
        charge = molecule.molecular_charge
        multiplicity = molecule.molecular_multiplicity
        num_atoms = len(symbols)
        heavy_atoms = sum(1 for s in symbols if s != "H")

        col1, col2, col3 = st.columns(3)
        col1.caption(f"**Formula:** {formula}")
        col2.caption(f"**Weight:** {weight:.2f} g/mol")
        col3.caption(f"**Atoms:** {num_atoms}")

        col4, col5, col6 = st.columns(3)
        col4.caption(f"**Heavy Atoms:** {heavy_atoms}")
        col5.caption(f"**Charge:** {charge}")
        col6.caption(f"**Multiplicity:** {multiplicity}")

        if smiles:
            st.caption(f"**SMILES:** `{smiles}`")


    def show_molecule(self):
        """Main UI for selecting and displaying an optimized molecule."""
        records_df = self.load_records()
        if records_df is None or records_df.empty:
            st.warning("No records available.")
            return

        st.markdown("### Select a Molecule (Optimization)")

        # STEP 1: Select Entry Name
        entry_names = records_df["Entry Name"].dropna().unique()
        selected_entry = st.selectbox(
            "Select entry name",
            sorted(entry_names),
            key=self.entry_selection_key
        )

        filtered_df = records_df[records_df["Entry Name"] == selected_entry]

        if filtered_df.empty:
            st.warning("No records found for selected entry.")
            return

        # STEP 2: Select Specification
        spec_names = filtered_df["Method/Basis Set"].dropna().unique()
        selected_spec = st.selectbox(
            "Select Method/Basis",
            sorted(spec_names),
            key=self.spec_key
        )

        record_row = filtered_df[filtered_df["Method/Basis Set"] == selected_spec]
        if record_row.empty:
            st.warning("No matching record found.")
            return

        try:
            entry_name = str(record_row["Entry Name"].values[0])
            spec_name = str(record_row["Method/Basis Set"].values[0])

            record = self.dataset.get_record(entry_name=entry_name, specification_name=spec_name)
            final_mol = record.final_molecule

            # Try getting SMILES if available
            smiles = final_mol.extras.get("smiles", "")

            tab2d, tab3d, tab_xyz = st.tabs(["2D Viewer", "3D Viewer", "XYZ Coordinates"])

            with tab2d:
                self.show_2d_viewer(smiles)

            with tab3d:
                self.show_3d_viewer(final_mol)

            with tab_xyz:
                # Prepare XYZ table
                symbols = final_mol.symbols
                coords = final_mol.geometry.reshape(-1, 3)

                xyz_df = pd.DataFrame({
                    "Element": symbols,
                    "X": coords[:, 0],
                    "Y": coords[:, 1],
                    "Z": coords[:, 2]
                })

                st.dataframe(xyz_df, hide_index=True)

                # Download button
                xyz_string = final_mol.to_string("xyz")
                st.download_button(
                    "Download Final XYZ",
                    xyz_string,
                    file_name=f"{selected_entry}_{selected_spec}_final.xyz",
                    mime="text/plain"
                )

            self.show_molecule_info(final_mol)

            # self.transfer_molecule(final_mol, selected_entry, selected_spec)

        except Exception as e:
            st.error(f"Error retrieving record: {e}")

        st.markdown("### Records for Selected Entry")
        st.dataframe(filtered_df, hide_index=True)

    def transfer_molecule(self, molecule, entry_name):
        """Transfer selected molecule to another dataset."""
        st.markdown("### Transfer Molecule to Another Dataset")

        metadata_file = f"./users/{st.session_state.username}/App.metadata/dataset_list.csv"
        if not os.path.exists(metadata_file):
            st.info("No dataset metadata found.")
            return

        df = pd.read_csv(metadata_file)

        # Load archive manager
        from archived import ArchiveManager
        archive_manager = ArchiveManager()
        dataset_list = df.to_dict(orient="records")
        active_datasets = archive_manager.filter_active_datasets(dataset_list)

        # Load recent list
        recent_path = os.path.join(f"./users/{st.session_state.username}/App.metadata", "recent_datasets.txt")
        recent = []
        if os.path.exists(recent_path):
            with open(recent_path, "r") as f:
                recent = f.read().splitlines()

        # Sort active datasets by recency
        def sort_key(ds):
            return recent.index(ds["dataset_name"]) if ds["dataset_name"] in recent else float('inf')

        active_datasets.sort(key=sort_key)

        # Build dropdown options
        dataset_options = [
            f"{row['dataset_name']} (Type: {row['dataset_type']}, Records: {row['record_count']})"
            for row in active_datasets
        ]

        selected_label = st.selectbox(
            "Select Target Dataset",
            dataset_options,
            key="opt_target_ds"
        )

        if selected_label:
            selected_index = dataset_options.index(selected_label)
            selected_dataset_info = active_datasets[selected_index]
            dataset_type = selected_dataset_info["dataset_type"]
            dataset_name = selected_dataset_info["dataset_name"]

            if st.button("Transfer Molecule"):
                try:
                    target_dataset = st.session_state.client.get_dataset(dataset_type, dataset_name)
                    target_dataset.add_entry(entry_name, molecule)
                    st.success("Molecule transferred successfully.")
                except Exception as e:
                    st.error(f"Transfer failed: {e}")


    def run(self):
        # Navigation and retrieval
        records_file, last_updated = self.get_latest_records_file()

        if records_file:
            st.write(f"**Last Updated:** {last_updated.strftime('%Y-%m-%d %H:%M:%S')}")
            self.compare_local_remote(records_file)
        else:
            st.write("No records found locally. Please retrieve them from QCArchive:")


        if st.button("Retrieve Records from QCArchive", key="3"):
            self.save_records()
            self.save_filters_metadata()
            st.rerun()

        self.show_molecule()
