import streamlit as st
import pandas as pd
import tempfile
import time
from streamlit_ketcher import st_ketcher
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw
import qcportal as ptl
from qcportal.singlepoint import QCSpecification, SinglepointDataset
from request_calculation_tab import DatasetManager, CalculationManager
from get_molecule_tab import GetMoleculeTab
# Import the unified request calculation menu from our new module.
from request_calculation_modes import UnifiedRequestCalculationMenu

# ----------------- AddProjectPage Class -------------------
class AddProjectPage:
    def __init__(self, client):
        """
        client: an instance of ptl.PortalClient.
        """
        self.client = client
        # Initialize (or reset) session keys for molecule storage.
        if "project_molecules" not in st.session_state:
            st.session_state.project_molecules = []  # For SMILES
        if "molecules" not in st.session_state:
            st.session_state.molecules = []  # For CSV/SDF
        # Initialize keys for project information.
        if "project_name" not in st.session_state:
            st.session_state.project_name = ""
        if "project_description" not in st.session_state:
            st.session_state.project_description = ""
        # We'll use "add_molecule_mode" to switch among the SMILES/CSV/SDF options in a unified tab.
        if "add_molecule_mode" not in st.session_state:
            st.session_state.add_molecule_mode = None

    def add_molecule_smiles(self):
        """
        Renders the "Add Molecule via SMILES" functionality.
        """
        st.subheader("Add Molecule via SMILES")


        # Allow the user to input a SMILES string (or draw with Ketcher).
        smiles_input = st.text_input("Enter SMILES string:", key="smiles_input")
        smiles_name = st.text_input("Enter Molecule name:", key="smiles_name")
        if st.button("Add Molecule SMILES", key="add_smiles_btn"):
            if not final_smiles:
                st.error("Please enter a valid SMILES string.")
            else:
                get_mol_tab = GetMoleculeTab()
                try:
                    qc_mol = get_mol_tab.convert_smiles_qcmol(final_smiles)
                    st.session_state.project_molecules.append((smiles_name, qc_mol))
                    st.success("Molecule added successfully!")
                except Exception as e:
                    st.error(f"Error converting SMILES: {e}")
        
        st.markdown("---")


        st.subheader("Add Molecule via Drawing")
        # Use drawn_smiles if available; otherwise, fallback to smiles_input.
        drawn_smiles = st_ketcher(smiles_input, key="ketcher_smiles")
        smiles_name = st.text_input("Enter Molecule name:", key="ketcher_name")
        final_smiles = drawn_smiles if drawn_smiles.strip() != "" else smiles_input.strip()
        st.write(f"SMILES: {final_smiles}")
        if st.button("Add Molecule Drawing", key="add_smiles_btn_draw"):
            if not final_smiles:
                st.error("Please enter a valid SMILES string.")
            else:
                get_mol_tab = GetMoleculeTab()
                try:
                    qc_mol = get_mol_tab.convert_smiles_qcmol(final_smiles)
                    st.session_state.project_molecules.append((smiles_name, qc_mol))
                    st.success("Molecule added successfully!")
                except Exception as e:
                    st.error(f"Error converting SMILES: {e}")

    def run_unified_add_molecule(self):
        """
        Unified interface for adding molecules via SMILES, CSV, or SDF.
        Displays three buttons with short descriptions and a back button.
        """
        st.header("Add Molecule")
        if st.session_state.get("add_molecule_mode") is None:
            st.write("Choose a method to add molecules:")
            col1, col2, col3 = st.columns(3)
            if col1.button("Add via SMILES", key="unified_add_smiles"):
                st.session_state.add_molecule_mode = "smiles"
                st.rerun()
            if col2.button("Add via CSV", key="unified_add_csv"):
                st.session_state.add_molecule_mode = "csv"
                st.rerun()
            if col3.button("Add via SDF", key="unified_add_sdf"):
                st.session_state.add_molecule_mode = "sdf"
                st.rerun()
            st.markdown("---")
            st.write("**Descriptions:**")
            st.write("- **SMILES:** Input or draw a molecule using a SMILES string.")
            st.write("- **CSV:** Upload a CSV file containing molecules.")
            st.write("- **SDF:** Upload an SDF file to extract molecules.")
        else:
            mode = st.session_state.add_molecule_mode
            if mode == "smiles":
                self.add_molecule_smiles()
            elif mode == "csv":
                from get_molecule_tab import GetMoleculeTab
                get_mol_tab = GetMoleculeTab()
                get_mol_tab.run_csv(key_suffix="project_update")
            elif mode == "sdf":
                from get_molecule_tab import GetMoleculeTab
                get_mol_tab = GetMoleculeTab()
                get_mol_tab.run_sdf(key_suffix="project_update")
            if st.button("Back", key="unified_add_mol_back"):
                st.session_state.add_molecule_mode = None
                st.rerun()
        # if st.button("Back to Home", key="back_home_unified_add"):
        #     st.session_state.page = "home"
        #     st.rerun()

    def run(self):
        st.title("Add New Project")

        # --------------------------- Project Info Input ---------------------------
        st.subheader("Project Information")
        project_name = st.text_input(
            "Enter Project Name:",
            st.session_state.get("project_name", "")
        )
        project_description = st.text_area(
            "Project Description (optional):",
            st.session_state.get("project_description", "")
        )

        # Persist in session
        if project_name:
            st.session_state.project_name = project_name.strip()
        if project_description:
            st.session_state.project_description = project_description

        # 1) enforce mandatory name
        if not st.session_state.project_name:
            st.error(" A project name is required to proceed.")
            return

        # 2) check uniqueness against QCArchive — only if we haven't already created one
        if "dataset_name" not in st.session_state:
            existing = [d["dataset_name"] for d in self.client.list_datasets()]
            if st.session_state.project_name in existing:
                st.error(
                    f" A project named '{st.session_state.project_name}' already exists. "
                    "Please choose a different name."
                )
                return
        # --------------------------- Tabs for Molecule Input & Request Calculation ---------------------------
        st.subheader("Molecule Input & Calculation Request")
        tab_add, tab_request = st.tabs(["Add Molecule", "Request Calculation"])

        with tab_add:
            self.run_unified_add_molecule()

        with tab_request:
            # Base name & description now come strictly from the user
            st.session_state.dataset_base_name = st.session_state.project_name.replace(" ", "_")
            st.session_state.dataset_description = (
                st.session_state.project_description
                or f"Project '{st.session_state.project_name}' QCArchive dataset"
            )

            from request_calculation_modes import UnifiedRequestCalculationAddMenu
            UnifiedRequestCalculationAddMenu(self.client).run()

            if st.button("Home 🏠", key="back_home_reqcalc"):
                # clear everything so we can start fresh next time
                for k in [
                    "dataset_created", "dataset_base_name", "dataset_description",
                    "dataset", "dataset_name", "dataset_type",
                    "qcarchive_objects", "molecules", "project_molecules",
                    "sp_dataset_created", "opt_dataset_created"
                ]:
                    st.session_state.pop(k, None)
                st.session_state.page = "home"
                st.rerun()

