import time
import streamlit as st
import pandas as pd
import qcportal as ptl
from get_molecule_tab import GetMoleculeTab
from qcportal.singlepoint import QCSpecification, SinglepointDataset
import inspect
import json
import qcportal

# ------------------------------
# DatasetManager Class
# ------------------------------
class DatasetManager:
    def __init__(self, client, base_name, description, dataset_type="singlepoint"):
        """
        client: ptl.PortalClient instance.
        base_name: dataset name.
        description: dataset description.
        dataset_type: e.g., "singlepoint" or "optimization".
        """
        self.client = client
        self.base_name = base_name
        self.description = description
        self.dataset_type = dataset_type
        self.dataset = None

    def dataset_exists(self):
        existing_datasets = self.client.list_datasets()
        existing_names = [ds.get("dataset_name") for ds in existing_datasets]
        base_norm = self.base_name.lower().strip()
        return any(name is not None and name.lower().strip() == base_norm for name in existing_names)

    def create_or_recreate_dataset(self, new_entry):
        if self.dataset_exists():
            # Retrieve the existing dataset.
            existing_dataset = self.client.get_dataset(self.dataset_type, self.base_name)
            existing_entries = []
            for entry_name in existing_dataset.entry_names:
                try:
                    entry_data = existing_dataset.get_entry(entry_name)
                    existing_entries.append(entry_data)
                except Exception as e:
                    print(f"Warning: Could not retrieve entry {entry_name}: {e}")
            # Ensure that new_entry is appended.
            all_entries = existing_entries + [new_entry]
            timestamp = int(time.time_ns())
            new_dataset_name = f"{self.base_name}_{timestamp}"

            self.dataset = self.client.add_dataset(
                name=self.base_name,
                description=self.description,
                dataset_type=self.dataset_type,                 
                default_compute_tag="",          
                default_compute_priority=0,      
                extras={}              
            )
            for current_entry in all_entries:
                try:
                    if isinstance(current_entry, dict):
                        entry_name = current_entry["name"]
                        entry_molecule = current_entry["molecule"]
                    else:
                        entry_name = current_entry.name
                        entry_molecule = current_entry.molecule
                    self.dataset.add_entry(entry_name, entry_molecule)
                except Exception as e:
                    print(f"Warning: Error adding entry {entry_name}: {e}")
        else:

            self.dataset = self.client.add_dataset(
                name=self.base_name,
                description=self.description,
                dataset_type=self.dataset_type,                 
                default_compute_tag="",          
                default_compute_priority=0,      
                extras={}              
            )
            self.dataset.add_entry(new_entry["name"], new_entry["molecule"])
        return self.dataset

    def add_specification(self, spec_name, spec):
        if self.dataset is not None:
            self.dataset.add_specification(spec_name, spec)
        else:
            raise ValueError("Dataset not created yet.")

    def refresh_cache(self):
        if self.dataset is not None:
            self.dataset.refresh_cache()

    def get_dataset(self):
        return self.dataset

    @staticmethod
    def build_dataframe(_dataset: SinglepointDataset) -> pd.DataFrame:
        rows = []
        for rec in _dataset.iterate_records():
            try:
                entry_name, spec_name, record = rec
            except Exception:
                continue
            if record is None:
                continue
            try:
                program = record.specification.program
                driver = record.specification.driver
                method = record.specification.method
                basis = record.specification.basis
            except Exception:
                program = driver = method = basis = None
            row = {
                "spec_name": spec_name,
                "entry_name": entry_name,
                "program": program,
                "driver": driver,
                "method": method,
                "basis": basis,
                "status": record.status,
                "return_energy": record.return_result
            }
            rows.append(row)
        return pd.DataFrame(rows)

# ------------------------------
# CalculationManager Class
# ------------------------------
class CalculationManager:
    def __init__(self, client, dataset_type, dataset_name, spec_name="test_calc_spec"):
        """
        client: the ptl.PortalClient instance.
        dataset_type: e.g., 'singlepoint'
        dataset_name: the name of the dataset to work with.
        spec_name: the specification name used for calculations.
        """
        self.client = client
        self.dataset_type = dataset_type
        self.dataset_name = dataset_name
        self.spec_name = spec_name

    def submit_calculations(self):
        ds = self.client.get_dataset(self.dataset_type, self.dataset_name)
        if ds is None:
            raise ValueError("No dataset available.")
        ds.submit(specification_names=self.spec_name)

    def check_status(self):
        ds = self.client.get_dataset(self.dataset_type, self.dataset_name)
        if ds is None:
            raise ValueError("No dataset available.")
        return ds.status()

# ------------------------------
# RequestCalculationTab Class
# ------------------------------
class RequestCalculationTab:
    def __init__(self, client):
        """
        client: an instance of ptl.PortalClient.
        """
        self.client = client

    def run(self):
        st.title("Request Calculation via QCArchive")
        
        # -------------------------------------------------
        # 1) Gather SMILES from session state (for "add molecule" mode)
        # -------------------------------------------------
        mode = st.session_state.get("request_calc_mode", "update_specification")
        if mode == "add_molecule":
            st.info("Add Molecule mode: Converting provided SMILES to QCArchive objects.")
            all_smiles_list = []
            if st.session_state.get("smile_code1"):
                all_smiles_list.append(st.session_state["smile_code1"])
            if st.session_state.get("smile_code2"):
                all_smiles_list.extend(st.session_state["smile_code2"])
            if st.session_state.get("smile_code3"):
                all_smiles_list.extend(st.session_state["smile_code3"])
            if not all_smiles_list:
                st.error("No SMILES found. Please go back and select a molecule.")
                return
            get_mol_tab = GetMoleculeTab()
            qcarchive_objects = []
            for s in all_smiles_list:
                try:
                    obj = get_mol_tab.convert_smiles_qcmol(s)
                    qcarchive_objects.append((s, obj))
                except Exception as e:
                    st.error(f"Error converting SMILES '{s}': {e}")
            if not qcarchive_objects:
                st.error("No valid molecule objects could be created.")
                return
            st.session_state.qcarchive_objects = qcarchive_objects
            
            # --- REMOVAL BUTTON ---
            st.write("### Molecules Added:")
            for idx, (smiles, _) in enumerate(qcarchive_objects):
                st.write(f"{idx+1}. {smiles}")
            if st.button("Remove Added Molecules"):
                keys_to_clear = ["project_molecules", "smile_code1", "smile_code2", "smile_code3", "molecules", "qcarchive_objects"]
                for key in keys_to_clear:
                    if key in st.session_state:
                        del st.session_state[key]
                st.success("Added molecules have been removed.")
                return  # Stop further processing after removal.
            # -------------------------
        else:
            st.info("Update Specification mode: Using existing dataset; no new SMILES conversion required.")
        
        st.markdown("---")
        
        # -------------------------------------------------
        # 2) Display Calculation Specifications
        # -------------------------------------------------
        st.header("Calculation Specifications")
        st.subheader("Select QM Program, Calculation Type, etc.")
        qmprogram = st.selectbox("Program:", ('psi4', 'nwchem', 'gamess', 'other'))
        calc_type = st.selectbox("Calculation Type:", ('energy', 'gradient', 'hessian', 'optimization', 'other'))
        calc_method = st.selectbox("Theory Method:", ('HF', 'UHF', 'B3LYP', 'M062X', 'MP2', 'MP4', 'G2', 'G2MP2', 'other'))
        basis_set = st.selectbox("Basis Set:", ('3-21G', '6-31G*', '6-31+g(d,p)', '6-311+g(d,p)', 'cc-pVDZ', 'cc-pVTZ', 'other'))
        qmprogram_custom = st.text_input("Other Program") if qmprogram == 'other' else None
        calc_type_custom = st.text_input("Other Calc Type") if calc_type == 'other' else None
        calc_method_custom = st.text_input("Other calc method") if calc_method == 'other' else None
        basis_set_custom = st.text_input("Other basis set") if basis_set == 'other' else None
        final_program = qmprogram_custom if (qmprogram == 'other' and qmprogram_custom) else qmprogram
        final_calc_type = calc_type_custom if (calc_type == 'other' and calc_type_custom) else calc_type
        final_calc_method = calc_method_custom if (calc_method == 'other' and calc_method_custom) else calc_method
        final_basis_set = basis_set_custom if (basis_set == 'other' and basis_set_custom) else basis_set
        
        # -------------------------------------------------
        # 3) Specification name generation
        # -------------------------------------------------
        if "new_spec_name" not in st.session_state:
            timestamp = int(time.time_ns())
            st.session_state.new_spec_name = f"{final_program}_{final_calc_type}_{final_calc_method}_{final_basis_set}_{timestamp}"
        new_spec_name = st.session_state.new_spec_name
        
        # -------------------------------------------------
        # 4) Create or update dataset with new specification.
        # -------------------------------------------------
        if not st.session_state.get("dataset_created", False):
            qcarchive_objects = st.session_state.get("qcarchive_objects", [])
            if not qcarchive_objects:
                st.error("No valid molecule objects found. Please go back and select a valid SMILES.")
                return
            base_name = st.session_state.get("dataset_base_name")
            if base_name is None:
                # If not set, default to the selected project's dataset name.
                base_name = st.session_state.dataset_info.get("dataset_name", "default_dataset")
                st.session_state.dataset_base_name = base_name

            ds_description = st.session_state.get("dataset_description")
            if ds_description is None:
                # Optionally, default to an empty string or a description from dataset_info.
                ds_description = (
                                    f"QCArchive dataset for project '{base_name}' with specifications: "
                                    f"Program={final_program}, Calculation Type={final_calc_type}, "
                                    f"Method={final_calc_method}, Basis Set={final_basis_set}."
                                )
                st.session_state.dataset_description = ds_description

            ds_manager = DatasetManager(self.client, base_name, ds_description)
            first_obj = qcarchive_objects[0][1]
            new_entry = {"name": qcarchive_objects[0][0], "molecule": first_obj}
            dataset = ds_manager.create_or_recreate_dataset(new_entry)
            for i, (smiles, obj) in enumerate(qcarchive_objects[1:], start=1):
                entry = {"name": f"{smiles}_{i}", "molecule": obj}
                try:
                    ds_manager.dataset.add_entry(entry["name"], entry["molecule"])
                except Exception as e:
                    st.warning(f"Error adding additional molecule entry: {e}")
            spec = QCSpecification(
                program=final_program,
                driver=final_calc_type,
                method=final_calc_method,
                basis=final_basis_set
            )
            ds_manager.add_specification(new_spec_name, spec)
            st.write(f"Dataset updated: {dataset.name}")
            st.write(f"Description: {st.session_state.get('dataset_description')}")
            st.write(f"QCSpecification: {spec}")
            st.session_state.dataset_name = dataset.name
            st.session_state.dataset_type = "singlepoint"
            st.session_state.dataset_created = True
            st.session_state.dataset_base_name = st.session_state.get("dataset_base_name")
            st.session_state.dataset_description = st.session_state.get("dataset_description")
        else:
            st.info(f"Using existing dataset: {st.session_state.dataset_name}")
            dataset = self.client.get_dataset(st.session_state.dataset_type, st.session_state.dataset_name)
            spec = QCSpecification(
                program=final_program,
                driver=final_calc_type,
                method=final_calc_method,
                basis=final_basis_set
            )
            dataset.add_specification(new_spec_name, spec)
            st.write(f"New specification '{new_spec_name}' added to dataset '{dataset.name}'.")
            st.session_state.new_spec_name = new_spec_name
        
        st.markdown("---")
        # -------------------------------------------------
        # 5) Calculation Management Section with Confirmation
        # -------------------------------------------------
        st.header("Calculation Management")
        from request_calculation_tab import CalculationManager
        calc_manager = CalculationManager(
            self.client,
            st.session_state.dataset_type,
            st.session_state.dataset_name,
            st.session_state.new_spec_name
        )

        # Initialize confirmation flag if not set.
        if "confirm_submission" not in st.session_state:
            st.session_state.confirm_submission = False

        # Determine the number of jobs (for display)
        if st.session_state.get("request_calc_mode", "update_specification") == "add_molecule":
            job_count = len(st.session_state.get("qcarchive_objects", []))
        else:
            job_count = 1

        if st.session_state.confirm_submission:
            st.warning(f"You are about to submit {job_count} job{'s' if job_count != 1 else ''}.")
            col1, col2 = st.columns(2)
            if col1.button("Confirm Submission"):
                try:
                    calc_manager.submit_calculations()
                    st.success("Calculations submitted!")
                except Exception as e:
                    st.error(f"Error submitting calculations: {e}")
                st.session_state.confirm_submission = False
            if col2.button("Cancel Submission"):
                st.session_state.confirm_submission = False
        else:
            if st.button("Submit calculations"):
                st.session_state.confirm_submission = True

        
        if st.button("Check Calculation Status"):
            try:
                status_info = calc_manager.check_status()
                st.write("Calculation Status:")
                st.json(status_info)
            except Exception as e:
                st.error(f"Error retrieving status: {e}")
        
        if st.button("Back to Projects", key="back_to_projects_reqcalc"):
            st.session_state.page = "project"
            st.session_state.request_calc_mode = None
