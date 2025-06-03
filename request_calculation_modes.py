import time
import streamlit as st
import qcportal as ptl
from qcportal.singlepoint import QCSpecification
from qcportal.optimization import OptimizationSpecification
from request_calculation_tab import DatasetManager, CalculationManager
from get_molecule_tab import GetMoleculeTab
from keyword_protocol_form import keyword_protocol_form


# SINGLEPOINT: Add Molecule
class RequestCalculationAddMoleculeTab:
    def __init__(self, client):
        self.client = client

    def get_specification_input(self):
        st.header("Calculation Specifications (Add Molecule)")
        qmprogram = st.selectbox("Program:", ('psi4', 'nwchem', 'gamess', 'other'), key="rcam_program")
        calc_type = st.selectbox("Calculation Type:", ('energy', 'gradient', 'hessian', 'optimization', 'other'), key="rcam_calc_type")
        calc_method = st.selectbox("Theory Method:", 
                                   ('HF', 'UHF', 'B3LYP', 'M062X', 'MP2', 'MP4', 'G2', 'G2MP2', 'other'),
                                   key="rcam_calc_method")
        basis_set = st.selectbox("Basis Set:", 
                                 ('3-21G', '6-31G*', '6-31+g(d,p)', '6-311+g(d,p)', 'cc-pVDZ', 'cc-pVTZ', 'other'),
                                 key="rcam_basis_set")
        qmprogram_custom = st.text_input("Other Program", key="rcam_program_custom") if qmprogram == 'other' else ""
        calc_type_custom = st.text_input("Other Calc Type", key="rcam_calc_type_custom") if calc_type == 'other' else ""
        calc_method_custom = st.text_input("Other Calc Method", key="rcam_calc_method_custom") if calc_method == 'other' else ""
        basis_set_custom = st.text_input("Other Basis Set", key="rcam_basis_set_custom") if basis_set == 'other' else ""
        final_program = qmprogram_custom if (qmprogram == 'other' and qmprogram_custom) else qmprogram
        final_calc_type = calc_type_custom if (calc_type == 'other' and calc_type_custom) else calc_type
        final_calc_method = calc_method_custom if (calc_method == 'other' and calc_method_custom) else calc_method
        final_basis_set = basis_set_custom if (basis_set == 'other' and basis_set_custom) else basis_set
        # Use the keywords/protocols form for singlepoint add:
        kw_dict, p_dict = keyword_protocol_form("Enter Keywords and Protocols for Singlepoint Calculation")
        return final_program, final_calc_type, final_calc_method, final_basis_set, kw_dict, p_dict

    def run(self):
        st.title("Request Calculation - Add Molecule (Singlepoint)")

        # 0’) introduce a view-specific “have I already made this dataset?” flag
        if "sp_dataset_created" not in st.session_state:
            st.session_state.sp_dataset_created = False
        qcarchive_objects = st.session_state.get("qcarchive_objects")


        if not qcarchive_objects:
            project_mols = st.session_state.get("project_molecules", [])
            additional_mols = st.session_state.get("molecules", []) # from csv
            sdf_mol_pairs = st.session_state.get("sdfmolecules", [])  # from SDF

            qcarchive_objects = []
            for tup in project_mols:
                qcarchive_objects.append(tup)

            get_mol_tab = GetMoleculeTab()
            # Convert SMILES from CSV
            for csv_name, csv_smiles in additional_mols:
                    try:
                        qc_obj = get_mol_tab.convert_smiles_qcmol(csv_smiles)
                        qcarchive_objects.append((csv_name, qc_obj))
                    except Exception as e:
                        st.error(f"Error converting molecule '{csv_name}': {e}") 
            # Convert RDKit Mols from SDF
            for name, smiles_, rdkit_mol in sdf_mol_pairs:
                try:
                    qc_mol_ = get_mol_tab.convert_mol_qcmol(smiles_, rdkit_mol)
                    qcarchive_objects.append((name, qc_mol_))   # change name to smiles_, so SMILES appears

                except Exception as e:
                    st.error(f"Failed to convert SDF molecule '{name}': {e}")

            if qcarchive_objects:
                st.session_state.qcarchive_objects = qcarchive_objects

        if not qcarchive_objects:
            st.error("No molecules to add. Go back to the Add Molecule tab first.")
            return

        # 2) Show list + deletion UI
        with st.expander("Molecules to be added"):
            for i, (smi, _) in enumerate(qcarchive_objects, 1):
                st.write(f"{i}. {smi}")

        to_delete = st.multiselect("(Optional) Delete any of these before proceeding", [s for s,_ in qcarchive_objects])
        if st.button("(Optional) Remove Selected Molecules"):
            kept = [t for t in qcarchive_objects if t[0] not in to_delete]
            st.session_state.qcarchive_objects = kept
            st.rerun()

        # 3) Get spec inputs
        final_program, final_calc_type, final_calc_method, final_basis_set, kw_dict, p_dict = self.get_specification_input()

        # 3a) Ask user for a spec name instead of auto-timestamp/counter
        spec_name = st.text_input(
            "Name your calculation specification:",
            key="rcam_spec_name"
        )
        if not spec_name:
            st.error("❗ Please enter a name for your specification to proceed.")
            return

        # 4) Base dataset name + description
        default_base = (
            st.session_state.get("dataset_base_name")
            or st.session_state.get("dataset_info", {})
                 .get("dataset_name", f"{st.session_state.get('project_name','')}_dataset")
                 .replace(" ", "_")
        )
        if not st.session_state.sp_dataset_created:
            # only editable before creation
            base = st.text_input("Dataset name:", value=default_base, key="rcam_base_input")
            st.session_state.dataset_base_name = base
        else:
            # locked after creation
            base = st.session_state.dataset_base_name
            st.write("Dataset name:", base)

        desc = (
             f"QCArchive dataset '{base}' — Program={final_program}, "
             f"Driver={final_calc_type}, Method={final_calc_method}, "
             f"Basis={final_basis_set}"
        )
        st.session_state.dataset_description = desc

        # 5) CREATE QCArchive DATASET ONCE
        if not st.session_state.sp_dataset_created:
            st.subheader("Create QCArchive Dataset")

            if st.button("Create QCArchive Dataset", key="rcam_create"):
                dm = DatasetManager(self.client, base, desc, dataset_type="singlepoint")

                # existence check only when they click
                if dm.dataset_exists():
                    st.error(f"A dataset named '{base}' already exists.")
                else:
                    # actually build the dataset
                    ds = dm.create_or_recreate_dataset({
                        "name": qcarchive_objects[0][0],
                        "molecule": qcarchive_objects[0][1]
                    })
                    for idx, (smi, obj) in enumerate(qcarchive_objects[1:], 1):
                        ds.add_entry(name=f"{smi}_{idx}", molecule=obj)

                    # add your spec
                    dm.dataset = ds
                    dm.add_specification(
                        spec_name,
                        QCSpecification(
                            program=final_program,
                            driver=final_calc_type,
                            method=final_calc_method,
                            basis=final_basis_set,
                            keywords=kw_dict,
                            protocols=p_dict
                        )
                    )

                    # flip the one‐time flag & record the new dataset
                    st.session_state.sp_dataset_created = True
                    st.session_state.dataset            = ds
                    st.session_state.dataset_name       = ds.name
                    st.session_state.dataset_type       = "singlepoint"

                    st.success(f"✅ Dataset '{ds.name}' created!")

            # **bail out** while the flag is still false; once you click, Streamlit will
            # rerun, see sp_dataset_created=True, and skip this entire block.
            if not st.session_state.sp_dataset_created:
                st.info("Please click the button above to create your QCArchive dataset.")
                return

        #  6) RE-HYDRATE dataset 
        ds = self.client.get_dataset("singlepoint", st.session_state.dataset_name)
        if ds is None:
            st.error("Could not load dataset—please recreate it.")
            return
        st.session_state.dataset = ds

        #  7) ADD SPECIFICATION (if not already added) 
        dm = DatasetManager(self.client, base, desc, dataset_type="singlepoint")
        dm.dataset = ds
        if spec_name not in ds.specifications:
            try:
                dm.add_specification(
                    spec_name,
                    QCSpecification(
                        program=final_program,
                        driver=final_calc_type,
                        method=final_calc_method,
                        basis=final_basis_set,
                        keywords=kw_dict,
                        protocols=p_dict
                    )
                )
            except Exception as e:
                st.error(f"Error adding specification: {e}")
                return
        else:
            st.info(f"Using existing dataset: {st.session_state.dataset_name}")

        st.write(f"Dataset updated: {ds.name}")
        st.write(f"New spec: **{spec_name}**")

        #  8) CALCULATION MANAGEMENT 
        st.header("Calculation Management")
        ds_name = st.session_state.dataset_name
        ds_type = st.session_state.dataset_type
        calc_mgr = CalculationManager(self.client, "singlepoint", ds_name, spec_name)

        if "confirm_submission" not in st.session_state:
            st.session_state.confirm_submission = False

        n_jobs = len(qcarchive_objects)
        if st.session_state.confirm_submission:
            st.warning(f"Submit **{n_jobs}** job{'s' if n_jobs>1 else ''}?")
            c1, c2 = st.columns(2)
            if c1.button("Confirm", key="rcam_confirm"):
                try:
                    calc_mgr.submit_calculations()
                    st.success("Calculations submitted!")
                except Exception as e:
                    st.error(f"Submit failed: {e}")
                st.session_state.confirm_submission = False
            if c2.button("Cancel", key="rcam_cancel"):
                st.session_state.confirm_submission = False
        else:
            if st.button(f"Submit {n_jobs} calculation{'s' if n_jobs>1 else ''}", key="rcam_start"):
                st.session_state.confirm_submission = True

        if st.button("Check Calculation Status", key="rcam_calc"):
            try:
                st.json(calc_mgr.check_status())
            except Exception as e:
                st.error(f"Status check failed: {e}")


# SINGLEPOINT: Update Specification
class RequestCalculationUpdateSpecificationTab:
    def __init__(self, client):
        self.client = client
    def get_specification_input(self):
        st.header("Calculation Specifications (Update Specification)")
        qmprogram = st.selectbox(
            "Program:",
            ('psi4', 'nwchem', 'gamess', 'other'),
            key="rcus_upd_program"
        )
        calc_type = st.selectbox(
            "Calculation Type:",
            ('energy', 'gradient', 'hessian', 'optimization', 'other'),
            key="rcus_upd_calc_type"
        )
        calc_method = st.selectbox(
            "Theory Method:",
            ('HF', 'UHF', 'B3LYP', 'M062X', 'MP2', 'MP4', 'G2', 'G2MP2', 'other'),
            key="rcus_upd_calc_method"
        )
        basis_set = st.selectbox(
            "Basis Set:",
            ('3-21G', '6-31G*', '6-31+g(d,p)', '6-311+g(d,p)', 'cc-pVDZ', 'cc-pVTZ', 'other'),
            key="rcus_upd_basis_set"
        )
        qmprogram_custom = st.text_input("Other Program", key="rcus_upd_program_custom") if qmprogram == 'other' else ""
        calc_type_custom = st.text_input("Other Calc Type",    key="rcus_upd_calc_type_custom")  if calc_type  == 'other' else ""
        calc_method_custom = st.text_input("Other Calc Method", key="rcus_upd_calc_method_custom") if calc_method == 'other' else ""
        basis_set_custom = st.text_input("Other Basis Set",    key="rcus_upd_basis_set_custom")   if basis_set  == 'other' else ""
        final_program     = qmprogram_custom   if (qmprogram    == 'other' and qmprogram_custom)   else qmprogram
        final_calc_type   = calc_type_custom   if (calc_type     == 'other' and calc_type_custom)     else calc_type
        final_calc_method = calc_method_custom if (calc_method   == 'other' and calc_method_custom)   else calc_method
        final_basis_set   = basis_set_custom   if (basis_set     == 'other' and basis_set_custom)     else basis_set
        kw_dict, p_dict = keyword_protocol_form("Enter Keywords and Protocols for Singlepoint Update Specification")
        return final_program, final_calc_type, final_calc_method, final_basis_set, kw_dict, p_dict
    def run(self):

        st.title("Request Calculation – Update Specification")
        #  grab name/type from session_state or fallback to dataset_info 
        ds_info = st.session_state.get("dataset_info", {})
        ds_name = st.session_state.get("dataset_name", ds_info.get("dataset_name"))
        ds_type = st.session_state.get("dataset_type", ds_info.get("dataset_type"))
        if not ds_name or not ds_type:
            st.error("No project loaded. Please go back to Project Search and select one first.")
            return

        #  load the dataset once 
        dataset = self.client.get_dataset(ds_type, ds_name)
        if dataset is None:
            st.error(f"Could not load dataset '{ds_name}'.")
            return
        
        st.session_state.dataset = dataset
        st.session_state.dataset_name = ds_name
        st.session_state.dataset_type = ds_type
         # 1) name your new spec
        new_spec_name = st.text_input("Name your new specification:", key="rcus_spec_name")
        if not new_spec_name:
            st.error("❗ Please enter a name for your new specification.")
            return

        # 2) collect the rest of the inputs
        final_program, final_calc_type, final_calc_method, final_basis_set, kw_dict, p_dict = \
            self.get_specification_input()

        # 3) show a summary & confirmation button
        st.subheader("Review Your Specification")
        st.markdown(f"""
        - **Specification Name:** `{new_spec_name}`
        - **Program:** `{final_program}`
        - **Driver:** `{final_calc_type}`
        - **Method:** `{final_calc_method}`
        - **Basis set:** `{final_basis_set}`
        """)

        if st.button("Add Specification", key="rcus_confirm_add"):
            try:
                dataset.add_specification(new_spec_name, QCSpecification(
                    program=final_program,
                    driver=final_calc_type,
                    method=final_calc_method,
                    basis=final_basis_set,
                    keywords=kw_dict,
                    protocols=p_dict
                ))
                st.success(f"Specification '{new_spec_name}' added to dataset '{dataset.name}'.")
            except Exception as e:
                st.error(f"Error adding specification: {e}")

        st.header("Calculation Management")
        calc_manager = CalculationManager(
            self.client,
            st.session_state.dataset_type,
            st.session_state.dataset_name,
            new_spec_name
        )
        calc_manager = CalculationManager(
            self.client,
            ds_type,
            ds_name,
            new_spec_name
        )
        if "confirm_submission" not in st.session_state:
            st.session_state.confirm_submission = False

        job_count = 1
        if st.session_state.confirm_submission:
            st.warning(f"You are about to submit {job_count} job{'s' if job_count != 1 else ''}.")
            col1, col2 = st.columns(2)
            if col1.button("Confirm Submission", key="rcus_confirm_submission"):
                try:
                    calc_manager.submit_calculations()
                    st.success("Calculations submitted!")
                except Exception as e:
                    st.error(f"Error submitting calculations: {e}")
                st.session_state.confirm_submission = False
            if col2.button("Cancel Submission", key="rcus_cancel_submission"):
                st.session_state.confirm_submission = False
        else:
            if st.button("Submit calculations", key="rcus_submit_calculations"):
                st.session_state.confirm_submission = True

        if st.button("Check Calculation Status", key="rcus_check_status"):
            try:
                status_info = calc_manager.check_status()
                st.write("Calculation Status:")
                st.json(status_info)
            except Exception as e:
                st.error(f"Error retrieving status: {e}")




class RequestCalculationOptimizationAddMoleculeTab:
    def __init__(self, client):
        self.client = client

    def get_specification_input(self):
        st.header("Optimization Calculation Specifications (Add Molecule)")
        # Inner‐level QCSpecification program
        inner_program = st.text_input(
            "QCSpecification Program (inner):",
            value="psi4",
            key="rcom_opt_inner_program"
        )
        st.info("Outer‐level optimization program is fixed to 'geometric'.")
        # Outer driver is fixed
        outer_program = "geometric"
        # Inner‐level method and basis
        inner_method = st.selectbox(
            "Theory Method (inner):",
            ('HF','UHF','B3LYP','M062X','MP2','MP4','G2','G2MP2','other'),
            key="rcom_opt_inner_method"
        )
        inner_basis = st.text_input(
            "Basis Set (inner):", 
            value="3-21G",
            key="rcom_opt_inner_basis"
        )
        # Keywords & protocols form (shared)
        kw_dict, p_dict = keyword_protocol_form(
            "Enter Keywords and Protocols for Optimization Calculation"
        )
        return inner_program, outer_program, inner_method, inner_basis, kw_dict, p_dict

    def run(self):
        st.title("Optimization Calculation – Add Molecule")

        # 0) Ensure we have a dataset_created flag 
        if "opt_dataset_created" not in st.session_state:
            st.session_state.opt_dataset_created = False

        # 1) Gather or re-use qcarchive_objects 
        qcarchive_objects = st.session_state.get("qcarchive_objects", [])
        if not qcarchive_objects:
            project_mols    = st.session_state.get("project_molecules", [])
            additional_mols = st.session_state.get("molecules", [])
            sdf_mol_pairs   = st.session_state.get("sdfmolecules", [])

            # start with any project‐level SMILES
            qcarchive_objects = list(project_mols)

            get_mol_tab = GetMoleculeTab()


            # Convert SMILES from CSV
            for csv_name, csv_smiles in additional_mols:
                    try:
                        qc_obj = get_mol_tab.convert_smiles_qcmol(csv_smiles)
                        qcarchive_objects.append((csv_name, qc_obj))
                    except Exception as e:
                        st.error(f"Error converting molecule '{csv_name}': {e}")
            # SDF entries
            for name, smiles_, rdkit_mol in sdf_mol_pairs:
                try:
                    qc_mol_ = get_mol_tab.convert_mol_qcmol(smiles_, rdkit_mol)
                    qcarchive_objects.append((smiles_, qc_mol_))
                except Exception as e:
                    st.error(f"Failed to convert SDF molecule '{smiles_}': {e}")

            # persist if we found anything
            if qcarchive_objects:
                st.session_state.qcarchive_objects = qcarchive_objects

        if not qcarchive_objects:
            st.error("No molecules to add. Go back to the Add Molecule tab first.")
            return

        # 2) Show list + deletion UI
        with st.expander("Molecules to be added"):
            for i, (smi, _) in enumerate(qcarchive_objects, 1):
                st.write(f"{i}. {smi}")

        to_delete = st.multiselect(
            "(Optional) Delete any of these before proceeding",
            [s for s,_ in qcarchive_objects]
        )
        if st.button("Remove Selected Molecules", key="rcom_remove"):
            kept = [t for t in qcarchive_objects if t[0] not in to_delete]
            st.session_state.qcarchive_objects = kept
            st.rerun()

        # 3) Get specification inputs
        inner_prog, outer_prog, inner_method, inner_basis, kw_dict, p_dict = (
            self.get_specification_input()
        )
        # 3a) Ask user to name the two specs
        inner_spec_name = st.text_input(
            "Name your inner‐level spec:",
            key="rcom_opt_inner_spec_name"
        )
        outer_spec_name = st.text_input(
            "Name your outer‐level (optimization) spec:",
            key="rcom_opt_outer_spec_name"
        )
        if not inner_spec_name or not outer_spec_name:
            st.error("❗ Please enter both inner and outer spec names to proceed.")
            return
        


        # 4) base & desc (once)
        #  allow the user to pick or edit the dataset name each run 
        default_base = (
            st.session_state.get("dataset_base_name")
            or st.session_state.get("dataset_info", {})
                .get("dataset_name",
                     f"{st.session_state.get('project_name','')}_dataset")
                .replace(" ", "_")
        )
        base = st.text_input(
            "Dataset name:",
            value=default_base,
            key="rcom_opt_base_input"
        )
        st.session_state.dataset_base_name = base
    
        # build/update the description too
        desc = (
            f"QCArchive dataset '{base}' — outer={outer_prog}, "
            f"inner={inner_prog}/{inner_method}/{inner_basis}"
        )
        st.session_state.dataset_description = desc

        #5) CREATE DATASET ONCE 
        if not st.session_state.opt_dataset_created:
            st.subheader("Create QCArchive Dataset")

            if st.button("Create QCArchive Dataset", key="rcom_create"):
                dm = DatasetManager(self.client, base, desc, dataset_type="optimization")

                if dm.dataset_exists():
                    st.error(f"🚨 A dataset named '{base}' already exists.")
                    return

                ds = dm.create_or_recreate_dataset({
                    "name": qcarchive_objects[0][0],
                    "molecule": qcarchive_objects[0][1]
                })
                for idx, (smi, obj) in enumerate(qcarchive_objects[1:], 1):
                    ds.add_entry(f"{smi}_{idx}", obj)

                # flip the flag and persist the dataset
                st.session_state.opt_dataset_created = True
                st.session_state.dataset           = ds
                st.session_state.dataset_name      = ds.name
                st.session_state.dataset_type      = "optimization"
                st.success(f"Dataset '{ds.name}' created!")

                # ensures next run skips this entire block
                st.rerun()

            st.info("Please click 'Create QCArchive Dataset' to proceed.")
            return

        # 6) RE-HYDRATE dataset
        ds = self.client.get_dataset("optimization", st.session_state.dataset_name)
        if ds is None:
            st.error("Could not load dataset—please recreate it.")
            return
        st.session_state.dataset = ds
        st.session_state.dataset_type = "optimization"

        # 7) ADD SPECS — guard against duplicates
        dm = DatasetManager(self.client, base, desc, dataset_type="optimization")
        dm.dataset = ds
        if outer_spec_name not in ds.specifications:
            dm.add_specification(outer_spec_name, OptimizationSpecification(
                program=outer_prog,
                qc_specification=QCSpecification(
                    program=inner_prog,
                    driver="deferred",
                    method=inner_method,
                    basis=inner_basis,
                    keywords=kw_dict,
                    protocols={"stdout": True, "wavefunction": "all"}
                ),
                keywords=kw_dict,
                protocols=p_dict
            ))
        else:
            st.info(f"Optimization spec '{outer_spec_name}' exists; skipping.")

        st.write(f"Dataset updated: {ds.name}")
        st.write(f"Outer spec: **{outer_spec_name}**")

        # 8) CALCULATION MANAGEMENT
        st.header("Calculation Management")
        calc_mgr = CalculationManager(
            self.client,
            st.session_state["dataset_type"],
            st.session_state["dataset_name"],
            outer_spec_name
        )

        if "confirm_submission" not in st.session_state:
            st.session_state.confirm_submission = False

        n_jobs = len(qcarchive_objects)
        if st.session_state.confirm_submission:
            st.warning(f"Submit **{n_jobs}** job{'s' if n_jobs>1 else ''}?")
            c1, c2 = st.columns(2)
            if c1.button("Confirm", key='rcom_confirm'):
                try:
                    calc_mgr.submit_calculations()
                    st.success("Calculations submitted!")
                except Exception as e:
                    st.error(f"Submit failed: {e}")
                st.session_state.confirm_submission = False
            if c2.button("Cancel", key='rcom_cancel'):
                st.session_state.confirm_submission = False
        else:
            if st.button(f"Submit {n_jobs} calculation{'s' if n_jobs>1 else ''}", key='rcom_submit'):
                st.session_state.confirm_submission = True

        if st.button("Check Calculation Status", key='rcom_status'):
            try:
                st.json(calc_mgr.check_status())
            except Exception as e:
                st.error(f"Status check failed: {e}")


# OPTIMIZATION: Update Specification 
class RequestCalculationOptimizationUpdateSpecificationTab:
    def __init__(self, client):
        self.client = client

    def get_specification_input(self):
        st.header("Optimization Calculation Specifications (Update Specification)")
        qmprogram = st.selectbox("Program:", ('psi4', 'nwchem', 'gamess', 'other'), key="rcou_upd_program")
        calc_method = st.selectbox("Theory Method:", ('HF', 'UHF', 'B3LYP', 'M062X', 'MP2', 'MP4', 'G2', 'G2MP2', 'other'), key="rcou_upd_calc_method")
        basis_set = st.selectbox("Basis Set:", ('3-21G', '6-31G*', '6-31+g(d,p)', '6-311+g(d,p)', 'cc-pVDZ', 'cc-pVTZ', 'other'), key="rcou_upd_basis_set")
        qmprogram_custom = st.text_input("Other Program", key="rcou_upd_program_custom") if qmprogram == 'other' else ""
        calc_method_custom = st.text_input("Other Calc Method", key="rcou_upd_calc_method_custom") if calc_method == 'other' else ""
        basis_set_custom = st.text_input("Other Basis Set", key="rcou_upd_basis_set_custom") if basis_set == 'other' else ""
        final_program = qmprogram_custom if (qmprogram == 'other' and qmprogram_custom) else qmprogram
        final_calc_method = calc_method_custom if (calc_method == 'other' and calc_method_custom) else calc_method
        final_basis_set = basis_set_custom if (basis_set == 'other' and basis_set_custom) else basis_set
        final_driver = "optimization"
        kw_dict, p_dict = keyword_protocol_form("Enter Keywords and Protocols for Optimization Update Specification")
        return final_program, final_driver, final_calc_method, final_basis_set, kw_dict, p_dict

    def run(self):

        st.title("Request Calculation – Update Specification")
        #  grab name/type from session_state or fallback to dataset_info 
        ds_info = st.session_state.get("dataset_info", {})
        ds_name = st.session_state.get("dataset_name", ds_info.get("dataset_name"))
        ds_type = st.session_state.get("dataset_type", ds_info.get("dataset_type"))
        if not ds_name or not ds_type:
            st.error("No project loaded. Please go back to Project Search and select one first.")
            return

        #  load the dataset once 
        dataset = self.client.get_dataset(ds_type, ds_name)
        if dataset is None:
            st.error(f"Could not load dataset '{ds_name}'.")
            return
        
        st.session_state.dataset = dataset
        st.session_state.dataset_name = ds_name
        st.session_state.dataset_type = ds_type
        # 1) name your new optimization spec
        new_spec_name = st.text_input(
            "Name your new optimization specification:",
            key="rcou_spec_name"
        )
        if not new_spec_name:
            st.error("❗ Please enter a name for your new optimization specification.")
            return

        # 2) collect the rest of the inputs
        final_program, final_driver, final_calc_method, final_basis_set, kw_dict, p_dict = \
            self.get_specification_input()

        # 3) show a summary & confirmation button
        st.subheader("Review Your Optimization Specification")
        st.markdown(f"""
        - **Specification Name:** `{new_spec_name}`
        - **Program (outer):** `{final_program}`
        - **Driver:** `{final_driver}`
        - **Method (inner):** `{final_calc_method}`
        - **Basis set (inner):** `{final_basis_set}`
        """)

        if st.button("➕ Add Optimization Specification", key="rcou_confirm_add"):
            try:
                # build the nested QCSpecification for the inner-level
                inner_spec = QCSpecification(
                    program=final_program,       
                    driver="deferred",
                    method=final_calc_method,
                    basis=final_basis_set,
                    keywords=kw_dict,
                    protocols=p_dict
                )
                # now wrap in an OptimizationSpecification
                opt_spec = OptimizationSpecification(
                    program=final_driver,         # e.g. "optimization"
                    qc_specification=inner_spec,
                    keywords=kw_dict,
                    protocols=p_dict
                )
                dataset.add_specification(new_spec_name, opt_spec)
                st.success(f"✅ Optimization specification '{new_spec_name}' added to dataset '{dataset.name}'.")
            except Exception as e:
                st.error(f"Error adding specification: {e}")

        st.header("Calculation Management")

        calc_manager = CalculationManager(
            self.client,
            ds_type,
            ds_name,
            new_spec_name
        )
        if "confirm_submission" not in st.session_state:
            st.session_state.confirm_submission = False

        job_count = 1
        if st.session_state.confirm_submission:
            st.warning(f"You are about to submit {job_count} job{'s' if job_count != 1 else ''}.")
            col1, col2 = st.columns(2)
            if col1.button("Confirm Submission", key="rcou_confirm_submission"):
                try:
                    calc_manager.submit_calculations()
                    st.success("Calculations submitted!")
                except Exception as e:
                    st.error(f"Error submitting calculations: {e}")
                st.session_state.confirm_submission = False
            if col2.button("Cancel Submission", key="rcou_cancel_submission"):
                st.session_state.confirm_submission = False
        else:
            if st.button("Submit calculations", key="rcou_submit_calculations"):
                st.session_state.confirm_submission = True

        if st.button("Check Calculation Status", key="rcou_check_status"):
            try:
                status_info = calc_manager.check_status()
                st.write("Calculation Status:")
                st.json(status_info)
            except Exception as e:
                st.error(f"Error retrieving status: {e}")


class RequestCalculationAddMoleculeOnlyTab:
    """
    For existing datasets: just append new molecule(s) and re‑submit
    to one or more already‑defined specifications.
    """
    def __init__(self, client):
        self.client = client

    def run(self):
        ds_type = st.session_state.dataset_info["dataset_type"]
        ds_name = st.session_state.dataset_info["dataset_name"]
        st.session_state.dataset = self.client.get_dataset(ds_type, ds_name)
        ds = st.session_state.dataset
        st.title("Add Molecule(s) to Existing Project")
        ds = st.session_state.get("dataset")
        if ds is None:
            st.error("No project loaded. Please select one first.")
            return

        # 1) gather any molecules the user queued in any of the three upload modes
        qcarchive_objects = st.session_state.get("qcarchive_objects", [])
        if not qcarchive_objects:
            project_mols    = st.session_state.get("project_molecules", [])
            additional_mols = st.session_state.get("molecules", [])
            sdf_pairs       = st.session_state.get("sdfmolecules", [])

            get_mol_tab = GetMoleculeTab()
            qcarchive_objects = []

            # from SMILES tab
            qcarchive_objects.extend(project_mols)

            # from CSV (string SMILES) tab
            for item in additional_mols:
                if isinstance(item, tuple):
                    qcarchive_objects.append(item)
                else:
                    try:
                        obj = get_mol_tab.convert_smiles_qcmol(item)
                        qcarchive_objects.append((item, obj))
                    except Exception as e:
                        st.error(f"Error converting molecule '{item}': {e}")

            # from SDF
            for name, smiles, rdkit_mol in sdf_pairs:
                try:
                    obj = get_mol_tab.convert_mol_qcmol(smiles, rdkit_mol)
                    qcarchive_objects.append((smiles, obj))
                except Exception as e:
                    st.error(f"Failed to convert SDF molecule '{smiles}': {e}")

            # stash back
            if qcarchive_objects:
                st.session_state["qcarchive_objects"] = qcarchive_objects

        if not qcarchive_objects:
            st.error("No new molecules found. Go to the ‘Add Molecule’ tab first.")
            return

        # 2) list them
        st.write("Loaded dataset:", ds.name)
        st.write("Specs on this dataset:", list(ds.specifications.keys()))
        with st.expander("Molecules to Add", expanded=True):
            for i, (smiles, _) in enumerate(qcarchive_objects, 1):
                st.write(f"{i}. {smiles}")

        # 3) pick existing specs
        specs = list(ds.specifications.keys())
        if not specs:
            st.error("This dataset has no specifications—nothing to run against.")
            return
        choice = st.multiselect("Select spec(s) to run:", options=["All"] + specs, default=["All"])
        spec_names = specs if "All" in choice else choice
        jobs = len(qcarchive_objects) * len(spec_names)

        # 4) confirm button
        if "addmol_confirm" not in st.session_state:
            st.session_state.addmol_confirm = False

        if st.session_state.addmol_confirm:
            st.warning(f"❗You are about to submit **{jobs}** job(s).")
            c1, c2 = st.columns(2)
            if c1.button("Confirm", key="addmol_yes"):
                # append entries
                for smiles, obj in qcarchive_objects:
                    ds.add_entry(name=smiles, molecule=obj)
                ds.submit(specification_names=spec_names)
                st.success(f"Submitted {jobs} calculation(s).")
                # clear
                st.session_state.qcarchive_objects = []
                st.session_state.addmol_confirm = False
            if c2.button("Cancel", key="addmol_no"):
                st.session_state.addmol_confirm = False

        else:
            if st.button(f"Add & Submit {len(qcarchive_objects)}×{len(spec_names)} jobs",
                         key="addmol_start"):
                st.session_state.addmol_confirm = True

        # 5) status check
        st.markdown("---")
        if st.button("Check Calculation Status"):
            st.json(ds.status())

class RequestCalculationOptimizationAddMoleculeOnlyTab:
    """
    For existing optimization datasets: 
    just append new molecule(s) and re‑submit to one or more already‑defined specs.
    """
    def __init__(self, client):
        self.client = client

    def run(self):
        #  1) Load the dataset 
        ds_type = st.session_state.dataset_info["dataset_type"]
        ds_name = st.session_state.dataset_info["dataset_name"]
        ds = self.client.get_dataset(ds_type, ds_name)
        st.title("Add Molecule(s) to Existing Optimization Project")
        if ds is None:
            st.error("No project loaded. Please select one first.")
            return

        # 2) Gather any molecules the user queued in any of the three upload modes
        qcarchive_objects = st.session_state.get("qcarchive_objects", [])
        if not qcarchive_objects:
            project_mols    = st.session_state.get("project_molecules", [])
            additional_mols = st.session_state.get("molecules", [])
            sdf_pairs       = st.session_state.get("sdfmolecules", [])

            get_mol_tab = GetMoleculeTab()
            qcarchive_objects = []

            # SMILES tab
            qcarchive_objects.extend(project_mols)

            # CSV tab (string SMILES)
            for item in additional_mols:
                if isinstance(item, tuple):
                    qcarchive_objects.append(item)
                else:
                    try:
                        obj = get_mol_tab.convert_smiles_qcmol(item)
                        qcarchive_objects.append((item, obj))
                    except Exception as e:
                        st.error(f"Error converting molecule '{item}': {e}")

            # SDF tab
            for name, smiles, rdkit_mol in sdf_pairs:
                try:
                    obj = get_mol_tab.convert_mol_qcmol(smiles, rdkit_mol)
                    qcarchive_objects.append((smiles, obj))
                except Exception as e:
                    st.error(f"Failed to convert SDF molecule '{smiles}': {e}")

            if qcarchive_objects:
                st.session_state["qcarchive_objects"] = qcarchive_objects

        if not qcarchive_objects:
            st.error("No new molecules found. Go to the ‘Add Molecule’ tab first.")
            return

        # 2) Show list + deletion UI
        with st.expander("Molecules to Add", expanded=True):
            for i, (smi, _) in enumerate(qcarchive_objects, 1):
                st.write(f"{i}. {smi}")


        #  3) Show molecules & allow removal 
        st.write("Loaded dataset:", ds.name)
        st.write("Existing specs:", list(ds.specifications.keys()))
        with st.expander("Molecules to Add", expanded=True):
            for i, (smi, _) in enumerate(qcarchive_objects, 1):
                st.write(f"{i}. {smi}")

        to_delete = st.multiselect(
            "(Optional) Remove any before proceeding",
            [s for s,_ in qcarchive_objects]
        )
        if st.button("(Optional) Remove Selected Molecules", key="opt_add_remove"):
            kept = [t for t in qcarchive_objects if t[0] not in to_delete]
            st.session_state.qcarchive_objects = kept
            st.rerun()

        #  4) Pick which specs to run 
        specs = list(ds.specifications.keys())
        if not specs:
            st.error("This dataset has no specifications—nothing to run.")
            return

        choice = st.multiselect("Select spec(s) to run:",
                                options=["All"] + specs,
                                default=["All"])
        spec_names = specs if "All" in choice else choice
        job_count = len(qcarchive_objects) * len(spec_names)

        #  5) Confirmation step 
        if "opt_add_confirm" not in st.session_state:
            st.session_state.opt_add_confirm = False

        if st.session_state.opt_add_confirm:
            st.warning(f"❗About to submit **{job_count}** job(s).")
            c1, c2 = st.columns(2)
            if c1.button("Confirm", key="opt_add_yes"):
                for smi, obj in qcarchive_objects:
                    ds.add_entry(name=smi, initial_molecule=obj)
                ds.submit(specification_names=spec_names)
                st.success(f"Submitted {job_count} calculation(s).")
                st.session_state.qcarchive_objects = []
                st.session_state.opt_add_confirm = False
            if c2.button("Cancel", key="opt_add_no"):
                st.session_state.opt_add_confirm = False
        else:
            if st.button(f"Add & Submit {len(qcarchive_objects)}×{len(spec_names)} jobs",
                         key="opt_add_start"):
                st.session_state.opt_add_confirm = True

        #  6) Status check 
        st.markdown("---")
        if st.button("Check Calculation Status", key="opt_add_status"):
            st.json(ds.status())

# Unified Request Calculation Menu for Add Project (Add Molecule Only)
class UnifiedRequestCalculationAddMenu:
    def __init__(self, client):
        self.client = client

    def run(self):
        st.header("Request Calculation")
        # For Add Project, only allow add molecule modes.
        if st.session_state.get("request_calc_mode") not in {"singlepoint_add", "optimization_add"}:
            st.session_state.request_calc_mode = None
        if st.session_state.get("request_calc_mode") is None:
            st.write("Please choose one of the following calculation modes:")
            col1, col2 = st.columns(2)
            if col1.button("Singlepoint Calculation - Add Molecule", key="menu_sp_add"):
                st.session_state.request_calc_mode = "singlepoint_add"
                st.rerun()
            if col2.button("Optimization Calculation - Add Molecule", key="menu_opt_add"):
                st.session_state.request_calc_mode = "optimization_add"
                st.rerun()
            st.markdown("---")
            st.write("**Descriptions:**")
            st.write("- **Singlepoint Add Molecule:** Create a new dataset by adding molecules for singlepoint calculations.")
            st.write("- **Optimization Add Molecule:** Create a new dataset by adding molecules for optimization calculations (driver set to 'optimization').")
        else:
            mode = st.session_state.request_calc_mode
            if mode == "singlepoint_add":
                req_calc = RequestCalculationAddMoleculeTab(self.client)
                req_calc.run()
            elif mode == "optimization_add":
                req_calc = RequestCalculationOptimizationAddMoleculeTab(self.client)
                req_calc.run()
            if st.button("Back to Request Calculation Menu", key="menu_back_from_add"):
                st.session_state.request_calc_mode = None
                st.rerun()


# Unified Request Calculation Menu for Project Search (Full Menu)
class UnifiedRequestCalculationMenu:
    def __init__(self, client):
        self.client = client

    def run(self):
        # Determine project dataset type
        ds_type = st.session_state.get("dataset_type") \
                  or st.session_state.get("dataset_info", {}).get("dataset_type")

        # Initialize the request‐calc mode flag
        if "request_calc_mode" not in st.session_state:
            st.session_state.request_calc_mode = None

        st.header("Request Calculation")

        # Build the button menu dynamically based on dataset type
        if ds_type == "singlepoint":
            col1, col2 = st.columns(2)
            if col1.button("Singlepoint – Add Molecule"):
                st.session_state.request_calc_mode = "sp_add_mol"
            if col2.button("Singlepoint – Update Specification"):
                st.session_state.request_calc_mode = "sp_update_spec"
        elif ds_type == "optimization":
            col1, col2 = st.columns(2)
            if col1.button("Optimization – Add Molecule"):
                st.session_state.request_calc_mode = "opt_add_mol"
            if col2.button("Optimization – Update Specification"):
                st.session_state.request_calc_mode = "opt_update_spec"
        else:
            st.info("No project loaded. Go back to Project Search to select a dataset.")
            return

        st.markdown("---")

        # Dispatch to the selected mode
        mode = st.session_state.request_calc_mode
        if mode == "sp_add_mol":
            RequestCalculationAddMoleculeOnlyTab(self.client).run()
        elif mode == "sp_update_spec":
            RequestCalculationUpdateSpecificationTab(self.client).run()
        elif mode == "opt_add_mol":
            RequestCalculationOptimizationAddMoleculeOnlyTab(self.client).run()
        elif mode == "opt_update_spec":
            RequestCalculationOptimizationUpdateSpecificationTab(self.client).run()
        else:
            st.info("Select an action above to begin.")

        # “Back to menu” button
        if mode is not None and st.button("Back to menu"):
            st.session_state.request_calc_mode = None
            st.rerun()
