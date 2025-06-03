import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import json
import os

class SmilesSearchPage:
    def __init__(self, client, archive_manager=None):
        self.client = client
        self.archive_manager = archive_manager
        self.lookup_path = "./users/shared/smiles_lookup.json"

    def run(self):
        st.title("Search by SMARTS Pattern")

        smarts_input = st.text_input("Enter SMARTS pattern (e.g., C=C, c1ccccc1, N-C=O)")

        if not os.path.exists(self.lookup_path):
            st.warning("No SMILES lookup file found.")
            return

        if not smarts_input:
            return

        try:
            smarts_mol = Chem.MolFromSmarts(smarts_input)
            if smarts_mol is None:
                st.error("Invalid SMARTS pattern.")
                return
        except Exception as e:
            st.error(f"Error parsing SMARTS: {e}")
            return

        with open(self.lookup_path, "r") as f:
            smiles_lookup = json.load(f)

        hits = []
        for smiles, dataset_names in smiles_lookup.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(smarts_mol):
                hits.append((smiles, dataset_names))

        if not hits:
            st.info("No matches found.")
            return

        st.markdown(f"### Matches for SMARTS: `{smarts_input}`")

        for i, (smiles, datasets) in enumerate(hits):
            st.markdown(f"#### Match {i+1}")
            st.image(Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(300, 300)), caption=smiles)

            if self.archive_manager:
                datasets = self.archive_manager.filter_active_datasets(
                    [{"dataset_name": name} for name in datasets]
                )
                datasets = [ds["dataset_name"] for ds in datasets]

            if not datasets:
                st.info("All datasets with this SMILES are archived.")
                continue

            selected_dataset = st.selectbox(
                "Select a dataset to view",
                datasets,
                key=f"select_ds_{i}"
            )

            if st.button(f"Go to Dataset: {selected_dataset}", key=f"go_button_{i}"):
                dataset = self.client.get_dataset("singlepoint", selected_dataset)
                st.session_state.dataset_info = {
                    "dataset_name": dataset.name,
                    "dataset_type": dataset.dataset_type,
                    "record_count": dataset.record_count
                }
                # st.write(dir(dataset))
                st.session_state.dataset = dataset
                st.session_state.page = "project"
                st.session_state["metadata_loaded"] = False
                st.rerun()
