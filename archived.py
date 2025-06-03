import json
import os
import streamlit as st

class ArchiveManager:
    def __init__(self):
        """
        Initializes the ArchiveManager using a JSON file to store archived dataset names.
        """
        if "username" in st.session_state:
            user_metadata_dir = f"./users/{st.session_state.username}/App.metadata"
            os.makedirs(user_metadata_dir, exist_ok=True)
            self.archive_file = os.path.join(user_metadata_dir, "archived_datasets.json")

        self.archived = self.load_archived_datasets()

    def load_archived_datasets(self):
        if os.path.exists(self.archive_file):
            try:
                with open(self.archive_file, "r") as f:
                    data = json.load(f)
                    return data  # expecting a list of dataset names
            except Exception:
                return []
        else:
            return []

    def save_archived_datasets(self):
        with open(self.archive_file, "w") as f:
            json.dump(self.archived, f, indent=4)

    def archive_datasets(self, dataset_names):
        for name in dataset_names:
            if name not in self.archived:
                self.archived.append(name)
        self.save_archived_datasets()

    def unarchive_dataset(self, dataset_name):
        if dataset_name in self.archived:
            self.archived.remove(dataset_name)
            self.save_archived_datasets()

    def is_archived(self, dataset_name):
        return dataset_name in self.archived

    def filter_active_datasets(self, all_datasets):
        """
        Given a list of dataset dictionaries (each with a "dataset_name" key),
        returns only those which are not archived.
        """
        return [ds for ds in all_datasets if ds.get("dataset_name") not in self.archived]


class ArchivePage:
    def __init__(self, client, archive_manager: ArchiveManager):
        self.client = client
        self.archive_manager = archive_manager

    def run(self):
        st.session_state.page = "archive"
        st.header("Archive Datasets")

        # Retrieve all datasets from the client.
        all_datasets = self.client.list_datasets()

        # Get active (non-archived) datasets
        active_datasets = self.archive_manager.filter_active_datasets(all_datasets)

        # --- Active Datasets Section ---
        with st.expander("Active Datasets (click to expand)", expanded=False):
            if active_datasets:
                for ds in active_datasets:
                    st.write(f"- {ds.get('dataset_name')} (Type: {ds.get('dataset_type')}, Records: {ds.get('record_count')})")
            else:
                st.info("No active datasets found.")

        # Multi-select active dataset names.
        active_names = [ds.get("dataset_name") for ds in active_datasets]
        selected_to_archive = st.multiselect("Select datasets to archive", options=active_names, key="archive_select")

                # Archive All Button
        if active_names:
            if st.button("Archive All Active Datasets", key="archive_all"):
                self.archive_manager.archive_datasets(active_names)
                st.success("All active datasets have been archived.")
                st.rerun()
        
        # Only show the Archive button when at least one active dataset is selected.
        if selected_to_archive:
            st.warning(f"You are about to archive {len(selected_to_archive)} dataset(s).")
            col1, col2 = st.columns(2)
            if col1.button("Confirm Archival", key="confirm_archive"):
                self.archive_manager.archive_datasets(selected_to_archive)
                st.success("Selected datasets have been archived.")
                st.rerun()
            if col2.button("Cancel", key="cancel_archive"):
                st.info("Archival cancelled.")
        
        # --- Archived Datasets Section ---
        with st.expander("Archived Datasets (click to expand)", expanded=False):
            if self.archive_manager.archived:
                for archived_ds in self.archive_manager.archived:
                    st.write(f"- {archived_ds}")
            else:
                st.info("No archived datasets.")

        # Multi-select for unarchiving.
        selected_to_unarchive = st.multiselect("Select archived dataset(s) to unarchive", options=self.archive_manager.archived, key="unarchive_select")
        if selected_to_unarchive:
            st.warning(f"You are about to unarchive {len(selected_to_unarchive)} dataset(s).")
            col1, col2 = st.columns(2)
            if col1.button("Confirm Unarchive", key="confirm_unarchive"):
                for ds_name in selected_to_unarchive:
                    self.archive_manager.unarchive_dataset(ds_name)
                st.success("Selected datasets have been unarchived.")
                st.rerun()
            if col2.button("Cancel", key="cancel_unarchive"):
                st.info("Unarchival cancelled.")

        # --- Navigation ---
        if st.button("Back to Home", key="back_from_archive"):
            st.session_state.page = "home"
            st.rerun()
