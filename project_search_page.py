import streamlit as st
import qcportal as ptl
import pandas as pd
import os
from settings_page import tooltip
from archived import ArchiveManager  # Adjust this import based on where your ArchiveManager class is

class ProjectSearchPage:
    def __init__(self, client, archive_manager=None):
        """Initialize the selector with a QCArchive client."""
        self.client = client
        self.archive_manager = archive_manager or ArchiveManager()
        self.metadata_file = os.path.join(f"./users/{st.session_state.username}/App.metadata", "dataset_list.csv")


    def initialize_metadata(self):
        """Fetch dataset metadata only once per session and store it in a CSV file."""
        status_box = st.empty()  # Placeholder for status message
        status_box.write("Fetching dataset metadata from QCArchive...")

        if not os.path.exists(f"./users/{st.session_state.username}/App.metadata"):
            os.makedirs(f"./users/{st.session_state.username}/App.metadata")  # Ensure directory exists

        dataset_list = self.client.list_datasets()
        df = pd.DataFrame(dataset_list)
        df.to_csv(self.metadata_file, index=False)  # Save to CSV

        user_dataset_dir = f"./users/{st.session_state.username}/datasets"

        df["downloaded"] = df["dataset_name"].apply(
            lambda name: os.path.exists(os.path.join(user_dataset_dir, name))
        )

        # Mark metadata as loaded
        st.session_state.metadata_loaded = True

        # Clear the status message
        status_box.empty()

    def update_recent_list(self, dataset_name):
        recent_path = os.path.join(f"./users/{st.session_state.username}/App.metadata", "recent_datasets.txt")

        recent = []
        if os.path.exists(recent_path):
            with open(recent_path, "r") as f:
                recent = f.read().splitlines()

        if dataset_name in recent:
            recent.remove(dataset_name)

        recent.insert(0, dataset_name)

        with open(recent_path, "w") as f:
            for name in recent:
                f.write(name + "\n")


    def run(self):
        """Displays dataset selection UI and dataset details using cached metadata."""
        # Check if metadata has already been loaded this session
        # if "metadata_loaded" not in st.session_state or not st.session_state.metadata_loaded:
        self.initialize_metadata()

        # Load dataset metadata from CSV
        if os.path.exists(self.metadata_file):
            df = pd.read_csv(self.metadata_file)
        else:
            st.error("Missing metadata. Please restart the app.")
            return

        # Convert to list of dicts and filter out archived datasets
        dataset_list = df.to_dict(orient="records")
        active_datasets = self.archive_manager.filter_active_datasets(dataset_list)

        if not active_datasets:
            st.info("All datasets are archived.")
            return
        
        with st.expander("Filters", expanded=False):
            filter_type = st.selectbox("Filter by dataset type", ["all", "singlepoint", "optimization"], key="filter_ds_type")

        # Reconstruct options list using filtered datasets
        if filter_type != "all":
            active_datasets = [ds for ds in active_datasets if ds["dataset_type"] == filter_type]

        # Prioritize by recency
        recent_path = os.path.join(f"./users/{st.session_state.username}/App.metadata", "recent_datasets.txt")
        recent = []
        if os.path.exists(recent_path):
            with open(recent_path, "r") as f:
                recent = f.read().splitlines()

        def sort_key(ds):
            return recent.index(ds["dataset_name"]) if ds["dataset_name"] in recent else float('inf')

        active_datasets.sort(key=sort_key)
        dataset_options = [f"{ds['dataset_name']} (Records: {ds['record_count']}  Type: {ds['dataset_type']})" for ds in active_datasets]

        selected_option = st.selectbox("Select Project:", dataset_options)

        if selected_option:
            selected_index = dataset_options.index(selected_option)
            selected_dataset_info = active_datasets[selected_index]

            st.write(f"**Type:** {selected_dataset_info['dataset_type']}")
            st.write(f"**Record Count:** {selected_dataset_info['record_count']}")

            st.session_state.dataset_info = selected_dataset_info
            st.session_state.dataset_type = selected_dataset_info.get("dataset_type", "singlepoint")

            if st.button("Go to Project Page",  help=tooltip("This takes you to the project you have selected")):
                self.update_recent_list(selected_dataset_info["dataset_name"])
                st.session_state.dataset_info = selected_dataset_info
                st.session_state.dataset = st.session_state.client.get_dataset(
                    selected_dataset_info["dataset_type"],
                    selected_dataset_info["dataset_name"]
                )
                from current_project_tab import UpdateCurrentProject
                updater = UpdateCurrentProject(selected_dataset_info, st.session_state.username)
                updater.run()

                st.session_state.page = "project"
                st.session_state["metadata_loaded"] = False
                st.rerun()

