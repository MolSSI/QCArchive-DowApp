import streamlit as st
import pandas as pd
import os

class UsageMetricsTab:
    def __init__(self):
        self.base_path = "./users"

    def load_user_metadata(self, username):
        """Load metadata for a specific user if the file exists."""
        metadata_path = os.path.join(self.base_path, username, "App.metadata", "dataset_list.csv")
        if os.path.exists(metadata_path):
            return pd.read_csv(metadata_path)
        return pd.DataFrame(columns=["dataset_type", "record_count"])  # Empty fallback

    def summarize_records(self, df):
        """Return a summary DataFrame showing counts per dataset_type and total."""
        if df.empty:
            return pd.DataFrame(columns=["Calculation Type", "Calculation Count"])

        summary = (
            df.groupby("dataset_type")["record_count"]
            .sum()
            .reset_index()
            .rename(columns={"dataset_type": "Calculation Type", "record_count": "Calculation Count"})
        )
        total = pd.DataFrame([["TOTAL", summary["Calculation Count"].sum()]], columns=["Calculation Type", "Calculation Count"])
        return pd.concat([summary, total], ignore_index=True)

    def get_all_users(self):
        """Get a list of usernames (directory names) under ./users"""
        return [
            name for name in os.listdir(self.base_path)
            if os.path.isdir(os.path.join(self.base_path, name)) and not name.startswith(".")
        ]

    def run(self):

        current_user = st.session_state.username

        # Load and summarize current user's data
        st.write(f"{current_user} Activity")
        user_df = self.load_user_metadata(current_user)
        user_summary = self.summarize_records(user_df)
        st.dataframe(user_summary, hide_index=True)

        # Load and summarize all users' data
        st.write("Total User Activity")

        all_users = self.get_all_users()
        combined_df = pd.concat([self.load_user_metadata(u) for u in all_users], ignore_index=True)
        combined_summary = self.summarize_records(combined_df)
        st.dataframe(combined_summary, hide_index=True)
