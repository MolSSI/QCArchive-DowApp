import streamlit as st
import pandas as pd

class SpecificationsTab:
    def __init__(self, dataset_name, dataset_type):
        """Initialize with dataset details."""
        self.dataset = st.session_state.client.get_dataset(dataset_type, dataset_name)

    def flatten_dict(self, d, parent_key=""):
        """
        Recursively flattens a nested dictionary.
        Converts {"specification": {"program": "psi4"}} into {"specification_program": "psi4"}
        """
        items = []
        for k, v in d.items():
            new_key = f"{parent_key}_{k}" if parent_key else k  # Create hierarchical keys
            if isinstance(v, dict):  # Recursively flatten dictionaries
                items.extend(self.flatten_dict(v, new_key).items())
            else:
                items.append((new_key, v))
        return dict(items)

    def display_specifications(self):
        """Displays dataset specifications in a structured table."""

        # Extract specifications
        specifications = self.dataset.specifications

        spec_list = []
        for name, details in specifications.items():
            details_dict = details.dict()  # Convert Pydantic model to dictionary
            flat_details = self.flatten_dict(details_dict)  # Flatten nested dictionary

            # Ensure "Specification Name" is included
            flat_details["Specification Name"] = name
            spec_list.append(flat_details)

        # Convert to DataFrame
        spec_df = pd.DataFrame(spec_list)

        # Display table
        st.dataframe(spec_df)

        csv = spec_df.to_csv(index=False).encode('utf-8')
        st.download_button("Download CSV", csv, f"{self.dataset.name}_specifications.csv", "text/csv")

    def run(self):
        self.display_specifications()
