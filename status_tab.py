import streamlit as st

class StatusTab:
    def __init__(self, dataset_name, dataset_type):
        """Initialize with dataset details."""
        self.dataset = st.session_state.client.get_dataset(dataset_type, dataset_name)

    def get_status_table(self):
            """Generates and prints a Markdown table summarizing record statuses per Method/Basis Set."""
            status_data = self.dataset.status()  # Retrieve status dictionary

            # Extract all unique Methods/Basis Sets
            method_basis_sets = list(status_data.keys())

            # Create table header
            table = "| Specification | Complete | Running | Waiting | Error |\n"
            table += "|---|---|---|---|---|\n"

            total_complete = total_running = total_waiting = total_error = 0

            # Iterate over each Method/Basis Set and fetch counts
            for method in method_basis_sets:
                complete = status_data[method].get("complete", 0)
                running = status_data[method].get("running", 0)
                waiting = status_data[method].get("waiting", 0)
                error = status_data[method].get("error", 0)

                # Append row to table
                table += f"| {method} | {complete} | {running} | {waiting} | {error} | \n"
                total_complete += complete
                total_running += running
                total_waiting += waiting
                total_error += error

            table += f"| **Total** | **{total_complete}** | **{total_running}** | **{total_waiting}** | **{total_error}** |\n"

            # Print Markdown table in Streamlit
            st.markdown(table)

    def run(self):
        self.get_status_table()
