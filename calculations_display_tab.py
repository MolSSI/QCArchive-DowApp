import streamlit as st
import pandas as pd
import os
from datetime import datetime
import json
import re

class CalculationsDisplayTab:
    def __init__(self, dataset_name, dataset_type):
        """Initialize with dataset details."""
        self.dataset = st.session_state.client.get_dataset(dataset_type, dataset_name)
        self.dataset_name = dataset_name
        self.record_count = self.dataset.record_count
        self.data_dir = os.path.join(f"./users/{st.session_state.username}/datasets", self.dataset_name)

    def _get_csv_path(self):
        """Returns the path to the dataset's CSV file."""
        return os.path.join(self.data_dir, "records.csv")

    def check_existing_records(self):
        """Checks if any records CSV exists in the dataset folder."""
        return any(f.startswith("records_") and f.endswith(".csv") for f in os.listdir(self.data_dir))

    def load_records(self):
        """Loads and displays the most recent dataset records."""
        records_file, _ = self.get_latest_records_file()

        if records_file:
            return pd.read_csv(records_file, dtype={"Entry Name": str})  # Return DataFrame instead of displaying it
        else:
            st.error("No records found. Please retrieve them.")
            return None

    def get_latest_records_file(self):
        """Finds the most recent records CSV in the dataset folder and extracts timestamp."""
        os.makedirs(self.data_dir, exist_ok=True)  # Ensure dataset folder exists

        files = [f for f in os.listdir(self.data_dir) if f.startswith("records_") and f.endswith(".csv")]
        
        if not files:
            return None, None  # No records file found

        # Sort by timestamp (most recent first)
        latest_file = max(files, key=lambda f: os.path.getctime(os.path.join(self.data_dir, f)))

        # Extract timestamp from filename
        timestamp_str = latest_file.replace("records_", "").replace(".csv", "")
        timestamp = datetime.strptime(timestamp_str, "%Y%m%d-%H%M%S")

        return os.path.join(self.data_dir, latest_file), timestamp
    
    def save_filters_metadata(self):
        """Extracts unique statuses, method/basis sets, and columns, then saves them to filters.json."""
        records_file, _ = self.get_latest_records_file()

        if not records_file:
            st.error("No records file found. Cannot generate filters metadata.")
            return

        df = pd.read_csv(records_file)

        filters_metadata = {
            "Unique Statuses": df["Status"].unique().tolist(),
            "Unique Method/Basis Sets": df["Method/Basis Set"].unique().tolist(),
            "All Columns": df.columns.tolist()
        }

        filters_path = os.path.join(self.data_dir, "filters.json")

        with open(filters_path, "w") as f:
            json.dump(filters_metadata, f, indent=4)

        #st.success(f"Filters metadata saved to {filters_path}")

    def save_records(self):
        """Retrieves records from QCArchive and saves them as a CSV file with a progress bar."""
        os.makedirs(self.data_dir, exist_ok=True)  # Create data folder if needed
        
        smiles_lookup_path = "./users/shared/smiles_lookup.json"
        os.makedirs(os.path.dirname(smiles_lookup_path), exist_ok=True)
        if os.path.exists(smiles_lookup_path):
            with open(smiles_lookup_path, "r") as f:
                smiles_lookup = json.load(f)
        else:
            smiles_lookup = {}
        
        records_list = []

        progress_placeholder = st.empty()  # Placeholder for progress bar
        status_placeholder = st.empty()  # Placeholder for counter display
        progress_bar = progress_placeholder.progress(0)  # Initialize progress bar
        counter = 0  # Initialize counter


        # First loop: Extract molecular formula & weight from dataset entries
        entry_molecule_data = {}
        for entry in self.dataset.iterate_entries():
            entry_name = entry.name  # Extract entry name
            if hasattr(entry, "molecule"):
                entry_molecule_data[entry_name] = {
                    "Molecular Formula": entry.molecule.get_molecular_formula(),
                    "Molecular Weight": entry.molecule.molecular_weight()
                }

        # Second loop: Extract other record details
        for entry_name, spec_name, record in self.dataset.iterate_records(include=["molecule"]):
            if hasattr(record, "molecule") and record.molecule:
                smiles = record.molecule.extras.get("smiles", "")
            elif hasattr(record, "final_molecule") and record.final_molecule:
                smiles = record.final_molecule.extras.get("smiles", "")
            else:
                smiles = ""
            record_data = {
                "ID": record.id,
                "Entry Name": str(entry_name),
                "Status": record.status,
                "SMILES": smiles
            }

            if entry_molecule_data:
                record_data["Molecular Formula"] = entry_molecule_data.get(entry_name, {}).get("Molecular Formula", "N/A")
                record_data["Molecular Weight"] = entry_molecule_data.get(entry_name, {}).get("Molecular Weight", "N/A")
                
            record_data["Method/Basis Set"] = spec_name
            record_data["Author"] = record.owner_user
            record_data["Properties"] = record.properties
        
                        
            records_list.append(record_data)

            smiles = record.molecule.extras.get("smiles", "") if hasattr(record, "molecule") else ""
            record_data["SMILES"] = smiles

            # Track dataset in smiles lookup
            if smiles and isinstance(smiles, str):
                if smiles not in smiles_lookup:
                    smiles_lookup[smiles] = []
                if self.dataset_name not in smiles_lookup[smiles]:
                    smiles_lookup[smiles].append(self.dataset_name)


            # Update progress and counter display
            counter += 1
            progress_bar.progress(int((counter / self.record_count) * 100))
            status_placeholder.write(f"Processing record {counter} / {self.record_count}")

        # Save updated smiles lookup table
        with open(smiles_lookup_path, "w") as f:
            json.dump(smiles_lookup, f, indent=4)

        # Convert records to DataFrame and flatten properties
        df = pd.DataFrame(records_list)
        df = self._flatten_properties(df)
        df["Entry Name"] = [rec["Entry Name"] for rec in records_list]
        
        # Clean up old CSVs
        for file in os.listdir(self.data_dir):
            if file.startswith("records_") and file.endswith(".csv"):
                os.remove(os.path.join(self.data_dir, file))
                
        # Generate timestamped filename
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        records_filename = f"records_{timestamp}.csv"
        records_path = os.path.join(self.data_dir, records_filename)

        df.to_csv(records_path, index=False)

            # Log update in log.csv
        log_path = os.path.join(self.data_dir, "log.csv")
        log_exists = os.path.exists(log_path)

        with open(log_path, "a") as log_file:
            if not log_exists:
                log_file.write("Timestamp,Records File,Record Count\n")  # Write header if new log
            log_file.write(f"{timestamp},{records_filename},{self.record_count}\n")


        # Clear progress bar and counter once completed
        progress_placeholder.empty()
        status_placeholder.empty()

        return records_path

    @staticmethod
    def _flatten_properties(df):
        """Expands the properties dictionary into separate columns."""
        properties_df = df["Properties"].apply(pd.Series)
        return pd.concat([df.drop(columns=["Properties"]), properties_df], axis=1)
    
    @staticmethod
    def parse_formula(formula):
        """Convert a molecular formula string into a dictionary of element counts."""
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        return {el: int(count) if count else 1 for el, count in elements}
    
    @staticmethod
    def parse_conditions(condition_str):
        """Parse a condition string like 'C>3H>5O<1' into a list of (element, operator, value) tuples."""
        conditions = re.findall(r'([A-Z][a-z]*)([><=]?)(\d+)', condition_str)
        parsed_conditions = []
        
        for el, op, value in conditions:
            if not op:  # If no operator is provided, assume '='
                op = '='
            if not value:  # If no value is provided, assume 1 (e.g., "C" -> "C=1")
                value = '1'
            parsed_conditions.append((el, op, int(value)))

        return parsed_conditions
    
    def _compute_checkbox_columns(self, options, num_columns=4):
        """Precomputes checkbox column distribution for efficiency."""
        num_rows = -(-len(options) // num_columns)  # Ceiling division
        return [options[i * num_rows:(i + 1) * num_rows] for i in range(num_columns)]
    
    def filter_molecules(self, df, condition_str):
        """Filter the DataFrame based on molecular formula constraints."""
        conditions = self.parse_conditions(condition_str)

        def check_formula(formula):
            parsed = self.parse_formula(formula)
            for el, op, value in conditions:
                el_count = parsed.get(el, 0)  # Default to 0 if element is absent
                if op == ">" and not (el_count > value):
                    return False
                elif op == "<" and not (el_count < value):
                    return False
                elif op == "=" and not (el_count == value):
                    return False
            return True

        return df[df['Molecular Formula'].apply(check_formula)]
    
    @staticmethod
    def filter_dataframe(df, columns, conditions):
        try:
            for col, cond in zip(columns, conditions):
                query_string = f"`{col}` {cond}"
                df = df.query(query_string)
            return df
        except Exception as e:
            st.error(f"Invalid condition: {e}")
            return df
        
    @staticmethod
    def get_column_types(df: pd.DataFrame) -> dict:
        column_types = {}
        for column in df.columns:
            non_null_values = df[column].dropna()
            if not non_null_values.empty:
                inferred_type = type(non_null_values.iloc[0])
                column_types[column] = inferred_type.__name__
            else:
                column_types[column] = "Unknown"  # or None, or keep as empty string
        return column_types


    def filter_records(self):
        """Displays filtering UI with checkboxes arranged in columns and applies filters to the dataset records."""
        filters_path = os.path.join(self.data_dir, "filters.json")

        if not os.path.exists(filters_path):
            # st.error("Dataset not found locally. Please retrieve records from QC Archive.")
            return

        # Load filters from JSON file
        with open(filters_path, "r") as f:
            filters_metadata = json.load(f)

        st.subheader("Filters")

        if st.session_state.coming_from_proj_select == True:
            if "selected_status" in st.session_state:
                del st.session_state["selected_status"] 
            if "selected_method" in st.session_state:
                del st.session_state["selected_method"]
            if "selected_column" in st.session_state:
                del st.session_state["selected_column"]
            st.session_state.conditional_columns = []
            st.session_state.column_conditions = []
            df = self.load_records()
            if df is None:
                return
            else:
                st.session_state.column_type_dict = self.get_column_types(df)
            st.session_state.coming_from_proj_select = False


        def create_checkbox_filter(label, options, key_prefix, num_columns=1):
            """Creates checkboxes for filtering with Select All and Clear All buttons, optionally divided into columns."""
            selected_options = st.session_state.get(f"selected_{key_prefix}", options.copy())

            col1, col2 = st.columns([0.5, 0.5])
            if col1.button(f"Select All {label}"):
                selected_options = options.copy()
            if col2.button(f"Clear {label}"):
                selected_options = []

            if num_columns > 1:
                # Split options into columns
                num_rows = -(-len(options) // num_columns)  # Ceiling division
                columns = st.columns(num_columns)
                checkboxes = []
                for col_idx in range(num_columns):
                    with columns[col_idx]:
                        for item in options[col_idx * num_rows : (col_idx + 1) * num_rows]:
                            checked = item in selected_options
                            if st.checkbox(item, value=checked, key=f"{key_prefix}_{item}"):
                                checkboxes.append(item)
            else:
                checkboxes = [
                    item for item in options
                    if st.checkbox(item, value=(item in selected_options), key=f"{key_prefix}_{item}")
                ]

            st.session_state[f"selected_{key_prefix}"] = checkboxes
            return checkboxes

        # Status Filter
        with st.expander("Status", expanded=False):
            selected_statuses = create_checkbox_filter("Statuses", filters_metadata["Unique Statuses"], "status")

        # Methods/Basis Sets Filter
        with st.expander("Methods/Basis Sets", expanded=False):
            selected_methods = create_checkbox_filter("Methods/Basis Sets", filters_metadata["Unique Method/Basis Sets"], "method", num_columns=4)

        # Column Selection
        with st.expander("Columns", expanded=False):
            selected_columns = create_checkbox_filter("Columns", filters_metadata["All Columns"], "column", num_columns=4)

        with st.expander("Column conditionals", expanded=False):
            selected_column = st.selectbox("Select a column to filter", filters_metadata["All Columns"])
            column_type = st.session_state.column_type_dict[selected_column]
            # st.write(st.session_state.column_type_dict)
            # st.write(selected_column)
            st.markdown("""
                <style>
                div[data-testid="stForm"] {
                border: 0px solid #444;
                border-radius: 8px;
                padding: 10px;
                box-shadow: none;
                }
                </style>
            """, unsafe_allow_html=True)
            with st.form("add_filter_form", clear_on_submit=False):
                # selected_column = st.selectbox("Select a column to filter", filters_metadata["All Columns"])
                column_type = st.session_state.column_type_dict[selected_column]
                # st.write(st.session_state.column_type_dict)
                # st.write(selected_column)
                # st.write(column_type)

                # Friendly type description + suggestion
                if column_type == "int64":
                    readable_type = "integers"
                    suggestion = "`> 1`, `< 1`, `!= 1` or `== 1`"
                elif column_type == "float64":
                    readable_type = "floating point numbers"
                    suggestion = "`> 0.0` or `< 0.0`"
                elif column_type == "str":
                    readable_type = "text values"
                    suggestion = '`== "sample_text"` or `!= "sample_text"`'
                elif column_type == "bool":
                    readable_type = "True/False values"
                    suggestion = "`== True` or `== False`"
                else:
                    readable_type = f"{column_type} values"
                    suggestion = ""

                if suggestion:
                    condition = st.text_input(f"Selected column contains **{readable_type}** (try {suggestion})")
                else:
                    condition = st.text_input("Selected column contains an unknown value type - filters may not work")

                valid_condition = True
                error_msg = ""

                # Validation based on type
                if column_type == "float64":
                    if "==" in condition:
                        valid_condition = False
                        error_msg = "Avoid using `==` for floating point numbers due to precision issues."
                    elif not re.match(r"^\s*[<>!]=?\s*-?\d+(\.\d+)?\s*$", condition):
                        valid_condition = False
                        error_msg = "Use a valid float comparison like `> 3.5` or `!= 0.0`."

                elif column_type == "int64":
                    if not re.match(r"^\s*[<>!=]=?\s*-?\d+\s*$", condition):
                        valid_condition = False
                        error_msg = "Use a valid integer comparison like `== 1`, `> 5`, or `!= 3`."

                elif column_type == "str":
                    if not re.match(r'^\s*(==|!=)\s*".+"\s*$', condition):
                        valid_condition = False
                        error_msg = 'Use the format `== "text"` or `!= "text"` (double quotes required).'

                elif column_type == "bool":
                    if condition.strip() not in ["== True", "== False"]:
                        valid_condition = False
                        error_msg = "Condition must be `== True` or `== False`."

                else:
                    valid_condition = False
                    error_msg = "Unsupported column type — filtering may not work."

                if condition and not valid_condition:
                    st.warning(f"Invalid condition: {error_msg}")

                submitted = st.form_submit_button("\+ Add Filter")
                if submitted and selected_column and condition and valid_condition:
                    st.session_state.conditional_columns.append(selected_column)
                    st.session_state.column_conditions.append(condition)
                elif submitted and not valid_condition:
                    st.warning("Fix condition before adding.")
            
            # --- Display Applied Filters ---
            st.markdown("### Applied Filters")
            if st.session_state.conditional_columns:
                for i, (col, cond) in enumerate(zip(st.session_state.conditional_columns, st.session_state.column_conditions)):
                    tag = f'<span style="background-color:#1f77b4; color:white; padding:4px 10px; margin-right:6px; border-radius:8px; display:inline-block;">{col} {cond}</span>'
                    col1, col2 = st.columns([5, 1])
                    with col1:
                        st.markdown(tag, unsafe_allow_html=True)
                    with col2:
                        if st.button("X", key=f"del_{i}"):
                            st.session_state.conditional_columns.pop(i)
                            st.session_state.column_conditions.pop(i)
                            st.rerun()
            else:
                st.write("None")

            # Delete all button
            if st.session_state.conditional_columns:
                if st.button("X Delete All Filters"):
                    st.session_state.conditional_columns = []
                    st.session_state.column_conditions = []
                    st.rerun()

        
        if self.dataset.dataset_type == 'singlepoint':
            with st.expander("Molecular formula (singlepoint only)", expanded=False):
                molecule_condition_string = st.text_input("Enter molecular formula")
                # st.markdown("C1H3 (or C=1H=3) Will show all molecules with exactly 1 carbon and 3 hydrogens (and any number of other atoms)  \n \
                #              C>2  Will show all molecules with more than 2 carbons  \n \
                #              C>2 C<5 Will show all molecules with more than 3 and less than 5 carbons  \n \
                #              N0 or N=0 will show all molecules with no nitrogen  \n \
                #             ")
        else:
            molecule_condition_string=""

        

        st.markdown("---")

        # Apply filters only when the button is clicked
        if st.button("Show Records w/ Filters"):
            df = self.load_records()
            if df is None:
                return

            # Apply filtering
            if selected_statuses:
                df = df[df["Status"].isin(selected_statuses)]
            if selected_methods:
                df = df[df["Method/Basis Set"].isin(selected_methods)]
            if molecule_condition_string:
                df = self.filter_molecules(df, molecule_condition_string.replace(" ", ""))
            if selected_column:
                df = self.filter_dataframe(df, st.session_state.conditional_columns, st.session_state.column_conditions)
            if selected_columns:
                df = df[selected_columns]

            # Display the filtered DataFrame
            st.dataframe(df, hide_index=True)

            if df.empty:
                st.warning("Table is empty, try fewer filters")

            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("Download CSV", csv, f"{self.dataset_name}.csv", "text/csv")

    def get_status_table(self):
        """Generates and prints a Markdown table summarizing record statuses per Method/Basis Set."""
        status_data = self.dataset.status()  # Retrieve status dictionary

        # Extract all unique Methods/Basis Sets
        method_basis_sets = list(status_data.keys())

        # Create table header
        table = "| Method/Basis Set | Complete | Waiting | Error |\n"
        table += "|---|---|---|---|\n"

        # Iterate over each Method/Basis Set and fetch counts
        for method in method_basis_sets:
            complete = status_data[method].get("complete", 0)
            waiting = status_data[method].get("waiting", 0)
            error = status_data[method].get("error", 0)

            # Append row to table
            table += f"| {method} | {complete} | {waiting} | {error} |\n"

        # Print Markdown table in Streamlit
        st.markdown("## Status Overview by Method/Basis Set")
        st.markdown(table)

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
    
    def compare_local_remote(self, records_file):
        # Compare local vs remote record counts
        local_record_count = 0
        if records_file:
            try:
                local_df = pd.read_csv(records_file, dtype={"Entry Name": str})
                local_record_count = len(local_df)

                # NEW: Check for running/waiting statuses
                if "Status" in local_df.columns:
                    incomplete_mask = local_df["Status"].str.lower().isin(["waiting", "running"])
                    num_incomplete = incomplete_mask.sum()

                    if num_incomplete > 0:
                        st.warning(
                            f"⚠️ Your downloaded data includes {num_incomplete} records that are still **running** or **waiting**. "
                            "You may want to check the **Status tab** for updates and **redownload** the dataset later to see final results."
                        )
            except Exception as e:
                st.warning(f"Error reading local records file: {e}")

        remote_record_count = self.record_count

        if local_record_count < remote_record_count:
            st.warning(f"""⚠️ The local copy has {local_record_count} records, but the QCArchive server has {remote_record_count} records.  
        You may want to click **Retrieve Records from QCArchive** below to update.""")


    def run(self):
        # Navigation button to go back to dataset selection
        records_file, last_updated = self.get_latest_records_file()

        if records_file:
            st.write(f"**Last Updated:** {last_updated.strftime('%Y-%m-%d %H:%M:%S')}")  # Display last update time
            self.compare_local_remote(records_file)
        else:
            st.write("No records found locally. Please retrieve them from QC Archive:")


        if st.button("Retrieve Records from QCArchive"):
            self.save_records()
            self.save_filters_metadata()  # New function to store filter metadata
            st.rerun()

        self.filter_records()  # Display filtering UI
