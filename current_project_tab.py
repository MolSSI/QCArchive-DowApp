import streamlit as st
import pandas as pd
import os

# current project tab
# This class reads the current project information from an external CSV file,
# displays its details, and provides a button to load the project.
class CurrentProjectTab:
    def __init__(self, username):
        """Initialize username, metadata path, and initializes load_current_project function when tab is clicked"""
        self.username = username # for retaining username, will be stored in csv file later
        #self.username = st.session_state.username 
        self.metadata_path = os.path.join(f"./users/{self.username}/App.metadata", "current_project.csv") #join path
        #self.metadata_path = os.path.join(f"./users/{st.session_state.username}/App.metadata", "current_project.csv") #join path
        self.current_project = self.load_current_project() # 

    def load_current_project(self):
        if os.path.exists(self.metadata_path): # check if current_project.csv file exists
            df = pd.read_csv(self.metadata_path)
            if not df.empty:                    # check if dataframe is empty
                #returning the second row
                return df.iloc[0].to_dict()
        
        return None
    

    
    def run(self):
        # if there is a current project
        
        if self.current_project:
            #st.write(f"**Project Name:** {self.current_project.get('dataset_name', 'N/A')}")
            #st.write(f"**Project Type:** {self.current_project.get('dataset_type', 'N/A')}")
            #st.write(f"**Record Count:** {self.current_project.get('record_count', 'N/A')}")
        
            if st.button("Current Project", help="This takes you to last Project you were working on"):
                # Get the dataset using the client stored information in session state
                dataset_type = self.current_project.get('dataset_type')
                dataset_name = self.current_project.get('dataset_name')
                st.session_state.dataset = st.session_state.client.get_dataset(dataset_type, dataset_name)
                st.session_state.page = "project"
                st.rerun()
        
        else:
            st.write("No current project set")

# This class will be initialized with the working st.session_state.dataset, and updates
# the current_project.csv with the project being viewed (i.e the datasets name, type, record and id)
# it should hold a singular record at any given time

class UpdateCurrentProject():
    def __init__(self, dataset, username):
        self.dataset = dataset
        self.username = username
        self.metadata_directory = os.path.join(f"./users/{self.username}/App.metadata")
        self.metadata_path = os.path.join(self.metadata_directory, "current_project.csv")
        if not os.path.exists(self.metadata_directory): # check if the metadata directory exists
            os.makedirs(self.metadata_directory) # if not, create it and all intermediary paths. os.makedirs is a super mkdir
                                                # that creates all intermediaries
    

    def run(self):
        "creates a dictionary of the current project, turn to a dataframe, then convert to csv file"
        #if st.session_state.dataset_info is None:
        if self.dataset is None:
            st.error("No project information available. Please select a project first.")
            return
    
        # data = {
        #    "dataset_name" : [self.dataset.dataset_name],
        #     "dataset_type" : [self.dataset.dataset_type],
        #      "record_count" : [self.dataset.record_count],
        #       "dataset_id" : [getattr(self.dataset, "dataset_id", "N/A")]
        # }
        # data = {
        #         "dataset_name": [st.session_state.dataset_info['dataset_name']],
        #         "dataset_type": [st.session_state.dataset_info['dataset_type']],
        #         "record_count": [st.session_state.dataset_info['record_count']],
        #         "dataset_id": [st.session_state.dataset_info.get('dataset_id', 'N/A')],
        # }

        data = {
           "dataset_name" : [self.dataset['dataset_name']],
            "dataset_type" : [self.dataset['dataset_type']],
             "record_count" : [self.dataset['record_count']],
              "dataset_id" : [self.dataset.get('dataset_id', 'N/A')],
        }
        df = pd.DataFrame(data)
        df.to_csv(self.metadata_path, index=False)
        st.write("Current project updated")
