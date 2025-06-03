# construct environment using environment.yml  --  may need to be modified for further QCArchive functionality
# conda env create -f environment.yml

# in terminal: streamlit run streamlit-app.py
from login_page import LoginPage
from project_search_page import ProjectSearchPage

# Classes for Project Page
from calculations_display_tab import CalculationsDisplayTab
from specifications_tab import SpecificationsTab
from status_tab import StatusTab
from request_calculation_tab import RequestCalculationTab
from request_calculation_modes import RequestCalculationAddMoleculeTab
from request_calculation_modes import RequestCalculationUpdateSpecificationTab
from request_calculation_modes import RequestCalculationOptimizationAddMoleculeTab
from request_calculation_modes import RequestCalculationOptimizationUpdateSpecificationTab
from request_calculation_modes import UnifiedRequestCalculationMenu
from get_molecule_tab import GetMoleculeTab
from add_project_page import AddProjectPage
from current_project_tab import CurrentProjectTab, UpdateCurrentProject
from home_page import HomePage
from settings_page import tooltip
from single_molecule_viewer_tab import SingleMoleculeViewerTab, OptimizationMoleculeViewerTab
from usage_metrics_tab import UsageMetricsTab
from archived import ArchiveManager, ArchivePage
import streamlit as st

st.set_page_config(page_title="QC Archive Dataset Viewer", layout="wide")

#help=tooltip("This takes you to the project you have selected")
# Adjust sidebar width
st.markdown(
    """
    <style>
    [data-testid="stSidebar"] {
        min-width: 250px;  /* Minimum width */
        max-width: 400px;  /* Maximum width */
        width: 250px;      /* Fixed width */
    }
    </style>
    """,
    unsafe_allow_html=True
)

if "client" in st.session_state and st.session_state.client is not None:
# ---- SIDEBAR ----
    st.sidebar.header("🔍 External Resources")
    st.sidebar.markdown("[📄 PubChem](https://pubchem.ncbi.nlm.nih.gov/)")
    

    st.sidebar.markdown("---")

    current_project_tab = CurrentProjectTab(st.session_state.username)
    if current_project_tab.current_project:
        st.sidebar.caption(f"**Current Project:** {current_project_tab.current_project.get('dataset_name', 'N/A')}")
        if st.sidebar.button("Current Project 📊", help="Return to your last active project"):
            dataset_type = current_project_tab.current_project.get('dataset_type')
            dataset_name = current_project_tab.current_project.get('dataset_name')
            # Load the dataset again into session
            st.session_state.dataset = st.session_state.client.get_dataset(dataset_type, dataset_name)
            st.session_state.dataset_info = {
                "dataset_name": dataset_name,
                "dataset_type": dataset_type,
                "record_count": current_project_tab.current_project.get('record_count', 'N/A'),
                "dataset_id": current_project_tab.current_project.get('dataset_id', 'N/A')
            }
            st.session_state.page = "project"
            st.rerun()
    else:
        st.sidebar.button("Current Project 📊", disabled=True, key="current_project_button_disabled") # No project → Grayed out button (disabled)

    st.sidebar.markdown("---")

    if st.sidebar.button("Reset 🔄", key="home1", help=tooltip("Clears unsaved progress, and returns to Home")):
        
        # Keys you want to keep
        keys_to_keep = {"client", "username", "page"}

        # Loop and delete everything else
        for key in list(st.session_state.keys()):
            if key not in keys_to_keep:
                del st.session_state[key]
        st.session_state.page = "home"

        st.rerun()

    if st.sidebar.button("Logout 🔙", key="login1"):
        keys_to_keep = {"page"}
        for key in list(st.session_state.keys()):
            if key not in keys_to_keep:
                del st.session_state[key]
        st.session_state.page = "login"
        
        st.rerun()

    #st.sidebar.markdown("---")
    with st.sidebar:
        st.markdown("---")
        st.markdown("### 📈 Usage Metrics")

        usage_metrics_tab = UsageMetricsTab()
        usage_metrics_tab.run()

        #st.markdown("---")

if "page" not in st.session_state:
    st.session_state.page = "login"

if st.session_state.page == "login":
    st.session_state.client = None
    st.session_state.username = None
    login_page = LoginPage()
    login_page.run()

elif st.session_state.page == "home":
    landing_page = HomePage()
    landing_page.run()

elif st.session_state.page == "archive":
    archive_manager = ArchiveManager()
    archive_page = ArchivePage(st.session_state.client, archive_manager)
    archive_page.run()
    st.stop()

elif st.session_state.page == "add_project":
    add_project_page = AddProjectPage(st.session_state.client)
    add_project_page.run()

elif st.session_state.page == "get_molecule":
    get_molecule_tab = GetMoleculeTab()
    get_molecule_tab.run()

elif st.session_state.page == "request_calculation":
    request_calculation_tab = RequestCalculationTab(st.session_state.client)
    request_calculation_tab.run()


elif st.session_state.page == "project":
    dataset = st.session_state.client.get_dataset(st.session_state.dataset_info['dataset_type'], 
                                                                   st.session_state.dataset_info['dataset_name'])
    
    print("GETTING DATASET FROM APP.PY")

    st.title(f"Project: {st.session_state.dataset_info['dataset_name']}")
    st.write(f"**Total Records:** {dataset.record_count}   |   \
            **Type:** {st.session_state.dataset_info['dataset_type']}")


    
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Display Calculations", "Display Calculations by Molecule", "Status", "Specifications", "Add Molecule", "Request Calculation"])

    with tab1:
        calculations_display_tab = CalculationsDisplayTab(st.session_state.dataset_info["dataset_name"], 
                                                          st.session_state.dataset_info["dataset_type"])
        calculations_display_tab.run()

    with tab2:
        if st.session_state["dataset_type"] == "singlepoint":
            SingleMoleculeViewerTab(st.session_state.dataset_info["dataset_name"], st.session_state.dataset_info["dataset_type"]).run()
        elif st.session_state["dataset_type"] == "optimization":
            OptimizationMoleculeViewerTab(st.session_state.dataset_info["dataset_name"], st.session_state.dataset_info["dataset_type"]).run()
        else:
            st.warning(f"Single molecule viewer not available for {st.session_state['dataset_type']} datasets")

    with tab3:
        status_tab = StatusTab(st.session_state.dataset_info["dataset_name"], 
                               st.session_state.dataset_info["dataset_type"])
        status_tab.run()

    with tab4:
        specification_tab = SpecificationsTab(st.session_state.dataset_info["dataset_name"], 
                                              st.session_state.dataset_info["dataset_type"])
        specification_tab.run()

    with tab5:
        add_page = AddProjectPage(st.session_state.client)
        add_page.run_unified_add_molecule()

    with tab6:
        req_calc_menu = UnifiedRequestCalculationMenu(st.session_state.client)
        req_calc_menu.run()



# for key, value in st.session_state.items():
#     st.write(f"**{key}**: {value}")
