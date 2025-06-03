import streamlit as st
import pandas as pd
import os
from current_project_tab import CurrentProjectTab, UpdateCurrentProject
from settings_page import SettingsPage
from getting_started import GettingStartedPage
from about_us import AboutUsPage
from streamlit_javascript import st_javascript # works
from project_search_page import ProjectSearchPage
from add_project_page import AddProjectPage
from settings_page import tooltip
import datetime
import base64
#import streamlit.components.v1 as components

from usage_metrics_tab import UsageMetricsTab

from archived import ArchiveManager
from smiles_search_tab import SmilesSearchPage


def get_img_as_base64(path):
    with open(path, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

class HomePage:
    def __init__(self):
        self.title = f"🚀 Welcome User: {st.session_state.username}"
        self.subtitle = "Corporate Sponsorship: Dow Inc and the Molssi Team."


    def run(self):
        # ---- LOGOS TOP RIGHT ----
        logo_paths = [
            "gif/dow.png",
            "gif/molssi.png"
        ]
        logos = [get_img_as_base64(path) for path in logo_paths]
        
        html_logos = []
        for path, img in zip(logo_paths, logos):
            if "dow" in path.lower():
                html_logos.append(f'<img src="data:image/png;base64,{img}" class="dow-logo">')
            elif "molssi" in path.lower():
                html_logos.append(
                    f'<a href="https://molssi.org" target="_blank">'
                    f'<img src="data:image/png;base64,{img}" alt="MolSSI Logo">'
                    f'</a>'
        )
            else:
                html_logos.append(f'<img src="data:image/png;base64,{img}">')

        # Insert into styled HTML
        st.markdown(
            f"""
            <style>
            .logo-bar {{
                position: absolute;
                top: 10px;
                right: 10px;
                display: flex;
                gap: 10px;
                align-items: flex-start;
                z-index: 100;
            }}
            .logo-bar img {{
                height: 40px;
            }}
            .dow-logo {{
                height: 100px !important;
                transform: translateY(-28px);  /* adjust upward */
            }}
            </style>
            <div class="logo-bar">
                {''.join(html_logos)}
            </div>
            """,
            unsafe_allow_html=True
        )
        # ---- MAIN TABS ----
        tabs = st.tabs(["🏠 Overview","📘 Getting Started", "🔍 Project Search",  "SMILES Search", "➕ Begin New Project", "👥 About Us", "⚙️ Settings"])

        # Overview (landing content)
        with tabs[0]:
            st.markdown(f"### Welcome {st.session_state.username} !")
            st.write("Use the tabs above to navigate through the different tools and pages of the application.")

         # Getting Started
        with tabs[1]:
            from getting_started import GettingStartedPage  
            start = GettingStartedPage()
            start.run()

        # Project Search
        with tabs[2]:
            st.markdown("## Project Search")
            st.session_state.dataset = None
            st.session_state.coming_from_proj_select = True
            project_search_page = ProjectSearchPage(st.session_state.client)
            project_search_page.run()
            if st.button("Archive Projects", help=tooltip("Archives datasets so you don’t see it in the project search")):
                st.session_state.page = "archive"
                st.rerun()
            
            # tab1, tab2 = st.tabs(["Project Search", "Usage Metrics"])
            # with tab1:
            #     st.session_state.dataset = None
            #     st.session_state.coming_from_proj_select = True
            #     project_search_page = ProjectSearchPage(st.session_state.client)
            #     project_search_page.run()
            #     if st.button("Archive Projects", help=tooltip("Archives datasets so you don’t see it in the project search")):
            #         st.session_state.page = "archive"
            #         st.rerun()
            #     #current_proj_tab.run()
            # with tab2:
            #     usage_metrics_tab = UsageMetricsTab()
            #     usage_metrics_tab.run()

         #Archive
        with tabs[3]:
            archive_manager = ArchiveManager()
            SmilesSearchPage(st.session_state.client, archive_manager).run()

        # Begin New Project
        with tabs[4]:
            add_project_page = AddProjectPage(st.session_state.client)
            add_project_page.run()

        # About Us
        with tabs[5]:
            about = AboutUsPage()
            about.run()

        # Settings
        with tabs[6]:
            settings = SettingsPage()
            settings.run()

       


        # ---------- App metrics -----------#
        
        #user_folder = '/users/{st.session_state.username}/App.metadata/datalist_list.csv'
        def count_users(path):
            return sum(os.path.isdir(os.path.join(path, name) for name in os.listdir(path)))
        
        def count_projects(path):
            total_projects = 0
            user_project_counts = {}
            for username in os.listdir(path):
                user_path =os.path.join(path, username)
                if os.path.isdir(user_path):
                    metadata_file = os.path.join(user_path, 'App.metadata')
                    new_metadata = os.path.join(metadata_file, 'dataset_list.csv')
                    if os.path.isfile(new_metadata):
                        df = pd.read_csv(new_metadata)
                        num_projects = len(df)
                        user_project_counts[username] = num_projects
                        total_projects += num_projects
                    else:
                        user_project_counts[username] = 0 

            return total_projects
        
        total = count_projects("users/")
        

                # Footer with metric
        year = datetime.datetime.now().year
        footer = f"""
        <style>
        .footer {{
            position: fixed;
            left: 0;
            bottom: 0;
            width: 100%;
            background-color: #f0f2f6;
            color: #333333;
            text-align: center;
            padding: 10px;
            font-size: 14px;
            z-index: 100;
        }}
        .footer a {{
            color: #0066cc;
            text-decoration: none;
            margin: 0 8px;
        }}
        .footer a:hover {{
            text-decoration: underline;
        }}
        </style>
        <div class="footer">
            📑 Total Projects: {total}
            &nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;
            🐞 <a href="mailto:jwstorer@dow.com?subject=Bug%20Report">Report Bug</a>
            &nbsp;&nbsp;&nbsp;&nbsp;
            💡 <a href="mailto:jwstorer@dow.com?subject=Feature%20Request">Request Feature</a>
            &nbsp;&nbsp;&nbsp;&nbsp;
            📧 <a href="mailto:jwstorer@dow.com?subject=Help%20Needed">Help</a>
            &nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;
            © {year}
        </div>
        """
        #👤
        st.markdown(footer, unsafe_allow_html=True)

        st.sidebar.markdown("---")
        # Optional back-to-login logic
        # if st.button("Logout 🔙"):
        #     st.session_state.page = "login"
        #     st.rerun()

