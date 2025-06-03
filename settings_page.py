import streamlit as st
import pandas as pd
import os


def tooltip(text):
    return text if st.session_state.get("show_tooltips", True) else None


# Allows user to choose their App Preferences
class SettingsPage:
    def __init__(self):
        self.initialize_session_defaults()

    def initialize_session_defaults(self):
        if "show_tooltips" not in st.session_state:
            st.session_state.show_tooltips = True

        # if "dark_mode" not in st.session_state:
        #     st.session_state.dark_mode = True
    
    # def apply_theme(self):
    #     if st.session_state.dark_mode:
    #         dark_css = """
    #         <style>
    #         .stApp {
    #             background-color: #0e1117;
    #             color: white;
    #         }
    #         body {
    #             background-color: #0e1117;
    #             color: white;
    #         }
    #         .css-10trblm, .css-1v3fvcr, .stButton>button, .stMarkdown, .stTextInput>div>div>input,
    #         .stSelectbox>div>div>div>div, .stSidebar, .css-1d391kg {
    #             color: white !important;
    #         }
    #         </style>
    #     """
    #         st.markdown(dark_css, unsafe_allow_html=True)
    #     else:
    #         light_css = """
    #         <style>
    #         .stApp {
    #             background-color: white;
    #             color: #1e1e1e;
    #         }
    #         body {
    #             background-color: white;
    #             color: #1e1e1e;
    #         }
    #         h1, h2, h3, h4, h5, h6, p, div {
    #             color: #1e1e1e !important;
    #             }
    #         .css-10trblm, .css-1v3fvcr, .stButton>button, .stMarkdown, .stTextInput>div>div>input,
    #         .stSelectbox>div>div>div>div, .stSidebar, .css-1d391kg {
    #             color: #1e1e1e !important;
    #         }
    #         </style>
    #     """
    #         st.markdown(light_css, unsafe_allow_html=True)

    def run(self):
        #self.apply_theme() #apply the current theme
        st.title("⚙️ Settings & Configuration")
        st.markdown("Choose App preferences")

        #----App display and behavior -------

        st.subheader("App Behavior")
        st.session_state.show_tooltips = st.checkbox(
            "Show tooltips", value = st.session_state.show_tooltips
        )

        # st.subheader("Display Options")
        # st.session_state.dark_mode = st.checkbox(
        #     "Enable dark mode",
        #     value=st.session_state.dark_mode
        # )

        st.markdown("---")

        if st.button("✅ Save Settings"):
            st.success("Settings saved! These will be active during this session.")

        st.markdown("#### Current Session Settings")
        st.json({
            "Show Tooltips": st.session_state.show_tooltips
            # "Dark Mode": st.session_state.dark_mode,
        })