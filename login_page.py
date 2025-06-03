import streamlit as st
import qcportal as ptl

class LoginPage():
    def login_page(self):
        st.title("Dow Molecular Property Calculator")
        st.subheader("Login")
        username = st.text_input("Username")
        password = st.text_input("Password", type="password")

        if st.button("Login"):
            st.session_state.username = username

            # =====================================
            # MAY NEED TO CHANGE PORTAL CLIENT LINK 
            # =====================================

            st.session_state.client = ptl.PortalClient("https://dowtest.qcarchive.molssi.org/",
                                      username=username,
                                      password=password)
            st.session_state.page = "home"
            st.rerun()

    def run(self):
        self.login_page()