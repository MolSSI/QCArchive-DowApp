import streamlit as st

class AboutUsPage:
    def __init__(self):
        self.title = "👥 About Us"
        self.team_name = "The Molecular Science and Software Engineering Capstone"
        self.institution = "University of California, Berkeley"
        self.contact_email = "kofi_mireku@berkeley.edu"
        self.contact_email_2 = "brandon_ton@berkeley.edu"
        self.contact_email_3 = "bkrohn-hansen@berkeley.edu"
        self.github_url = "https://github.com/MirekuKofi"
        self.funding_note = "This project was sponsored by Dow Inc, Molssi and the UC Berkeley College of Chemistry"

    def run(self):
        st.title(self.title)

        st.markdown(f"### 🧬 {self.team_name}")
        st.markdown(f"**Institution:** {self.institution}")
        st.markdown(f"**Sponsorship:** {self.funding_note}")
        st.markdown("---")
 
        st.subheader("🔍 Mission")
        st.write("""
        We aim to make quantum chemistry and molecular property prediction more accessible to students,
        researchers, and professionals. Our tools bridge user-friendly interfaces with powerful computational chemistry engines. Computational chemistry 
        capability was made possible through QCArchive.
        """)

        st.subheader("🛠 Tech Stack")
        st.write("""
        - Python + Streamlit for UI
        - RDKit for molecular processing
        - Quantum chemistry backends (QCArchive)
        """)

        st.subheader("🤝 Contributors")
        st.markdown("""
        - **Kofi Mireku** – Lead Developer  
        - **Brandon Ton** – Lead Developer 
        - **Ben Krohn-Hansen** – Lead Developer
        - **Dr. Storer** – Project Sponsor
        - **Dr. Benjamin Pritchard** – Faculty Guide
        - **Dr. Jessica Nash** – Faculty Guide 
        """)

        st.subheader("📫 Contact")
        st.markdown(f"Email: [{self.contact_email}](mailto:{self.contact_email})")
        st.markdown(f"Email: [{self.contact_email_2}](mailto:{self.contact_email_2})")
        st.markdown(f"Email: [{self.contact_email_3}](mailto:{self.contact_email_3})")
        st.markdown(f"GitHub: [Visit Repository]({self.github_url})")

        st.markdown("---")
        st.success("Thank you for using the DOW Molecular Property Calculator!")
