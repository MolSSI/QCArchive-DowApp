import streamlit as st

class GettingStartedPage:
    def __init__(self):
        self.title = "📘 Getting Started"
        self.subtitle = "Welcome to The DOW Molecular Property Calculator!"

    def run(self):
      st.title(self.title)
      st.markdown(f"### {self.subtitle}")

      st.write("""
      This platform is designed to help you:
      - Manage and analyze molecular project data.
      - Perform quantum chemistry calculations.
      - Visualize results and generate insights in real time.
      """)

      st.markdown("----")

      st.subheader("🧭 How to Use the App")
      THUMBNAIL_WIDTH = 250  
      gif_data = [
      {"name": "Adding a Project", "path": "gif/mp4/add_project.mp4"},
      {"name": "Requesting Calculations", "path": "gif/mp4/request_calc.mp4"},
      {"name": "Filtering", "path": "gif/mp4/filter_resized.mp4"}
      ]

      cols = st.columns(len(gif_data))

      for col, gif in zip(cols, gif_data):       
         with col:
            st.video(gif["path"])
            st.caption(gif["name"])
            #with st.expander(f"View {gif['name']}"):
                  #st.video(gif["path"])

      st.markdown("""
      1. **Project Dashboard**  
         View a summary of all projects. Quickly access recent data, results, and activity.
         
      2. **Add Project**  
         Create a new project by providing molecular structure input (SMILES, MOL files, etc.), 
         and configure calculation settings.
         
      3. **Search for Projects**  
         Use filters to find projects by name, molecule, tags, or properties.
         
      4. **Settings**  
         Toggle UI preferences like tooltips, layout, and (soon) themes.

      5. **Tooltips**  
         Hover over buttons and inputs to get quick hints—can be disabled in Settings.

      6. **Sidebar Navigation**  
         Access quick links like GitHub or PubChem.

      """)

      st.markdown("----")
      st.subheader("🛠 Supported Features")

      st.markdown("""
      - 💡 Real-time calculation results  
      - 🧪 SMILES input + 3D structure generation  
      - 📈 Visualization of properties and optimization paths  
      - 🧬 Integration with quantum chemistry backends
      """)

      st.markdown("----")
      st.info("Still have questions? Reach out via the **Contact Us** link in the sidebar!")
      


