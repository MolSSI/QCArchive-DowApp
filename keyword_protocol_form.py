import streamlit as st

def auto_convert(val):
    try:
        return int(val)
    except:
        try:
            return float(val)
        except:
            return val

def keyword_protocol_form(form_label="Enter Keywords and Protocols"):
    st.subheader(form_label)

    # Initialize storage lists if not present.
    if "keywords" not in st.session_state:
        st.session_state.keywords = []  # List of tuples: (key, value)
    if "protocols" not in st.session_state:
        st.session_state.protocols = []  # List of tuples: (key, value)

    # --- Keywords Section ---
    st.markdown("#### Keywords")
    with st.expander("Existing Keywords (click to expand)", expanded=False):
        if st.session_state.keywords:
            for idx, (k, v) in enumerate(st.session_state.keywords, start=1):
                st.write(f"{idx}. {k}: {v}")
        else:
            st.write("No keywords added yet.")

    # Do not attempt to clear the widget values after instantiation.
    new_kw_key = st.text_input("New Keyword Key", key="new_kw_key")
    new_kw_value = st.text_input("New Keyword Value", key="new_kw_value")
    if st.button("Add Keyword", key="btn_add_keyword"):
        if new_kw_key.strip():
            st.session_state.keywords.append((new_kw_key.strip(), new_kw_value.strip()))
            st.success("Keyword added!")
            # Removed the lines that reset the keys:
            # st.session_state.new_kw_key = ""
            # st.session_state.new_kw_value = ""
            st.rerun()  # Use st.rerun() to refresh the UI
        else:
            st.error("Please provide a keyword key.")

    # --- Protocols Section ---
    st.markdown("#### Protocols")
    with st.expander("Existing Protocols (click to expand)", expanded=False):
        if st.session_state.protocols:
            for idx, (k, v) in enumerate(st.session_state.protocols, start=1):
                st.write(f"{idx}. {k}: {v}")
        else:
            st.write("No protocols added yet.")
    new_p_key = st.text_input("New Protocol Key", key="new_p_key")
    new_p_value = st.text_input("New Protocol Value", key="new_p_value")
    if st.button("Add Protocol", key="btn_add_protocol"):
        if new_p_key.strip():
            st.session_state.protocols.append((new_p_key.strip(), new_p_value.strip()))
            st.success("Protocol added!")
            # Removed: st.session_state.new_p_key = ""
            # Removed: st.session_state.new_p_value = ""
            st.rerun()
        else:
            st.error("Please provide a protocol key.")

    # Build dictionaries from the lists, with type conversion.
    kw_dict = {k: auto_convert(v) for (k, v) in st.session_state.keywords} if st.session_state.keywords else {}
    p_dict = {k: auto_convert(v) for (k, v) in st.session_state.protocols} if st.session_state.protocols else {}

    return kw_dict, p_dict
