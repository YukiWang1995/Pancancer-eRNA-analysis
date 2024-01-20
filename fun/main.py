import streamlit as st
import pandas as pd
import joblib

import enhancer_gene_search

st.header("Cancer Enhancer database")
gene=st.text_input("Enter Gene name")
if st.button("Submit"):
    gene=str(gene)
    enhancer=enhancer_gene_search.search(gene)
    st.text(f"This instance is a {enhancer}")

