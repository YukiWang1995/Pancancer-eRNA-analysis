import streamlit as st
import time
import numpy as np

import fun.enhancer_gene_corr

st.header("Enhancer expression")

enhancer=st.text_input("Enter enhancer id")
user_cancer = st.multiselect('Select Cancer name:',("ACC","BLCA","BRCA","CESC","CHOL","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"))
cancer_list=user_cancer
enhancer_tumor_normal=0
options = ['Only Tumor', 'Tumor and Normal']
selected_option = st.radio('Only Tumor: ', options)
if selected_option == 'Only Tumor':
    st.write("You have chosen Only Tumor")
    enhancer_tumor_normal=0
elif selected_option == 'Tumor and Normal':
    st.write("You have chosen Tumor and Normal")
    enhancer_tumor_normal=1


if st.button("Submit"):
    enhancer_expression=fun.enhancer_gene_corr.enhancer_expression_boxplot(enhancer,cancer_list,enhancer_tumor_normal)
    img=enhancer_expression
    st.image(img)
    