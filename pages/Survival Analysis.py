import streamlit as st
import pandas as pd

import fun.enhancer_gene_corr




st.header("Cancer Enhancer database")
enhancer=st.text_input("Enter Enhancer ID")

user_cancer = st.multiselect('Select Cancer name:',("ACC","BLCA","BRCA","CESC","CHOL","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"))
cancer_list=user_cancer

enhancer_choose_time=0
options = ['Overall Survival', 'Disease Free Survival (RFS)']
selected_option = st.radio('Select an option: ', options)
if selected_option == 'Overall Survival':
    st.write("You have chosen Overall Survival")
    enhancer_choose_time=0
elif selected_option == 'Disease Free Survival (RFS)':
    st.write("You have chosen Disease Free Survival (RFS)")
    enhancer_choose_time=1

survival_group = st.slider("Group Cutoff", min_value=0, max_value=100,value=50)

image_list=[]
if st.button("Submit"):
    img_path=fun.enhancer_gene_corr.survival(enhancer,cancer_list,enhancer_choose_time,survival_group)
    img=img_path.split('\n')
    for mm in img:
        if mm !='':
            st.image(mm)
