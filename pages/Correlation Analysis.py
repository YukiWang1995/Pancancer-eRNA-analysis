import streamlit as st
import pandas as pd
import joblib

import fun.enhancer_gene_corr


st.header("Cancer Enhancer database")
gene=st.text_input("Enter Gene name")
enhancer=st.text_input("Enter enhancer id")

user_cancer = st.multiselect('Select Cancer name:',("ACC","BLCA","BRCA","CESC","CHOL","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"))
cancer_list=user_cancer

if st.button("Submit"):
    corr=fun.enhancer_gene_corr.enhancer_gene_corr(gene,enhancer,cancer_list)   
    corr=corr.split('\n')
    i=0
    gene_list=[]
    enhancer_list=[]
    corr_list=[]
    p_list=[]
    cancer_list=[]
    for line1 in corr:
        i=i+1
        if i >1:
            a=line1.split('\n')[0]
            a=a.split('\t')
            if len(a)>2:
                gene_list.append(a[0])
                enhancer_list.append(a[1])
                corr_list.append(str(a[2]))
                p_list.append(str(a[3]))
                cancer_list.append(str(a[4]))
    
    
    data = {
        'Gene':gene_list,
        'Enhancer':enhancer_list,
        'Pearson_correlation': corr_list,
        'p_value': p_list,
        'cancer': cancer_list
    }
    df = pd.DataFrame(data)
    df.index = df.index + 1
    st.write('结果展示：')
    st.dataframe(df)
    

        
