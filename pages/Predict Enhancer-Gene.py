import streamlit as st
import pandas as pd
import joblib

import fun.enhancer_gene_corr

st.header("Cancer Enhancer database")
gene=st.text_input("Enter Gene name")

user_cancer = st.multiselect('Select Cancer name:',("ACC","BLCA","BRCA","CESC","CHOL","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"))
cancer_list=user_cancer

user_distance = st.slider("增强子与基因距离(Kb)", min_value=0, max_value=2000,value=2000)
#st.write("你的评分是：", user_distance) 

enhancer_choose=0
options = ['Strict', 'Tolerant']
selected_option = st.radio('Select an option: ', options)
if selected_option == 'Strict':
    st.write("You have chosen Strict.")
    enhancer_choose=0
elif selected_option == 'Tolerant':
    st.write("You have chosen Tolerant.")
    enhancer_choose=0


if st.button("Submit"):
    gene=str(gene)
    enhancer=fun.enhancer_gene_corr.search(gene,enhancer_choose,cancer_list,user_distance)
    predict_enhacner=enhancer.split('\n')
    gene_list=[]
    enhancer_list=[]
    corr_list=[]
    p_list=[]
    cancer_list=[]
    for mm in predict_enhacner:
        if 'enhancerPearson_correlation' not in mm:
            a=mm.split('\t')
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
        


    
    
    
