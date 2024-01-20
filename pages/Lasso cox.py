import streamlit as st
import pandas as pd
import altair as alt

import fun.enhancer_gene_corr

st.set_page_config(page_title="Lasso-cox", page_icon="📊")
txt = st.text_area('请输入不超过50个基因的名字')

enhancer_choose=0
enhancer_choose_pre=0
option = st.selectbox('是否对基因相关增强子进行风险回归分析',('Yes', 'No'))
if option == 'Yes':
    enhancer_choose=0
    option2 = st.selectbox('采用何种标准',('严格', '非严格'))
    if option2 == '严格':
        enhancer_choose_pre=0
    else:
        enhancer_choose_pre=1
else:
    enhancer_choose=1

option_cancer = st.selectbox('癌症类型',("ACC","BLCA","BRCA","CESC","CHOL","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"))


if st.button("Submit"):
    cox_img=fun.enhancer_gene_corr.lasso_cox(txt,enhancer_choose,enhancer_choose_pre,option_cancer)
    img=cox_img.split('\n')
    for mm in img:
        if mm !='':
            nn=mm.split('\t')
            for nm in nn:
                if nm !='':
                    st.image(nm)
