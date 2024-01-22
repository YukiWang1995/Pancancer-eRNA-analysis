import os
import streamlit as st
import pandas as pd
import joblib

import fun.enhancer_gene_corr

aaa=(os.getcwd())

st.header("Cancer Enhancer database")
st.header(aaa)


chr_number = ['chr1', 'chr2', 'chr3','chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9','chr10', 'chr11', 'chr12','chr13', 'chr14', 'chr15','chr16', 'chr17', 'chr18','chr19', 'chr20', 'chr21','chr22', 'chrX', 'chrY']
selected_option = st.selectbox('Select an option', chr_number)
start=st.text_input("Enter loc start") 
end=st.text_input("Enter loc end") 

if st.button("Submit"):
    ID=fun.enhancer_gene_corr.enhancer_loc_search(selected_option,start,end)  
    enhancer_ID=[]
    enhancer_chrnumber=[]
    enhancer_start=[]
    enhancer_end=[]
    ID=ID.split('\n')
    for mm in ID:
        if 'enhancer_id' not in mm:
            a=mm.split('\t')
            if len(a)>2:
                    enhancer_ID.append(a[0])
                    enhancer_chrnumber.append(a[1])
                    enhancer_start.append(str(a[2]))
                    enhancer_end.append(str(a[3]))
    
    data = {
                'enhancer_ID':enhancer_ID,
                'chrnumber':enhancer_chrnumber,
                'start': enhancer_start,
                'end': enhancer_end
    }
    df = pd.DataFrame(data)
    df.index = df.index + 1
    st.write('结果展示：')
    st.dataframe(df)
