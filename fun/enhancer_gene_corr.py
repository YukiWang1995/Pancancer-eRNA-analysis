import numpy as np
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from lifelines.datasets import load_regression_dataset
from lifelines import NelsonAalenFitter, CoxPHFitter, KaplanMeierFitter
from lifelines.datasets import load_waltons
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from lifelines.statistics import logrank_test
from matplotlib.backends.backend_pdf import PdfPages 
from IPython.core.pylabtools import figsize
from collections import Counter
from lifelines import CoxPHFitter


def enhancer_loc_search(chr_number,statr,end):
    f1=open('./rawdata/enhancer.id.out','r')
    write='enhancer_id'+'\t'+'chrumber'+'\t'+'start'+'\t'+'end'+'\n'
    enhancer_exist=0
    for line1 in f1:
        a=line1.split('\n')[0]
        a=a.split('\t')
        enhancer_id=a[-1]
        enhancer_chr=a[0]
        if enhancer_chr==str(chr_number):
            enh_st_b=int(a[1])+1
            enh_ed_b=int(a[2])-1
            if enh_st_b > int(statr) and enh_ed_b < int(end):
                write=write+enhancer_id+'\t'+str(chr_number)+'\t'+str(a[1])+'\t'+str(a[2])+'\n'
                enhancer_exist=1
    if enhancer_exist==1:
        return(write)
    else:
        write='Na'
        return(write)
                

def enh_id_loc(name):
    name=name.split('\n')[0]
    if name[0:2]=='eh':
        type_id=0
        f1=open('./rawdata/enhancer.id.out','r')
        for line1 in f1:
            a=line1.split('\n')[0]
            a=a.split('\t')
            enhancer_id=a[-1]
            if enhancer_id == name:
                type_id=1
                write=line1.split('\n')[0]
                return(write)
        if type_id==0:
            write='Wrong ID'
            return(write)
            
    elif name[0:2]=='ch':
        type_id=0
        f1=open('./rawdata/enhancer.id.out','r')
        for line1 in f1:
            a=line1.split('\n')[0]
            a=a.split('\t')
            enhancer_loc=str(a[0])+':'+str(a[1])+'-'+str(a[2])
            if enhancer_loc == name:
                type_id=1
                write=line1.split('\n')[0]
                return(write)
        if type_id==0:
            write='Wrong ID'
            return(write)
    else:
        write='Wrong ID'
        return(write)

def enhancer_gene_distance(gene,enhancer):
    f1=open('./rawdata/gene.hg19.position','r')
    for line1 in f1:
        a=line1.split('\n')[0]
        a=a.split('\t')
        genename=a[-1]
        if genename == gene:
            geneloc_chr=a[0]
            geneloc_st=int(a[1])
            geneloc_ed=int(a[2])
    f1.close()
    
    enhancer_true=enh_id_loc(enhancer)
    enhancer_true=enhancer_true.split('\t')
    enhancerloc_chr=enhancer_true[0]
    enhancerloc_st=int(enhancer_true[1])
    enhancerloc_ed=int(enhancer_true[2])
    
    if geneloc_chr==enhancerloc_chr:
        dis1=abs(enhancerloc_st-geneloc_st)
        dis2=abs(enhancerloc_st-geneloc_ed)
        dis3=abs(enhancerloc_ed-geneloc_st)
        dis4=abs(enhancerloc_ed-geneloc_ed)
        dis=min(dis1,dis2,dis3,dis4)
        return(dis)
    else:
        dis=-100
        #write='不在同一染色体'
        return(dis)
    
def enhancer_expression_boxplot(enhancer,cancer_type,tumor_normal):  
    file_png='enhancer_box.png'
    
    enhancer_true=enh_id_loc(enhancer)
    enhancer_true=enhancer_true.split('\t')
    enhancer_true=enhancer_true[0]+':'+str(enhancer_true[1])+'-'+str(enhancer_true[2])
    
    dict_box={}
    for mm in cancer_type:
        cancer=mm
        file1='./rawdata/enhancers/TCGA_'+str(cancer)+'_FANTOM5_60k_eRNA_v2.tsv'
        f1=open(file1,'r')
        list_tumor=[]
        list_normal=[]
        dict1_tumor={}
        dict2_normal={}
        
        i=0
        for line1 in f1:
            i=i+1
            if i ==1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                t=0
                for name in a[0:]:
                    t=t+1
                    b=name.split('_')
                    if b[-1]=='tumor':
                        b=b[0].split('-')
                        name=b[0]+'-'+b[1]+'-'+b[2]     
                        dict1_tumor[t]=name
                    else:
                        b=b[0].split('-')
                        name=b[0]+'-'+b[1]+'-'+b[2]     
                        dict2_normal[t]=name

            if i >1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                enhancer_name=a[0].split('_')[-1]
                if enhancer_name == enhancer_true:
                    enhancer_exist=1
                    t=0
                for enhancer_expression in a[1:]:
                    t=t+1
                    if t in dict1_tumor.keys():
                        if float(enhancer_expression) >0:
                            list_tumor.append(float(enhancer_expression))
                    elif t in dict2_normal.keys():
                        if float(enhancer_expression) >0:
                            list_normal.append(float(enhancer_expression))
       
        f1.close()
        
        if len(list_tumor) <2:
            list_tumor.append(0)
        if len(list_normal)<2:
            list_normal.append(0)
        
        cancer_tumor=str(cancer)+'_tumor'
        cancer_normal=str(cancer)+'_normal'
        dict_box[cancer_tumor]=list_tumor
        dict_box[cancer_normal]=list_normal
    
    list_plot=[]
    if tumor_normal==0:
        for mm in dict_box:
            if '_tumor' in mm:
                list_plot.append(dict_box[mm])
    else:
        for mm in dict_box:
            list_plot.append(dict_box[mm])
    
    plt.boxplot(list_plot,patch_artist=True,widths=0.4,meanprops={'marker':'+','markerfacecolor':'k','markeredgecolor':'k','markersize':5})        
    plt.savefig(file_png)
    write=file_png
    return(write)
                
    
def enhancer_gene_corr(gene,enhancer,cancer_type):  
    write=str('gene')+'\t'+str('enhancer')+str('Pearson_correlation')+'\t'+str('P_value')+'\t'+str('cancer')+'\n'
    for mm in cancer_type:
        cancer=mm
        file1='./rawdata/Gene/'+str(cancer)
        file2='./rawdata/enhancers/TCGA_'+str(cancer)+'_FANTOM5_60k_eRNA_v2.tsv'
        f1=open(file1,'r')
        i=0
        dict1_name={}
        dict2_exp={}
        
        gene_exist=0
        enhancer_exist=0
        
        for line1 in f1:
            i=i+1
            if i==1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                t=0
                for name in a[1:]:
                    t=t+1
                    b=name.split('-')
                    if b[-1]=='01':
                        name=b[0]+'-'+b[1]+'-'+b[2]                        
                        dict1_name[t]=name
            if i >1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                genename=a[0]
                if genename == gene:
                    gene_exist=1
                    t=0
                    for gene_expression in a[1:]:
                        t=t+1
                        if t in dict1_name.keys():
                            sample=dict1_name[t]
                            dict2_exp[sample]=str(gene_expression)
        f1.close()
        
        enhancer_true=enh_id_loc(enhancer)
        enhancer_true=enhancer_true.split('\t')
        enhancer_true=enhancer_true[0]+':'+str(enhancer_true[1])+'-'+str(enhancer_true[2])

        f2=open(file2,'r')
        i=0
        dict1_enahcner_name={}
        for line1 in f2:
            i=i+1
            if i ==1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                t=0
                for name in a[0:]:
                    t=t+1
                    b=name.split('_')
                    if b[-1]=='tumor':
                        b=b[0].split('-')
                        name=b[0]+'-'+b[1]+'-'+b[2]       
                        if name in dict2_exp.keys():
                            dict1_enahcner_name[t]=name

            if i >1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                enhancer_name=a[0].split('_')[-1]
                if enhancer_name == enhancer_true:
                    enhancer_exist=1
                    t=0
                    for enhancer_expression in a[1:]:
                        t=t+1
                        if t in dict1_enahcner_name.keys():
                            sample=dict1_enahcner_name[t]
                            dict2_exp[sample]=dict2_exp[sample]+'\t'+str(enhancer_expression)        
        f2.close()
        
        list_gene=[]
        list_enhancer=[]
        for mm in dict2_exp.keys():
            list_gene.append(float(dict2_exp[mm].split('\t')[0]))
            list_enhancer.append(float(dict2_exp[mm].split('\t')[-1]))
        

        pccs = pearsonr(list_gene, list_enhancer)      
        if gene_exist==1 and enhancer_exist==1:
            write=write+str(gene)+'\t'+str(enhancer)+'\t'+str(pccs[0])+'\t'+str(pccs[1])+'\t'+str(cancer)+'\n'
        else:
            write=write+str(gene)+'\t'+str(enhancer)+'\t'+'Na'+'\t'+'Na'+'\t'+str(cancer)+'\n'
    
    return(write)
   
     
def search(gene_search,enhancer_choose,cancer_list,user_distance):
    user_distance=int(user_distance)*1000
    if enhancer_choose==0:
        write=str('gene')+'\t'+str('enhancer')+str('Pearson_correlation')+'\t'+str('P_value')+'\t'+str('cancer')+'\n'
        f1=open('../../rawdata/enhancer-gene.txt','r')
        for line1 in f1:
            a=line1.split('\n')[0]
            a=a.split('\t')
            gene=a[8]
            if gene_search==gene:
                enhancer=str(a[0])+':'+str(a[1])+'-'+str(a[2])
                enhancer_true=enh_id_loc(enhancer)
                enhancer_id=enhancer_true=enhancer_true.split('\t')[-1]
                distance=enhancer_gene_distance(gene_search,enhancer_id)
                if distance >0 and  distance< user_distance:
                    predict_enhancer=enhancer_gene_corr(gene_search,enhancer_id,cancer_list)
                    write=write+predict_enhancer+'\n'
    #此处需要重新构建表格
    if enhancer_choose==1:
        write=str('gene')+'\t'+str('enhancer')+str('Pearson_correlation')+'\t'+str('P_value')+'\t'+str('cancer')+'\n'
        f1=open('./rawdata/enhancer-gene.txt','r')
        for line1 in f1:
            a=line1.split('\n')[0]
            a=a.split('\t')
            gene=a[8]
            if gene_search==gene:
                enhancer=str(a[0])+':'+str(a[1])+'-'+str(a[2])
                enhancer_true=enh_id_loc(enhancer)
                enhancer_id=enhancer_true=enhancer_true.split('\t')[-1]
                distance=enhancer_gene_distance(gene_search,enhancer_id)
                if distance >0 and  distance< user_distance:
                    predict_enhancer=enhancer_gene_corr(gene_search,enhancer_id,cancer_list)
                    write=write+predict_enhancer+'\n'
                    
    return(write)
    
def survival(enhancer_choose,cancer_list,enhancer_choose_time,survival_group):
    
    enhancer_true=enh_id_loc(enhancer_choose)
    enhancer_true=enhancer_true.split('\t')
    enhancer_true=enhancer_true[0]+':'+str(enhancer_true[1])+'-'+str(enhancer_true[2])
    enhancer_exist=0
    write=''
    for mm in cancer_list:
        plt.clf()
        cancer=mm
        file_png=str(cancer)+'exam.png'
        file1='./rawdata/enhancers/TCGA_'+str(cancer)+'_FANTOM5_60k_eRNA_v2.tsv'
        f1=open(file1,'r')
        i=0
        dict1_enahcner={}
        dict2_enhancer={}
        list_enhancer_expression=[]
        for line1 in f1:
            i=i+1
            if i ==1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                t=0
                for name in a[0:]:
                    t=t+1
                    b=name.split('_')
                    if b[-1]=='tumor':
                        b=b[0].split('-')
                        name=b[0]+'-'+b[1]+'-'+b[2]  
                        dict1_enahcner[t]=name

            
            if i >1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                enhancer_name=a[0].split('_')[-1]
                if enhancer_name == enhancer_true:
                    enhancer_exist=1
                    t=0
                    for enhancer_expression in a[1:]:
                        t=t+1
                        if t in dict1_enahcner.keys():
                            sample=dict1_enahcner[t]
                            dict2_enhancer[sample]=float(enhancer_expression)
                            list_enhancer_expression.append(float(enhancer_expression))
        f1.close()
            

        list_enhancer_expression.sort()
        m1=int(len(list_enhancer_expression)*float(survival_group/100))
        middle=list_enhancer_expression[m1]
            
        x=[]
        y=[]
        E=[]
            
        file2='./rawdata/Survival.txt'
        f2=open(file2,'r')
        i=0
        for line1 in f2:
            i=i+1
            if i >1:
                a=line1.split('\t')
                cancer_name=a[2]
                if cancer_name==mm:
                    sample_name=a[1]
                    if sample_name in dict2_enhancer.keys():
                        if enhancer_choose_time ==0:
                            if a[-9] != '' and a[-8] != '':
                                if dict2_enhancer[sample_name] > middle:
                                    x.append('high')
                                    y.append(int(a[-8]))
                                    E.append(int(a[-9]))
                                else:
                                    x.append('low')
                                    y.append(int(a[-8]))
                                    E.append(int(a[-9]))                        
                        else:
                            if a[-7] != '' and a[-6] != '':
                                if dict2_enhancer[sample_name] > middle:
                                    x.append('high')
                                    y.append(int(a[-6]))
                                    E.append(int(a[-7]))
                                else:
                                    x.append('low')
                                    y.append(int(a[-6]))
                                    E.append(int(a[-7]))  
        f2.close()

        data={'T':y,'E':E,'group':x}
        labels=range(1,len(y)+1)       
        df = pd.DataFrame(data, index=labels)
        kmf = KaplanMeierFitter()
        kmf.fit(df['T'], event_observed=df['E'])
        median_ = kmf.median_survival_time_
        median_confidence_interval_ = median_survival_times(kmf.confidence_interval_)
        groups = df['group']
        ix2 = (groups == 'high')
        kmf.fit(df['T'][ix2], df['E'][ix2], label='high')#共享一个画布
        ax = kmf.plot(color='red',ci_show=False)
        control_median_confidence_interval_ = median_survival_times(kmf.confidence_interval_)
        ix = (groups == 'low')
        kmf.fit(df['T'][ix], df['E'][ix], label='low')
        ax = kmf.plot(ax=ax,color='blue',ci_show=False)
        treatment_median_confidence_interval_ = median_survival_times(kmf.confidence_interval_)
        lr = logrank_test(df['T'][ix], df['T'][ix2],df['E'][ix], df['E'][ix2], alpha=.95)
        plt.savefig(file_png)

        
        write=write+'\n'+str(file_png)
    return(write)



def cox_dataframe(name_list,type_list,cancer_name):
    if type_list==0:
        file1='./rawdata/Gene/'+str(cancer_name)
        f1=open(file1,'r')
        gene_list=[]
        dict1_name={}
        dict2_exp={}
        i=0
        for line1 in f1:
            i=i+1
            if i==1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                t=0
                for name in a[1:]:
                    t=t+1
                    b=name.split('-')
                    if b[-1]=='01':
                        name=b[0]+'-'+b[1]+'-'+b[2]                        
                        dict1_name[t]=name
                
            if i >1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                genename=a[0]
                if genename in name_list:
                    t=0
                    gene_list.append(genename)
                    for gene_expression in a[1:]:
                        t=t+1
                        if t in dict1_name.keys():
                            sample=dict1_name[t]
                            if sample not in dict2_exp.keys():
                                dict2_exp[sample]=str(gene_expression)
                            else:
                                dict2_exp[sample]=dict2_exp[sample]+'\t'+str(gene_expression)
        f1.close()  
    else:
        enhancer_name_list=[]
        for enhancer_mm in name_list:
            enhancer_true=enh_id_loc(enhancer_mm)
            enhancer_true=enhancer_true.split('\t')
            enhancer_true=enhancer_true[0]+':'+str(enhancer_true[1])+'-'+str(enhancer_true[2])
            enhancer_name_list.append(enhancer_true)
        
        file1='./rawdata/enhancers/TCGA_'+str(cancer_name)+'_FANTOM5_60k_eRNA_v2.tsv'
        f1=open(file1,'r')
        i=0
        gene_list=[]
        dict1_name={}
        dict2_exp={}
        for line1 in f1:
            i=i+1
            if i ==1:
                a=line1.split('\n')[0]
                a=a.split('\t')
                t=0
                for name in a[0:]:
                    t=t+1
                    b=name.split('_')
                    if b[-1]=='tumor':
                        b=b[0].split('-')
                        name=b[0]+'-'+b[1]+'-'+b[2]  
                        dict1_name[t]=name
                    
            if i >1:
                 a=line1.split('\n')[0]
                 a=a.split('\t')
                 enhancer_name=a[0].split('_')[-1]
                 if enhancer_name in enhancer_name_list:
                     gene_list.append(enhancer_name)
                     enhancer_exist=1
                     t=0
                     for enhancer_expression in a[1:]:
                         t=t+1
                         if t in dict1_name.keys():
                             sample=dict1_name[t]
                             if sample not in dict2_exp.keys():
                                dict2_exp[sample]=str(enhancer_expression)
                             else:
                                dict2_exp[sample]=dict2_exp[sample]+'\t'+str(enhancer_expression)
        f1.close()  

    
    dict3_cox={}
    file2='./rawdata/Survival.txt'
    f2=open(file2,'r')
    i=0
    for line1 in f2:
        i=i+1
        if i >1:
            a=line1.split('\t')
            cancer_name_sur=a[2]
            if cancer_name==cancer_name_sur:
                sample_name=a[1]
                if sample_name in dict2_exp.keys():
                    if a[-9] != '' and a[-8] != '':
                        dict3_cox[sample_name]=str(a[-9])+'\t'+str(a[-8])+'\t'+dict2_exp[sample_name]
    f2.close()
    
    cox_labels={0:'Event',1:'Time'}
    t=1
    for cox_mm in gene_list:
        t=t+1
        cox_labels[t]=cox_mm
        
    data=[]
    for cox_mm_2 in dict3_cox.keys():
        list_cox_mm=[]
        list_cox_mm=dict3_cox[cox_mm_2].split('\t')
        data.append(list_cox_mm)
        list_cox_mm=[]
        
    df = pd.DataFrame(data)
    df.rename(columns=cox_labels, inplace=True)
    return(df)
    

def lasso_cox(txt,enhancer_choose,enhancer_choose_pre,option_cancer):
    gene_list=[]
    enhancer_list=[]
    cancer_list=[]
    cancer_list.append(option_cancer)
    
    file_gene_png=str(option_cancer)+'cph.png'
    file_enhancer_png=str(option_cancer)+'_enhancer_cph.png'
    write=''
    
    gene_txt=txt.split('\n')
    for mm in gene_txt:
        gene_mm=mm.split(',')
        for nn in gene_mm:
            if nn !='':
                gene_list.append(nn)
                
    user_distance=2000
    if enhancer_choose==0:
        for gene in gene_list:
            enhancer=search(gene,enhancer_choose_pre,cancer_list,user_distance)
            enh_1=(enhancer.split('\n'))
            for enh_mm in enh_1:
                enh_2=enh_mm.split('\t')
                if len(enh_2) >2:
                    if 'eh' in enh_2[1]:
                        if enh_2[1] not in enhancer_list:
                            enhancer_list.append(enh_2[1])
    else:
        enhancer_list=[]
                        
    
    data_cox=cox_dataframe(gene_list,0,str(option_cancer))
    cph=CoxPHFitter()
    cph.fit(data_cox, duration_col='Time', event_col='Event')
    cph.plot()
    plt.savefig(file_gene_png)
    write=write+'\t'+file_gene_png
    
    if len(enhancer_list) >1 :   
        data_cox=cox_dataframe(enhancer_list,1,str(option_cancer))
        cph=CoxPHFitter()
        cph.fit(data_cox, duration_col='Time', event_col='Event')
        cph.plot()
        plt.savefig(file_enhancer_png)
        write=write+'\t'+file_enhancer_png
        
    return(write)
