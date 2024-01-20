import enhancer_gene_corr

cancer_list=['ACC','BRCA']     
gene_search='SDF4'  
enhancer_choose=0
user_distance='100'

enhancer=enhancer_gene_corr.search(gene_search,enhancer_choose,cancer_list,user_distance)
predict_enhacner=enhancer.split('\n')
gene_list=[]
enhancer_list=[]
corr_list=[]
p_list=[]
cancer_list=[]
for mm in predict_enhacner:

    if 'enhancerPearson_correlation' not in mm:
        a=mm
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
print(data)