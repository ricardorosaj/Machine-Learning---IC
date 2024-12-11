import pandas as pd
import csv
import util
import os 

"""
Title: Description
Filter genes from necroptosis pathway or protein coding
"""

"""
Function: filter_genes

 Description:
    This function selects only the genes that will be further analysed on the next steps, being protein coding genes
    and genes that are present in te necroptosis pathway.

 Parameters:
 	df - path for dataframe with expression of genes
    filter_by - filter by necroptosis or by protein coding (pc) genes
 	
 Returns:
 	dictionary of expression data for each gene with necroptosis pathway genes removed.
"""

def filter_genes(df, filter_by):
    # "Protein coding genes" were previously selected as interesting genes to be analized
    if filter_by == 'pc':
        genes_to_filter = pd.read_csv(util.path__+'../../../data/external/genes_protein_conding.tsv'.replace('/',os.sep),sep='\t')
        genes_to_filter = genes_to_filter['Ensembl gene ID'].values
        dict_df = {}
        with open(df,'r') as gct:
            read=csv.reader(gct, delimiter=',')
            for line in read:
                if line[0]=='Gene_ID':
                    tissue=line
                    for i in range(len(line)):
                        dict_df[line[i]]=[]
                else:
                    for i in range(len(line)):
                        if line[0].split('.')[0] in genes_to_filter:
                            dict_df[tissue[i]].append(line[i])
        return dict_df                    
    # "Necroptosis genes" are genes that are part of the death of the cell pathway and could affect expression of other genes 
    elif filter_by == 'necrop':
        genes_to_filter = pd.read_csv(util.path__+'../../../data/external/genes_necroptose.csv'.replace('/',os.sep))
        genes_to_filter = genes_to_filter['x'].values
        dict_df = {}
        with open(df,'r') as gct:
            read=csv.reader(gct, delimiter=',')
            for line in read:
                if line[0]=='Gene_ID':
                    tissue=line
                    for i in range(len(line)):
                        dict_df[line[i]]=[]
                else:
                    for i in range(len(line)):
                        if line[0].split('.')[0] in genes_to_filter:
                            dict_df[tissue[i]].append(line[i])
        return dict_df                     
    else:
        print('Error: Specified filtering method not valid')                             
                       

# dict_df_aorta = filter_genes(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_aorta_gtex.csv'.replace('/',os.sep))
# dict_df_coronary = filter_genes(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_coronary_gtex.csv'.replace('/',os.sep))
# pd.DataFrame.from_dict(dict_df_aorta).to_csv(util.path__+'../../../data/interim/gtex_data/resultados/filtered_necrop_log2_tpm_aorta.csv'.replace('/',os.sep), index=False)
# pd.DataFrame.from_dict(dict_df_coronary).to_csv(util.path__+'../../../data/interim/gtex_data/resultados/filtered_necrop_log2_tpm_coronary.csv'.replace('/',os.sep), index=False)
# df_atrial = filter_genes(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_atrial_appendage_gtex.csv'.replace('/',os.sep))
# df_ventricle = filter_genes(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_left_ventricle_gtex.csv'.replace('/',os.sep))
# pd.DataFrame.from_dict(df_atrial).to_csv(util.path__+'../../../data/interim/gtex_data/filtered_protein_coding_log2_tpm_data_atrial_appendage.csv'.replace('/',os.sep), index=False)
# pd.DataFrame.from_dict(df_ventricle).to_csv(util.path__+'../../../data/interim/gtex_data/filtered_protein_coding_log2_tpm_data_left_ventricle.csv'.replace('/',os.sep), index=False)
df_liver = filter_genes(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_Liver_gtex.csv'.replace('/',os.sep), filter_by='pc')
pd.DataFrame.from_dict(df_liver).to_csv(util.path__+'../../../data/interim/gtex_data/filtered_protein_coding_log2_tpm_data_Liver_gtex.csv'.replace('/',os.sep), sep=',', index=False)