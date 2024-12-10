import numpy as np
import pandas as pd
import util 
import csv
import os 
import math

"""
Title: Descrition
Code for filtering samples from specified tissue and normalizing tpm data to log2(tpm+1)
"""

"""
Function: create_log2tpm_for_tissue_from_gtex

Description: Reads which samples are from specified tissue and filter them, creating a dataframe with normalized tpm expression.

Parameters: 
	tissue - Specific tissue from GTEx data (needs to be exactly written like on GTEx)

Returns:
	Saves dataframe in interim/gtex_data folder with specified tissue name
"""


def create_log2tpm_for_tissue_from_gtex(tissue):
    """Creates a dataframe of the expression, in log2tpm, of genes from samples from specified tissue"""
    df = pd.read_csv(util.path__+'../../../data/raw/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'.replace('/',os.sep), usecols=['SAMPID', 'SMTSD'], sep='\t') #File with id and tissue info
    list_tissue = []
    for i in range(len(df)):
        if df['SMTSD'][i]==tissue:
            list_tissue.append(df['SAMPID'][i])
    dict_for_list_tissue={'Gene_ID':[]}
    with open(util.path__+'../../../data/raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'.replace('/',os.sep),'r') as gct: #Open big file and separates the expression data
        read=csv.reader(gct,delimiter='\t')
        next(read)
        next(read)
        for line in read:
            if line[0]=='Name': #Creates keys on dictionaries with the sample name
                tissues=line
                for i in range(len(line)):
                    if line[i] in list_tissue:
                        dict_for_list_tissue[line[i]]=[]
            else:
                dict_for_list_tissue['Gene_ID'].append(line[0])
                for i in range(len(tissues)):
                    if tissues[i] in list_tissue:
                        dict_for_list_tissue[tissues[i]].append(math.log(float(line[i])+1,2))
    pd.DataFrame.from_dict(dict_for_list_tissue).to_csv(util.path__+f'../../../data/interim/gtex_data/log2_tpm_data_{tissue}_gtex.csv'.replace('/',os.sep), index=False)
                    
create_log2tpm_for_tissue_from_gtex('Liver')