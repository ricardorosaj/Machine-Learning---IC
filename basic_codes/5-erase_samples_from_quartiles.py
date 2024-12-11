import pandas as pd
import util

"""
Title: Description
Erases samples from quartiles of time of death that will not be used as input to the algorithm. Fixes column names for running 
"""

"""
Function: erase_quartiles_and_fix_cols_names
 
 Description:
    As a final preparation step to run the Decision Tree algorithm, this function erases the samples that are 
    within 24h of the time of death for the extraction of the sample from the patient.
 
 Parameters:
 	df - dataframe with expression of genes
 	quartiles - list of quartiles that you wish to exclude. Example: [0,1] (starts from 0)
 
 Returns:
 	Dataframe transformed with time of death of all samples and the quartile of time of death of each sample
"""

def erase_quartiles_and_fix_cols_names(df,quartiles):
    df.drop('SMTSISCH', axis=1, inplace=True)
    df.rename(columns={'sample_id':'Unnamed: 0'}, inplace=True)
    df.set_index('', inplace=True)
    samples_to_drop = []
    for i in range(len(df)):
        quartil = df['quartiles'].values[i]
        if quartil not in quartiles:
            samples_to_drop.append(df.index.values[i])
    df.drop(samples_to_drop, inplace=True)
    df.drop('quartiles', axis=1, inplace=True)

    new_cols = []
    for i in range(len(df.columns)):
        new_cols.append(df.columns.values[i].split('.')[0])
    df.columns = new_cols    
    return df

df = pd.read_csv(util.path__+'../../../data/interim/gtex_data/quartiles_log2_tpm_coronary.csv')
erase_quartiles_and_fix_cols_names(df,quartiles=[2,3])