import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import util 
import os

"""
Title: Description
Plots the expression of genes from necroptosis pathway, by removing stable genes in quartiles and plotting their exp or by the mean of all genes by sample
"""

def transform_data(df):
    # Save data in readable format
    df.columns=df.iloc[0]
    df = df[1:]
    df = df.reset_index()
    df.to_csv('temporario.csv', index=False)
    df = pd.read_csv('temporario.csv')
    df = df.rename(columns={'index':'sample_id'})
    df = df.set_index('sample_id')
    return df

"""
Function: create_df_for_plot

 Description: 
    This function creates the dataframe for plotting time of death of each patient separated by quartile 

 Parameters:
 	df - dataframe with expression of genes
 	tissue - from which tissue the samples of df are

 Returns:
 	Dataframe transformed with time of death of all samples and the quartile of time of death of each sample
"""

def create_df_for_plot(df,tissue):
    time_of_death = pd.read_csv(util.path__+'../../../data/raw/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t', usecols=['SAMPID','SMTSISCH'])
    samples = df.columns.values
    dict_df = {'SAMPID':[],'SMTSISCH':[]}
    for i in range(len(time_of_death)):
        sample = time_of_death['SAMPID'][i]
        if sample in samples:
            dict_df['SAMPID'].append(sample)
            dict_df['SMTSISCH'].append(time_of_death['SMTSISCH'][i])
    for i in range(len(dict_df['SAMPID'])):
        dict_df['SMTSISCH'][i] = dict_df['SMTSISCH'][i]
    dict_df['quartiles'] = pd.qcut(dict_df['SMTSISCH'], q=4, labels=False)

    time_of_death_df = pd.DataFrame(dict_df, index=dict_df['SAMPID'])
    time_of_death_df.drop('SAMPID', inplace=True, axis=1)
    df = df.T
    df = transform_data(df)
    quartile = df.join(time_of_death_df, how='inner')
    quartile.to_csv(util.path__+f'../../../data/interim/gtex_data/quartiles_log2_tpm_{tissue}.csv') #if want to save file, uncomment this line
    return quartile

"""
Function: remove_stable_genes

 Description: 
    As we want to analyse genes that affect the expression of other by being present in the necroptosis pathway, this functions selects
    only the genes that change its expression by varying the time of death of the patients.

 Parameters:
 	df - dataframe with expression of genes
 	
 Returns:
 	Dataframe removing the genes from necrop pathway that present stable expression between all quartiles
"""

def remove_stable_genes(df):
    genes_to_drop = []
    df = df.drop(['quartiles','SMTSISCH'], axis=1)
    q0 = df[df['quartiles']==0].mean(axis=0).values[:-2]
    q1 = df[df['quartiles']==1].mean(axis=0).values[:-2]
    q2 = df[df['quartiles']==2].mean(axis=0).values[:-2]
    q3 = df[df['quartiles']==3].mean(axis=0).values[:-2]
    std = np.std(df)[:-2]
    for i in range(len(std)):
        if abs(q0[i]-q3[i]) <= std[i]: # Checks if the difference of the expression values of the genes in the first and 4th quartile is 
            genes_to_drop.append(std.index.values[i]) # lower than the standard deviation of the expression value of this gene
    new_df = df.drop(genes_to_drop,axis=1)
    return new_df

"""
Function: df_mean_exp_by_time_of_death

 Description: 
    This function calculates the mean expression of genes from a single patient for each quartile.
 
 Parameters:
 	df - dataframe with expression of genes
 	
 Returns:
 	Dataframe with the mean expression of all genes for a single sample, time of death and quartile of each sample
"""

def df_mean_exp_by_time_of_death(df):
    quartile = df['quartiles']
    time = df['SMTSISCH']
    q0 = df[df['quartiles']==0].drop(['SMTSISCH','quartiles'],axis=1)
    q0 = q0.mean(axis=1)
    q1 = df[df['quartiles']==1].drop(['SMTSISCH','quartiles'],axis=1)
    q1 = q1.mean(axis=1)
    q2 = df[df['quartiles']==2].drop(['SMTSISCH','quartiles'],axis=1)
    q2 = q2.mean(axis=1)
    q3 = df[df['quartiles']==3].drop(['SMTSISCH','quartiles'],axis=1)
    q3 = q3.mean(axis=1)    
    exp = pd.concat([q0,q1,q2,q3])
    new_df = pd.DataFrame({'exp':exp,'SMTSISCH':time,'quartiles':quartile}, index=exp.index)
    return new_df

"""
Function: save_plot

 Description:
    This function plots, for each tissue that is being analysed, the scatter plot of the expression of genes
    separated by quartile of the time of deat of each patient.

 Parameters:
 	df - dataframe with expression of genes
    tissue - values 'ventricle', 'atrial' or 'coronary'
 	
 Returns:
 	Saves the plots of expression by time of death of each significant gene from specified tissue. Prints error message if tissue is not valid.
"""

def save_plot(df,tissue):
    if tissue == 'ventricle':
        for gene in df.columns[:-2].values:
            plot = df[[gene,'SMTSISCH','quartiles']]
            sns.scatterplot(x=plot['SMTSISCH'],y=plot[gene],hue=plot['quartiles'],palette='colorblind')
            plt.title(f'Expression of {gene} gene by time of collection')
            plt.xlabel('Time of collection post-mortem (min)')
            plt.ylabel(r'Expression in $log_2{(tpm+1)}$')
            plt.savefig(util.path__+f'../../../analysis/results/gtex_results/genes_from_necrop/left_ventricle/lef_ventricle_{gene}.png',format='png')
            plt.close()
    elif tissue == 'atrial':
        for gene in df.columns[:-2].values:
            plot = df[[gene,'SMTSISCH','quartiles']]
            sns.scatterplot(x=plot['SMTSISCH'],y=plot[gene],hue=plot['quartiles'],palette='colorblind')
            plt.title(f'Expression of {gene} gene by time of collection')
            plt.xlabel('Time of collection post-mortem (min)')
            plt.ylabel(r'Expression in $log_2{(tpm+1)}$')
            plt.savefig(util.path__+f'../../../analysis/results/gtex_results/genes_from_necrop/atrial_appendage/atrial_appendage_{gene}.png',format='png')
            plt.close()
    elif tissue == 'coronary':
        for gene in df.columns[:-2].values:
            plot = df[[gene,'SMTSISCH','quartiles']]
            sns.scatterplot(x=plot['SMTSISCH'],y=plot[gene],hue=plot['quartiles'],palette='colorblind')
            plt.title(f'Expression of {gene} gene by time of collection')
            plt.xlabel('Time of collection post-mortem (min)')
            plt.ylabel(r'Expression in $log_2{(tpm+1)}$')
            plt.savefig(util.path__+f'../../../analysis/results/gtex_results/genes_from_necrop/coronary/coronary_{gene}.png',format='png')
            plt.close()        
    else:
        print('Error: wrong tissue inserted')                



# df_ventricle = pd.read_csv(util.path__+'../../../data/interim/gtex_data/filtered_necrop_log2_tpm_data_left_ventricle.csv'.replace('/',os.sep))
# df_atrial = pd.read_csv(util.path__+'../../../data/interim/gtex_data/filtered_necrop_log2_tpm_data_atrial_appendage.csv'.replace('/',os.sep))

# quartile_ventricle = create_df_for_plot(df_ventricle,time_of_death)
# quartile_ventricle_filtered = remove_stable_genes(quartile_ventricle)
# save_plot(quartile_ventricle_filtered,'ventricle')

# quartile_atrial = create_df_for_plot(df_atrial, time_of_death)
# quartile_atrial_filtered = remove_stable_genes(quartile_atrial)
# save_plot(quartile_atrial_filtered,'atrial')

df_coronary = pd.read_csv(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_coronary_gtex.csv')
df_coronary = create_df_for_plot(df_coronary,'coronary')
# df_coronary_filtered = df_mean_exp_by_time_of_death(df_coronary)
# sns.scatterplot(x=df_coronary_filtered['SMTSISCH'], y=df_coronary_filtered['exp'], hue=df_coronary_filtered['quartiles'], palette='colorblind')
# plt.title('Mean expression of necroptosis genes by time of collection - Coronary')
# plt.xlabel('Time of collection post-mortem (min)')
# plt.ylabel(r'Mean expression in $log_2{(tpm+1)}$')
# plt.savefig(util.path__+'../../../analysis/results/gtex_results/genes_from_necrop/coronary/mean_coronary.png',format='png')
# plt.close()   

df_atrial = pd.read_csv(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_atrial_appendage_gtex.csv'.replace('/',os.sep))
df_atrial = create_df_for_plot(df_atrial,'atrial_appendage')
# df_atrial_filtered = df_mean_exp_by_time_of_death(df_atrial)
# sns.scatterplot(x=df_atrial_filtered['SMTSISCH'], y=df_atrial_filtered['exp'], hue=df_atrial_filtered['quartiles'], palette='colorblind')
# plt.title('Mean expression of necroptosis genes by time of collection - Atrial Appendage')
# plt.xlabel('Time of collection post-mortem (min)')
# plt.ylabel(r'Mean expression in $log_2{(tpm+1)}$')
# plt.savefig(util.path__+'../../../analysis/results/gtex_results/genes_from_necrop/atrial_appendage/mean_atrial.png',format='png')
# plt.close()   

df_ventricle = pd.read_csv(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_left_ventricle_gtex.csv'.replace('/',os.sep))
df_ventricle = create_df_for_plot(df_ventricle,'left_ventricle')
# df_ventricle_filtered = df_mean_exp_by_time_of_death(df_ventricle)
# sns.scatterplot(x=df_ventricle_filtered['SMTSISCH'], y=df_ventricle_filtered['exp'], hue=df_ventricle_filtered['quartiles'], palette='colorblind')
# plt.title('Mean expression of necroptosis genes by time of collection - Left Ventricle')
# plt.xlabel('Time of collection post-mortem (min)')
# plt.ylabel(r'Mean expression in $log_2{(tpm+1)}$')
# plt.savefig(util.path__+'../../../analysis/results/gtex_results/genes_from_necrop/left_ventricle/mean_ventricle.png',format='png')
# plt.close() 

#save_plot(df_coronary_filtered,'coronary')

# q0 = df_coronary[df_coronary['quartiles']==0].mean(axis=0)
# q1 = df_coronary[df_coronary['quartiles']==1].mean(axis=0)
# q2 = df_coronary[df_coronary['quartiles']==2].mean(axis=0)
# q3 = df_coronary[df_coronary['quartiles']==3].mean(axis=0)

# fig, ax = plt.subplots()
# plt.bar(x=q3.index[:-2],height=q3[:-2], color='green', label=f'3 - 60 samples')
# plt.bar(x=q2.index[:-2],height=q2[:-2], color='orange', label=f'2 - 60 samples')
# plt.bar(x=q1.index[:-2],height=q1[:-2], color='blue', label=f'1 - 58 samples')
# plt.bar(x=q0.index[:-2],height=q0[:-2], color='red', label=f'0 - 62 samples')
# ax.set_xticklabels([])
# plt.title('Mean expression of each necrop gene for all samples - Coronary')
# plt.xlabel(f'{len(q0)} genes from necroptosis pathway')
# plt.ylabel(r'Mean expression in $log_2{(tpm+1)}$')
# plt.legend()
# plt.show()