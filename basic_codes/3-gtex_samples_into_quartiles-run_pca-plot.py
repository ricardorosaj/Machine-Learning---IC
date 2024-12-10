import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import util 
import os
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

"""
Title: Description
Applies PCA or tSNE on expression data (for necrop, protein coding or all genes) for dimensionality reduction and plots it
"""

def transform_data(df):
    df.columns=df.iloc[0]
    df = df[1:]
    df = df.reset_index()
    df.to_csv('temporario.csv', index=False)
    df = pd.read_csv('temporario.csv')
    df = df.rename(columns={'index':'sample_id'})
    df = df.set_index('sample_id')
    return df

def transform_data_pca(df, pca_model):
    pca_list = []
    list_index = []
    for i in range(len(df)):
        line = df.iloc[i]
        pca_list.append(pca_model.transform([line])[0])
        list_index.append(line.name) #saves names of genes for indexing with quarties and zscores
    return pd.DataFrame(np.array(pca_list),index=list_index)

"""
Function: create_df_for_method_with_time_of_death_and_applies_it
 Parameters:
 	df - dataframe with expression of genes
    method - string with values 'pca' or 'tsne' for dimensionality reduction
 	
 Returns:
 	dictionary with dataframe transformed, the mean and std of expression of genes by quartile, the samples names and the number of genes
"""

def create_df_for_method_with_time_of_death_and_applies_it(df, method):
    n_genes = len(df)
    samples = df.columns.values

    time_of_death = pd.read_csv(util.path__+'../../../data/raw/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t', usecols=['SAMPID','SMTSISCH'])

    dict_df = {'SAMPID':[],'SMTSISCH':[]}
    for i in range(len(time_of_death)):
        sample = time_of_death['SAMPID'][i]
        if sample in samples:
            dict_df['SAMPID'].append(sample)
            dict_df['SMTSISCH'].append(time_of_death['SMTSISCH'][i])

    for i in range(len(dict_df['SAMPID'])):
        dict_df['SMTSISCH'][i] = dict_df['SMTSISCH'][i]/60 #from minutes to hours
    dict_df['quartiles'] = pd.qcut(dict_df['SMTSISCH'], q=4, labels=False)

    time_of_death_df = pd.DataFrame(dict_df,index=dict_df['SAMPID'])
    mean = time_of_death_df.groupby(time_of_death_df['quartiles']).mean()
    mean = mean['SMTSISCH'].values
    std = time_of_death_df.groupby(time_of_death_df['quartiles']).std()
    std = std['SMTSISCH'].values

    df = df.T
    df = transform_data(df)
    index = df.index

    if method=='pca':
        pca = PCA(n_components=2)
        pca.fit(df)
        df = transform_data_pca(df,pca)
        quartile = df.join(time_of_death_df, how='inner')
        return {'quartile':quartile,'mean':mean,'std':std,'samples':samples,'n_genes':n_genes,'pca':pca}
    elif method=='tsne':
        tsne = TSNE(n_components=2)
        df = tsne.fit_transform(df)
        df = pd.DataFrame(df,index=index,columns=[0,1])
        quartile = df.join(time_of_death_df,how='inner')
        return {'quartile':quartile,'mean':mean,'std':std,'samples':samples,'n_genes':n_genes}
    else:
        print('Error: invalid method')    

"""
Function: plot_df_transformed
 Parameters:
 	quartile - dataframe created by 'create_df_for_method_with_time_of_death_and_applies_it' function
    method - string with values 'pca' or 'tsne'
    tissue - which tissue is the data from
 	
 Returns:
 	Plot of tSNE or PCA for the data inputed. Prints error message if method is not 'tsne' or 'pca'.
"""

def plot_df_transformed(quartile,method,tissue):
    df = quartile['quartile']
    mean = quartile['mean']
    std = quartile['std']
    samples = quartile['samples']
    n_genes = quartile['n_genes']
    if method == 'pca':
        pca = quartile['pca']
        sns.scatterplot(x=df[0],y=df[1],hue=df['quartiles'],palette='bright')
        plt.xlabel(f'PCA 1 - {round(pca.explained_variance_ratio_[0]*100,1)}%')
        plt.ylabel(f'PCA 2 - {round(pca.explained_variance_ratio_[1]*100,1)}%')
        plt.title(f'PCA - {tissue} - {len(samples)} samples - 2 components - Quartile - {n_genes} genes')
        plt.show()
    elif method == 'tsne':
        sns.scatterplot(x=df[0],y=df[1],hue=df['quartiles'],palette='bright')
        plt.xlabel('x1')
        plt.ylabel('x2')
        plt.title(f'tSNE {tissue} - {len(samples)} samples - 2 components - Quartile - {n_genes} genes')
        plt.show()
    else:
        print('Error: invalid method')    

# filtered_df_atrial = pd.read_csv(util.path__+'../../../data/interim/gtex_data/filtered_protein_coding_log2_tpm_data_atrial_appendage.csv'.replace('/',os.sep))
# filtered_df_ventricle = pd.read_csv(util.path__+'../../../data/interim/gtex_data/filtered_protein_coding_log2_tpm_data_left_ventricle.csv'.replace('/',os.sep))

# df_atrial = pd.read_csv(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_atrial_appendage_gtex.csv'.replace('/',os.sep))
# df_ventricle = pd.read_csv(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_left_ventricle_gtex.csv'.replace('/',os.sep))

# filtered_atrial = create_df_for_method_with_time_of_death_and_applies_it(filtered_df_atrial,'pca')
# plot_df_transformed(filtered_atrial,'pca','Atrial Appendage')
# atrial = create_df_for_method_with_time_of_death_and_applies_it(df_atrial,'tsne')
# plot_df_transformed(atrial,'tsne','Atrial Appendage')

# filtered_ventricle = create_df_for_method_with_time_of_death_and_applies_it(filtered_df_ventricle,'pca')
# plot_df_transformed(filtered_ventricle,'pca','Left Ventricle')
# ventricle = create_df_for_method_with_time_of_death_and_applies_it(df_ventricle,'tsne')
# plot_df_transformed(ventricle,'tsne','Left Ventricle')

# df_coronary = pd.read_csv(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_coronary_gtex.csv'.replace('/',os.sep), index_col=['Gene_ID'])
# filtered_df_coronary = pd.read_csv(util.path__+'../../../data/interim/gtex_data/filtered_necrop_log2_tpm_coronary.csv'.replace('/',os.sep), index_col=['Gene_ID'])
# atrial = create_df_for_method_with_time_of_death_and_applies_it(df_atrial,'pca')
# plot_df_transformed(atrial,'pca','Atrial Appendage')

# filtered_coronary = create_df_for_method_with_time_of_death_and_applies_it(filtered_df_coronary,'tsne')
# plot_df_transformed(filtered_coronary,'tsne','Coronary')

# coronary = create_df_for_method_with_time_of_death_and_applies_it(df_coronary,'tsne')
# plot_df_transformed(coronary,'tsne','Coronary')

df_liver = pd.read_csv(util.path__+'../../../data/interim/gtex_data/log2_tpm_data_Liver_gtex.csv'.replace('/',os.sep), index_col=['Gene_ID'])
liver = create_df_for_method_with_time_of_death_and_applies_it(df_liver, 'tsne')
plot_df_transformed(liver, 'tsne','Liver')
df_necrop_liver = pd.read_csv(util.path__+'../../../data/interim/gtex_data/filtered_necrop_log2_tpm_data_Liver_gtex.csv'.replace('/',os.sep), index_col=['Gene_ID'])
necrop_liver = create_df_for_method_with_time_of_death_and_applies_it(df_necrop_liver, 'tsne')
plot_df_transformed(necrop_liver, 'tsne','Liver')
df_pc_liver = pd.read_csv(util.path__+'../../../data/interim/gtex_data/filtered_protein_coding_log2_tpm_data_Liver_gtex.csv'.replace('/',os.sep), index_col=['Gene_ID'])
pc_liver = create_df_for_method_with_time_of_death_and_applies_it(df_pc_liver, 'tsne')
plot_df_transformed(pc_liver, 'tsne','Liver')