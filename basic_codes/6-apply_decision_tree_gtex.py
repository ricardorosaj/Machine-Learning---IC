import pandas as pd
from sklearn import tree
import graphviz
import numpy as np
import random
import sys 
import sqlite3
import os
import util

"""
Title:Descrition
Main code which runs a desicion tree algorithm on GTEx data and saves information on database.
Functions used in this script are the same used on scripts in this folder and on /decision_tree_iguaracy scripts.
"""

"""
Function: change_position_randomly

Description: Changes positions of columns randomly to avoid bias on algorithm.

Parameters: 
	condition - column of a dataframe.

Returns:
	A list of the expression values randomly shuffled.
"""
def change_position_randomly(condition):
    list_ = []
    for i in range(len(condition)):
        list_.append(condition[i])
    for i in range(100):
        rand_num1 = random.randint(0,len(condition)-1)
        rand_num2 = random.randint(0,len(condition)-1)
        var = 0
        var = list_[rand_num1]
        list_[rand_num1] = list_[rand_num2]
        list_[rand_num2] = var
    return list_

"""
Function: create_x_y_from_data_frame

 Description:
    This function creates the x and y for decision tree algorithm with decoder of positions to name and vice-versa.

 Parameters:
 	df - dataframe with expression data and classes of samples.
 	
  Returns:
 	Dictionary with x and y values to use as input of algorithm and decoders linking positions of values in list and vice-versa. 
"""
def create_x_y_from_data_frame(df):
    
    x = df[df.columns[:-1]] #Values of expression excluding the last column, containing the conditions (targets) 
    x_copy = np.array(x)
    y = df['Pacient_Condition'] #Targets
    y_copy = np.array(y)
    list_of_unique_y = y.unique()
    decoderNamesToPositionsOfValues = {}
    decoderPositionsOfValuesToNames = []
    size_of_list = len(list_of_unique_y)
    for i in range(size_of_list):
        decoderNamesToPositionsOfValues[list_of_unique_y[i]]=i #Saves which position is the value
        decoderPositionsOfValuesToNames.append(list_of_unique_y[i]) #Saves which value is linked to a given position

    y_class_are_numbers= []
    for i in range(len(y_copy)):
        y_class_are_numbers.append(decoderNamesToPositionsOfValues[y_copy[i]])

    return {'x':x_copy, 
            'y':np.array(y_class_are_numbers),
            'decoderPositionsOfValuesToNames': decoderPositionsOfValuesToNames,
            'decoderNamesToPositionsOfValues':decoderNamesToPositionsOfValues}

"""
Function: get_lineage

Description: Function to print a decision tree using information from sklearn.tree. It's useful for saving information to the database.

 Parameters:
 	tree - decision tree model fitted
 	feature_names - name of genes
 	
  Returns:
 	List containing useful information of each node in the generated tree
"""
def get_lineage(tree, feature_names):
    left      = tree.tree_.children_left
    right     = tree.tree_.children_right
    threshold = tree.tree_.threshold
    features  = [feature_names[i] for i in tree.tree_.feature]
    values = tree.tree_.value #values of all vectors of conditions
    idx = np.argwhere(left == -1)[:,0] # get ids of child nodes 
    def recurse(left, right, child, lineage=None):          
        if lineage is None:
            lineage = [(child,)]
        if child in left:
            parent = np.where(left == child)[0].item()
            split = 'l'
        else:
            parent = np.where(right == child)[0].item()
            split = 'r'
        lineage.append((parent, split, threshold[parent], features[parent]))
        if parent == 0:
            lineage.reverse()
            return lineage
        else:
            return recurse(left, right, parent, lineage)
    global_nodes_list = []          
    for child in idx:
        nodes_list = []
        for node in recurse(left, right, child):
            if(len(node)==1):
                #print(values[node[0]])
                nodes_list.append(values[node[0]])
            #else:
                #print(node)
            nodes_list.append(node)
        global_nodes_list.append(nodes_list)    
    return global_nodes_list

"""
Function: print_tree_in_pdf

Description: Saves ther tree generated in the iteration to a pdf file.

 Parameters:
 	clf - decision tree model fitted
 	genes_symbols - name of genes
 	
  Returns:
 	PDF render of generated tree.
"""

def print_tree_in_pdf(clf, genes_symbols):
    dot_data = tree.export_graphviz(clf,  feature_names=genes_symbols, out_file=None)
    graph = graphviz.Source(dot_data)
    print_tree = graph.render(util.path__+'../../../results/gtex_results/test_gtex_tree')
    return print_tree

df = pd.read_csv(util.path__+'../../../../data/processed/gtex_data/coronary_gtex_log2_tpm_pacient_condition_only_protein_coding_genes_quartiles_1_2.csv'.replace('/',os.sep), index_col=['Unnamed: 0']) #Data read
importance_for_features={}
for feature in list(df.columns[:-1]):
    importance_for_features[feature]=[]

db = sqlite3.connect(util.path__+'../../../../data/interim/sql/bd_ic_v2.sqlite3') #Conn to database
cur = db.cursor()

exper = "Teste"
cur.execute(f"""INSERT INTO experiencia (description) VALUES ({exper});""") #Initiate an experiment
exp = cur.lastrowid
db.commit()

genes_name = list(df.columns)[:-1]
print("")
n_times_run = 1 #Number of trees generated
for iteration in range(n_times_run):

    cur.execute("""INSERT INTO tree (experiencia_id) VALUES ({exp});""".format(exp=exp)) #Insert tree generation at tree    
    cursor_id = cur.lastrowid
    db.commit()

    new_postion_for_columns=change_position_randomly(genes_name) #Randomizes the positions of columns
    df=df[new_postion_for_columns+list([df.columns[-1]]) ] #Creates new df with these new positions
    
    inputs_for_fit=create_x_y_from_data_frame(df) #Creates x and y inputs
    clf = tree.DecisionTreeClassifier()
    clf = clf.fit(inputs_for_fit['x'], inputs_for_fit['y']) #Apply classifier
    
    importancia = clf.feature_importances_
    for j in range(len(importancia)):
        importance_for_features[df.columns[j]].append(importancia[j])

    print_tree_in_pdf(clf, genes_name) #Saves figure
    
    #insert data to place_of_genes_in_trees
    global_tree_nodes = get_lineage(clf, genes_name)
    dict_for_repeated_genes = {}
    values_of_cond = clf.tree_.value
    for i in range(len(global_tree_nodes)):
        for j in range(len(global_tree_nodes[i])-2):
            if not(global_tree_nodes[i][j][3] in dict_for_repeated_genes):
                dict_for_repeated_genes[global_tree_nodes[i][j][3]]=[]
            dict_for_repeated_genes[global_tree_nodes[i][j][3]].append({'deep':j ,'vectorConditions': values_of_cond[global_tree_nodes[i][j][0]]})
    for key in dict_for_repeated_genes.keys():
        for i in range(len(dict_for_repeated_genes[key])):
            depth = dict_for_repeated_genes[key][i]['deep']
            cvda=int(dict_for_repeated_genes[key][i]['vectorConditions'][0][0])
            ca=int(dict_for_repeated_genes[key][i]['vectorConditions'][0][1])
            c=int(dict_for_repeated_genes[key][i]['vectorConditions'][0][2])
            cvd=int(dict_for_repeated_genes[key][i]['vectorConditions'][0][3])
            sql = """INSERT INTO place_of_genes_in_tree (gene_name,tree_generation,experiencia_id,depth,
CVDA,CA,C,CVD) VALUES ("{id}",{tree},{exp},{dep},{cvda},{ca},{c},{cvd});""".format(id=key,tree=cursor_id,exp=exp,dep=depth,cvda=cvda,ca=ca,c=c,cvd=cvd)
            cur.execute(sql)  
        db.commit()

    #insert data to father_and_son_nodes    
    dict_father_and_son={}
    dict_father_and_son_cond={}
    for i in range(len(global_tree_nodes)):
        for j in range(len(global_tree_nodes[i])-2):
            if not(global_tree_nodes[i][j][3] in dict_father_and_son):
                dict_father_and_son[global_tree_nodes[i][j][3]]=[]
            if j==0:
                cond = global_tree_nodes[i][j][3]+'|'+global_tree_nodes[i][j][3]+'|'+str(j)+'|'+str(values_of_cond[global_tree_nodes[i][j][0]])+'|'+str(values_of_cond[global_tree_nodes[i][j][0]])
                if not(cond in dict_father_and_son_cond):
                    dict_father_and_son[global_tree_nodes[i][j][3]].append({'father':global_tree_nodes[i][j][3],'deep':j,
    'vectorConditions': values_of_cond[global_tree_nodes[i][j][0]],'vectorConditionsFather':values_of_cond[global_tree_nodes[i][j][0]], 'expression_value':global_tree_nodes[i][j][2]})
                    dict_father_and_son_cond[cond] = None
            else:
                cond = global_tree_nodes[i][j][3]+'|'+global_tree_nodes[i][j-1][3]+'|'+str(j)+'|'+str(values_of_cond[global_tree_nodes[i][j][0]])+'|'+str(values_of_cond[global_tree_nodes[i][j-1][0]])
                if not(cond in dict_father_and_son_cond):
                #vectorConditionsFather=dict_father_and_son[global_tree_nodes[i][j-1][3]][-1]['vectorConditions']
                    dict_father_and_son[global_tree_nodes[i][j][3]].append({'father':global_tree_nodes[i][j-1][3],'deep':j,
    'vectorConditions': values_of_cond[global_tree_nodes[i][j][0]],'vectorConditionsFather':values_of_cond[global_tree_nodes[i][j-1][0]], 'expression_value':global_tree_nodes[i][j][2]})
                    dict_father_and_son_cond[cond] = None

    for key in dict_father_and_son.keys():
        for i in range(len(dict_father_and_son[key])):
            depth = dict_father_and_son[key][i]['deep']
            father = dict_father_and_son[key][i]['father']
            expression_value = dict_father_and_son[key][i]['expression_value']
            cvda_father=int(dict_father_and_son[key][i]['vectorConditionsFather'][0][0])
            ca_father=int(dict_father_and_son[key][i]['vectorConditionsFather'][0][1])
            c_father=int(dict_father_and_son[key][i]['vectorConditionsFather'][0][2])
            cvd_father=int(dict_father_and_son[key][i]['vectorConditionsFather'][0][3])
            cvda_son=int(dict_father_and_son[key][i]['vectorConditions'][0][0])
            ca_son=int(dict_father_and_son[key][i]['vectorConditions'][0][1])
            c_son=int(dict_father_and_son[key][i]['vectorConditions'][0][2])
            cvd_son=int(dict_father_and_son[key][i]['vectorConditions'][0][3])
            sql = """INSERT INTO father_and_son_nodes (gene_id_parent,gene_id_son,experiencia_id,depth,
CVDA_father,CA_father,C_father,CVD_father,CVDA_son,CA_son,C_son,CVD_son,expression_value)
VALUES ("{id_f}","{id_s}",{exp},{dep},{cvda_father},{ca_father},{c_father},{cvd_father},{cvda_son},{ca_son},
{c_son},{cvd_son},{expression_value});""".format(id_f=father,id_s=key,exp=exp,dep=depth,cvda_father=cvda_father,ca_father=ca_father,
c_father=c_father,cvd_father=cvd_father,cvda_son=cvda_son,ca_son=ca_son,c_son=c_son,cvd_son=cvd_son,expression_value=expression_value)
            cur.execute(sql)
        db.commit()
            
    #insert data to tree_nodes, esta_em, genes_in_path and path
    for i in range(len(global_tree_nodes)):
        path_id_lastrowid = None 
        for j in range(len(global_tree_nodes[i])):
            j0=len(global_tree_nodes[i])-j-2
            if len(global_tree_nodes[i][j0])!=1:
                gene_id_from_node = global_tree_nodes[i][j0][3]
                sql_insert_on_genes_in_path = """INSERT INTO genes_in_path (gene_id,depth,path_id,experiencia_id)
                VALUES ("{gene_id}",{depth},{path_id},{exp})""".format(gene_id = gene_id_from_node, depth=j0, path_id = path_id_lastrowid, exp=exp)
                cur.execute(sql_insert_on_genes_in_path)
                db.commit()
            elif(j==0):
                list_of_condition = list(global_tree_nodes[i][j0][0])
                condition_decoded = inputs_for_fit['decoderPositionsOfValuesToNames']
                for k in range(len(list_of_condition)):
                    if list_of_condition[k]!=0:
                        condition_to_insert = condition_decoded[k]
                        sql_insert_condition_on_path = """INSERT INTO path (cond,nodes,experiencia_id)
                        VALUES ("{condition}",{nodes},{exp})""".format(condition = condition_to_insert, nodes = len(global_tree_nodes[i])-2,exp=exp)
                        cur.execute(sql_insert_condition_on_path)
                        path_id_lastrowid = cur.lastrowid
                    db.commit()   
        
        for j in range(len(global_tree_nodes[i])):        
            if len(global_tree_nodes[i][j])==1:                
                node=None
                for m in range(j-1):
                    nodeT=global_tree_nodes[i][m]                    
                    idForVariableNode=None

                    sql_search_esta_em = """SELECT gene_id FROM esta_em WHERE tree_id = {tree_id} AND gene_id= "{gene_id}" """.format(tree_id = cursor_id,gene_id=nodeT[3])
                    cur.execute(sql_search_esta_em)
                    id_for_gene_in_esta_em = cur.fetchall() 
                    if len(id_for_gene_in_esta_em)==0:
                            sql_insert_on_esta_em = """INSERT INTO esta_em (tree_id,gene_id,experiencia_id)
                    VALUES ({tree_id}, "{gene_id}",{exp});""".format(tree_id=cursor_id ,gene_id=nodeT[3],exp=exp)
                            cur.execute(sql_insert_on_esta_em)

                    if(node!=None):
                        sql = """SELECT tree_node_id FROM tree_node WHERE gene_id="{gene_id}" AND tree_id = {tree_id};""".format(gene_id=node[3], tree_id=cursor_id)
                        cur.execute(sql)
                        idForVariableNode=cur.fetchall()[0][0]
                        sql_insert_on_tree_node = """INSERT INTO tree_node (depth,father_node_id,tree_id,gene_id,experiencia_id)
                    VALUES ({depth}, {father_id}, {tree_id}, "{gene_id}",{exp});""".format(depth=m, father_id=idForVariableNode, tree_id=cursor_id ,gene_id=nodeT[3],exp=exp)
                        cur.execute(sql_insert_on_tree_node)

                    else:
                        sql_insert_on_tree_node = """INSERT INTO tree_node (depth, tree_id, gene_id,experiencia_id) VALUES ({depth}, {tree_id}, "{gene_id}",{exp});""".format(depth=m, tree_id=cursor_id ,gene_id=nodeT[3],exp=exp)
                        cur.execute(sql_insert_on_tree_node)
                        
                    node=nodeT
    sys.stdout.write('tree generated:' +str(iteration+1) +':' +str(n_times_run)+'\r') 
print("")

db.commit()
cur.close()
db.close()
