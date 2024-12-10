Codes in this folder execute the basic procedures for extracting samples from GTEx data, filtering the genes that will be considered, such as necroptosis genes and protein coding, and inputing this data to the decision tree model, saving all information in the SQL database.
Names of codes are in order of execution.

After running code 3, analyse which quartiles cluster and use this information to erase quartiles in code 4.

NOTE: processed data should have a final column with the tissue classes (Example: CVD+A, CVD, C+A, C), created by reading info of death causes, medication, medical report, etc and deciding which sample goes to which group. No code creates this column.
