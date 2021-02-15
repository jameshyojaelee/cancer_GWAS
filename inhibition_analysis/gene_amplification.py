import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#import copy number dataset and cell line info dataset
cn = pd.read_csv('C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_gene_cn.csv')
info = pd.read_csv('C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv')
cn.head()
info.head()

#change the column names and subset desired columns
cn = cn.rename(columns = {'Unnamed: 0':'ID'})
cn = cn.rename(columns = {'MYC (4609)':'MYC'})
myc_df = cn[['ID', 'MYC']]
myc_df.head()

info = info.rename(columns = {'DepMap_ID':'ID'})
info = info[['ID', 'lineage', 'stripped_cell_line_name']]
info = info.rename(columns = {'stripped_cell_line_name':'cell'})
info.head()

myc_cn = info.merge(myc_df, on='ID')
myc_cn.head()
myc_cn['MYC'].describe()

#export myc copy number dataset to csv
myc_cn.to_csv("C:/Users/james/Desktop/Sumoylation_Analysis/project/inhibition_analysis/myc_copy_number.csv")

#import inhibition dataset
inhibition = pd.read_csv("C:/Users/james/Desktop/Sumoylation_Analysis/project/inhibition_analysis/inhibition_data.csv")
inhibition = inhibition[['cell', 'AUC']]

#delete all whitespace
inhibition['cell'] = inhibition['cell'].str.replace(' ', '')
#capitalize all letters
inhibition['cell'] = inhibition['cell'].str.upper()
inhibition.head()

final_df = pd.merge(myc_cn, inhibition, on='cell')
final_df.head()
final_df.describe()


#plot AUC vs MYC Copy number
MYC_plot = sns.scatterplot(data=final_df, x="AUC", y="MYC", color='red')
MYC_plot.set(xlabel='AUC', ylabel='MYC copy number in log2(x+1)')

#plt.show()
plt.savefig('C:/Users/james/Desktop/Sumoylation_Analysis/project/inhibition_analysis/MYC_plot.png')
