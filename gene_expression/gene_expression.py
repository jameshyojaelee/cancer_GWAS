import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects

import seaborn as sns
sns.set_style('darkgrid')
sns.set_palette('muted')
sns.set_context("notebook", font_scale=1.5,
                rc={"lines.linewidth": 2.5})
RS = 123

exp = pd.read_csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_expression.csv") 
exp.shape
#1376 cell lines and 19177 genes
exp.head()
exp = exp.rename(columns = {'Unnamed: 0':'DepMap_ID'})

#import inhibition dataset
inhibition = pd.read_csv("C:/Users/james/Desktop/Sumoylation_Analysis/project/inhibition_analysis/inhibition_data.csv")
inhibition = inhibition[['cell', 'AUC', 'IC50']]
inhibition.shape
#delete all whitespace
inhibition['cell'] = inhibition['cell'].str.replace(' ', '')
#capitalize all letters
inhibition['cell'] = inhibition['cell'].str.upper()
ccl = pd.read_csv("C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv")
ccl = ccl[['DepMap_ID', 'stripped_cell_line_name']]
ccl = ccl.rename(columns = {'stripped_cell_line_name':'cell'})
ccl.head()
in2 = pd.merge(ccl, inhibition, on = 'cell')
in2.shape
in2.head()

final_df = pd.merge(in2, exp, on= 'DepMap_ID')
final_df.shape
#COLO205 (ACH-001039) is missing in the expression dataset

y = final_df[['AUC']]
y['AUC'] = pd.qcut(y.AUC, 10, labels=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
y = y.astype(int)
x = final_df.drop(['DepMap_ID', 'cell', 'AUC'], axis=1)
x.shape

from sklearn.decomposition import PCA
pca = PCA(n_components=4)
pca_result = pca.fit_transform(x)
pca.explained_variance_ratio_
# array([0.16234262, 0.12325893, 0.06339874, 0.05216982, 0.04676306, 0.04323975, 0.03880136, 0.03493272, 0.03401673, 0.03058193])
# The first 4 components take up to about 40% of the variance 

pca_df = pd.DataFrame(columns = ['pca1','pca2','pca3','pca4'])
pca_df['pca1'] = pca_result[:,0]
pca_df['pca2'] = pca_result[:,1]
pca_df['pca3'] = pca_result[:,2]
pca_df['pca4'] = pca_result[:,3]
pca_df

group = np.array(y)



plt.scatter(pca_df.iloc[:, 0], pca_df.iloc[:, 1],
            c= y, edgecolor='none', alpha=0.5,
            cmap=plt.cm.get_cmap('Accent', 10))
plt.xlabel('component 1')
plt.ylabel('component 2')
plt.show()



# Utility function to visualize the outputs of PCA and t-SNE
def fashion_scatter(x, colors):
    # choose a color palette with seaborn.
    num_classes = len(np.unique(colors))
    palette = np.array(sns.color_palette("hls", num_classes))

    # create a scatter plot.
    f = plt.figure(figsize=(100, 100))
    ax = plt.subplot(aspect='equal')
    sc = ax.scatter(x[:,0], x[:,1], lw=0, s=40, c=palette)
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    ax.axis('off')
    ax.axis('tight')

    # add the labels for each digit corresponding to the label
    txts = []
    for i in range(num_classes):
        # Position of each label at median of data points.
        xtext, ytext = np.median(x[colors == i, :], axis=0)
        txt = ax.text(xtext, ytext, str(i), fontsize=24)
        txt.set_path_effects([
            PathEffects.Stroke(linewidth=5, foreground="w"),
            PathEffects.Normal()])
        txts.append(txt)
    return f, ax, sc, txts

fashion_scatter(pca_df.values, y)