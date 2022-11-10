import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

f1 = pd.read_csv('expression_genes.csv',index_col=[0])
print(f1)
f1 = f1.applymap(np.log2)
print(f1.min(),f1.max())
cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)
sns.heatmap(f1, cmap=cmap, square=True, linewidth=0.3, cbar_kws={"shrink": .8}, vmin=-7, vmax=6,mask=f1.isnull())
plt.savefig('genes_log2.pdf')
# plt.show()
