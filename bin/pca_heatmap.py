# Figure S1

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



x = pd.read_table('/Users/nate/Projects/EBV_interactome/pca/rsquared_pca_logisticreg_BL_ff_only.tsv', index_col=0)
y = x[['ebvstatus', 'sex', 'age', 'genome', 'translocation', 'isotype', 'site', 'source','variant']]
c = sns.clustermap(y,  row_cluster=False,col_cluster=False, cmap='GnBu',lw=.5,linecolor='.75',square=True, vmin=0, vmax=1)
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/pca_heatmap_bl_ff.svg')



x = pd.read_table('/Users/nate/Projects/EBV_interactome/pca/rsquared_pca_logisticreg_BL.tsv', index_col=0)
y = x[['ebvstatus', 'sex', 'age', 'genome','biopsy', 'translocation', 'isotype', 'site', 'source','variant']]
c = sns.clustermap(y,  row_cluster=False,col_cluster=False, cmap='GnBu',lw=.5,linecolor='.75',square=True,  vmin=0, vmax=1)
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/rsquared_pca_logisticreg_BL.svg')



x = pd.read_table('/Users/nate/Projects/EBV_interactome/pca/bl_only_ebvpos_and_ff.tsv', index_col=0)
y = x[['sex', 'age', 'genome', 'translocation', 'isotype', 'site', 'source','variant']]
c = sns.clustermap(y,  row_cluster=False,col_cluster=False, cmap='GnBu',lw=.5,linecolor='.75',square=True,  vmin=0, vmax=1)
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/pca_heatmap_ebvonly_and_ff.svg')




x = pd.read_table('/Users/nate/Projects/EBV_interactome/pca/stad_pca_r2.tsv', index_col=0)
x=x.iloc[:15]
try:
    x=x.drop('cdkn1a',1)   
except:
    pass  
c = sns.clustermap(x,  row_cluster=False,col_cluster=False, cmap='GnBu',lw=.5,linecolor='.75',square=True,  vmin=0, vmax=1)
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/pca_heatmap_stad_all_genes.svg')

x = pd.read_table('/Users/nate/Projects/EBV_interactome/pca/stad_mss_only.tsv', index_col=0)
x=x.drop('cdkn1a',1) 
x=x.iloc[:15]

c = sns.clustermap(x,  row_cluster=False,col_cluster=False, cmap='GnBu',lw=.5,linecolor='.75',square=True,  vmin=0, vmax=1)
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/pca_heatmap_stad_mss_genes.svg')


x = pd.read_table('/Users/nate/Projects/EBV_interactome/pca/stad_mss_only_ebvpos_only.tsv', index_col=0)
x=x.iloc[:15]

try:
    x=x.drop('cdkn1a',1)     
except:
    pass
c = sns.clustermap(x,  row_cluster=False,col_cluster=False, cmap='GnBu',lw=.5,linecolor='.75',square=True,  vmin=0, vmax=1)
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/pca_heatmap_stad_mss_ebvpos_genes.svg')
