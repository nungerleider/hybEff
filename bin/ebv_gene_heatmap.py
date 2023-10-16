import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


bl_tpm = pd.read_table('/Users/nate/Projects/EBV_interactome/bl_gene_tpms_primary_tumors_only.tsv', index_col=0)
discovery = pd.read_table('discovery_cohort_ebvpos_bl_cpm.tsv', index_col=0)

#bl_tpm.columns = [i.replace('-','.') for i in bl_tpm.columns]

#bl_patients = open('/Users/nate/Projects/EBV_interactomebl/bl_patients').read().strip('\n').split('\t')
bl_tpm = bl_tpm[discovery.columns]


stad_tpm = pd.read_table('/Users/nate/Projects/EBV_interactome/stad_mss_only_211neg24pos_tpms.tsv', index_col=0)
stad_tpm = stad_tpm[stad_tpm.columns[-24:]]

dlbcl_tpm = pd.read_table('/Users/nate/Projects/EBV_interactome/dlbcl_ebv_rna.tsv', index_col=0)

ebv_genes = set(dlbcl_tpm.index) & set(bl_tpm.index) & set(stad_tpm.index)


bl_tpm = bl_tpm.loc[ebv_genes]
dlbcl_tpm = dlbcl_tpm.loc[ebv_genes]
stad_tpm = stad_tpm.loc[ebv_genes]

order = np.mean(bl_tpm,1).sort_values().index
order = [i for i in order if 'EBER' not in i]

bl_tpm = bl_tpm.loc[order]

stad_tpm = stad_tpm.loc[order]
dlbcl_tpm = dlbcl_tpm.loc[order]

dlbcl_plot = sns.clustermap(np.log2(dlbcl_tpm+1),cmap='Greens',method='median',lw=.1,row_cluster=False,linecolor='0.8',vmin=0,vmax=10)
plt.savefig('/Users/nate/Projects/EBV_interactome/dlbcl_ebv_gene_heatmap.svg')
stad_plot = sns.clustermap(np.log2(stad_tpm+1),cmap='Greens',method='median',lw=.1,row_cluster=False,linecolor='0.8',vmin=0,vmax=10)
plt.savefig('/Users/nate/Projects/EBV_interactome/stad_ebv_gene_heatmap.svg')
bl_plot = sns.clustermap(np.log2(bl_tpm+1),cmap='Greens',method='median',lw=.1,row_cluster=False,linecolor='0.8',vmin=0,vmax=10)
plt.savefig('/Users/nate/Projects/EBV_interactome/bl_ebv_gene_heatmap.svg')


bl_mir = pd.read_table('/Users/nate/Projects/EBV_interactome/bl_mir_cpm_primary_tumors_only.tsv', index_col=0)
bl_mir = bl_mir.loc[[i for i in bl_mir.index if 'ebv' in i]]
#l = sorted(bl_mir.index, key=lambda x:int(x[12:].split('-')[0]))
#l = ['ebv-miR-BHRF1-1',
#  'ebv-miR-BHRF1-2-5p',
#  'ebv-miR-BHRF1-2-3p',
#  'ebv-miR-BHRF1-3','ebv-miR-BART1-3p',
#  'ebv-miR-BART1-5p'] + l[6:]
#bl_mir = bl_mir.loc[l]



tpm = set(bl_tpm.columns)
mir = set(bl_mir.columns)
cols = mir & tpm
bl_mir = bl_mir[cols]
bl_tpm = bl_tpm[cols]
c = sns.clustermap(np.log2(bl_mir+1),cmap='Blues',method='ward',lw=.1,linecolor='0.8')
plt.savefig('/Users/nate/Projects/EBV_interactome/bl_ebv_mir_heatmap.svg')

column_order = [bl_mir.columns[i] for i in c.dendrogram_col.reordered_ind ]
row_order = [bl_mir.index[i] for i in c.dendrogram_row.reordered_ind ]

bl_tpm = bl_tpm[column_order]
bl_plot = sns.clustermap(np.log2(bl_tpm+1),cmap='Greens',method='median',lw=.1,row_cluster=False,col_cluster=False,linecolor='0.8',vmin=0,vmax=10)
plt.savefig('/Users/nate/Projects/EBV_interactome/bl_ebv_gene_heatmap.svg')

stad_mir = pd.read_table('/Users/nate/Projects/EBV_interactome/stadmirs_new.ebvposonly.cpm.tsv', index_col=0)
cols = set(stad_mir.columns) & set(stad_tpm.columns)
stad_mir = stad_mir[cols]
stad_mir = stad_mir.loc[[i for i in stad_mir.index if 'ebv' in i]]

ne = set(bl_mir.index) - set(stad_mir.index)
stad_mir=stad_mir.T
for i in ne:
    stad_mir[i] = 0
stad_mir = stad_mir.T


stad_mir = stad_mir.loc[row_order]
stad_tpm = stad_tpm[cols]

c = sns.clustermap(np.log2(stad_mir+1), cmap='Blues', method='ward', lw=.1, linecolor='0.8',row_cluster=False)
plt.savefig('/Users/nate/Projects/EBV_interactome/stad_ebv_mir_heatmap.svg')
column_order = [stad_mir.columns[i] for i in c.dendrogram_col.reordered_ind ]
stad_tpm = stad_tpm[column_order]

# Rerun stomach cancer microRNA because some are missing - use the same index as bl 
stad_plot = sns.clustermap(np.log2(stad_tpm+1),cmap='Greens',method='median',lw=.1,row_cluster=False,col_cluster=False,linecolor='0.8',vmin=0,vmax=10)
plt.savefig('/Users/nate/Projects/EBV_interactome/stad_ebv_gene_heatmap.svg')

