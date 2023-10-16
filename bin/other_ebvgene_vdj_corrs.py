import pandas as pd
import glob
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

#bl
files = glob.glob('/Users/nate/Projects/EBV_interactome/bl/vdj/*vdj.tsv')
vlen = {}
vcount= {}
for i in files:
    x = pd.read_table(i, index_col=0)
    barcode = i.split('/')[-1].split('.')[0]
    x = x[x['v'].str.startswith('T')]
    vlen[barcode] = len(x.index)
    vcount[barcode] = np.sum(x.index)
vdj= pd.DataFrame.from_dict(vlen, orient='index')
vdj.columns = ['unique']
vdj['counts'] = pd.DataFrame.from_dict(vcount,orient='index')

b = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/revised_bl_ebvneg14ebvpos59_fordeseq.tsv', index_col=0)
ffpe = ['BLGSP-71-22-00339', 'BLGSP-71-22-00347']
b.columns = [i.replace('.', '-') for i in b.columns]
b = b[[i for i in b.columns if i not in ffpe]]
ebv_pos = b[b.columns[12:]]
mir = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/bl_mir_cpm_primary_tumors_only.tsv', index_col=0)
mir = mir.loc[[i for i in mir.index if 'ebv' in i]]
both = set(ebv_pos.columns) & set(mir.columns)
ebv_pos = ebv_pos[both]
mir = mir[both]

vdj = vdj.loc[ebv_pos.columns]

ebvgenes = pd.read_table('/Users/nate/Projects/EBV_interactome/ebvgenes',header=None)
ebvgenes = set(ebvgenes[0]) & set(ebv_pos.index)
ebvgenes = [i for i in ebvgenes if 'miR' not in i]

genes = []
for i in ebvgenes:
    if i in ebv_pos.index:
        cor = spearmanr(ebv_pos.loc[i], vdj['unique'])[0]
        genes.append((cor,i))
mir = mir.loc[np.mean(mir,1) > 100]
for i in mir.index:
    cor = spearmanr(mir.loc[i], vdj['unique'])[0]
    genes.append((cor,i))
genes.sort(key=lambda x:x[0])

y = range(len(genes))
x= [i[0] for i in genes]

c = ['r' if 'ebv' in i[1] else 'y' for i in genes]
fig = plt.figure(figsize=(8,6))
ax = plt.subplot()
plt.scatter(x,y,c=c,alpha=0.5,lw=0,s=100) 
plt.savefig('/Users/nate/Projects/EBV_interactome/bl/bl_vdj_corr_with_all_ebv_genes_and_mirs.svg')
plt.close()

#stad

files = glob.glob('/Users/nate/Projects/EBV_interactome/stad/vdj/*vdj.tsv')
vlen = {}
vcount= {}
for i in files:
    x = pd.read_table(i, index_col=0)
    barcode = i.split('/')[-1].split('.')[0]
    x = x[x['v'].str.startswith('T')]
    vlen[barcode] = len(x.index)
    vcount[barcode] = np.sum(x.index)
vdj= pd.DataFrame.from_dict(vlen, orient='index')
vdj.columns = ['unique']
vdj['counts'] = pd.DataFrame.from_dict(vcount,orient='index')

mir = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stad_mirs.tsv', index_col=0)
mir = mir.loc[[i for i in mir.index if 'ebv' in i]]
stad = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stad_mss_only_211neg_24pos.fordeseq.tsv', index_col=0)
st = pd.read_table('stad_gene_est_counts.tsv', index_col=0)   
st.columns = [i[:12] for i in st.columns]                                                      
st=st[stad.columns[-24:]]                                                                      
mir = mir[stad.columns[-24:]]
mir = mir.loc[np.mean(mir,1) > 100]
vdj.index = [i.rsplit('-',1)[0] for i in vdj.index] 
vdj=vdj.loc[mir.columns]

genes = []
for i in ebvgenes:
    if i in st.index:
        cor = spearmanr(st.loc[i], vdj['unique'])[0]
        if cor == cor:
            genes.append((cor,i))
for i in mir.index:
    cor = spearmanr(mir.loc[i], vdj['unique'])[0]
    genes.append((cor,i))
genes.sort(key=lambda x:x[0])


y = range(len(genes))
x= [i[0] for i in genes]

fig = plt.figure(figsize=(8,6))
ax = plt.subplot()
c = ['r' if 'ebv' in i[1] else 'y' for i in genes]
plt.scatter(x,y,c=c,alpha=.5,lw=0, s=100) 
plt.savefig('/Users/nate/Projects/EBV_interactome/stad/stad_vdj_corr_with_all_ebv_genes_and_mirs.svg')
