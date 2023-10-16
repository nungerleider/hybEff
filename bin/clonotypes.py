import glob
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import numpy as np
import matplotlib.pyplot as plt

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

b55 = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/bl_cpm_55ebvpos.tsv', index_col=0)
b55.columns = b55.columns.map(lambda x:x.replace('.', '-'))
vdj2 = vdj.loc[b55.columns]
bl = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/bl_gene_counts_barcode.tsv', index_col=0)
cib = pd.read_csv('/Users/nate/Projects/EBV_interactome/bl/CIBERSORTx_Job187_Results.csv', index_col=0)
wgs = pd.read_table('/Users/nate/Projects/EBV_interactome/wgs_bl_tumors.tsv', index_col=0,sep=',') 

cl = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/blood871418-suppl2.csv', sep=',')
cl = cl.set_index('Patient barcode')

discovery = set(cl[cl['Cohort']=='Discovery'].index)
ebv_pos = set(bl.loc['RPMS1'].sort_values().index[18:])
mir = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/bl_mir_cpm_primary_tumors_only.tsv', index_col=0)
mirs = set(mir.columns)

colnames = list(discovery & ebv_pos & mirs)

bl = bl[colnames]
# bl = 1000000*bl/np.sum(bl)
mir = mir[colnames]
ebv = np.sum(mir.loc[[i for i in mir.index if 'ebv' in i]])

#both = set(wgs.index) & set(bl.columns)

vdj2['mapped'] = np.sum(bl)
vdj2['unique_norm'] = 1000000*(vdj2['unique'] / vdj2['mapped'])
vdj2['mir'] = ebv

cibv = cib.loc[vdj2.index] 

fig = plt.figure(figsize=(8,8))
ax = plt.subplot()

ax.scatter(vdj2['mir'], vdj2['unique'], s=50, c='orange')
#plt.scatter(vdj2['mir'], cibv['T cells CD8'],s=50,c='orange')

m, b = np.polyfit(vdj2['mir'], vdj2['unique'], 1) 
plt.plot(vdj2['mir'],vdj2['mir']*m  + b, c='k')
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/bl_clonotype_v_mir.svg')
print("ebv mir vs cd3 clonotypes", spearmanr(vdj2['mir'], vdj2['unique']))


# for i in cibv.columns:
#     print(i,spearmanr( cibv[i], vdj2['unique_norm']))
# for i in cibv.columns:
#     print(i,pearsonr( cibv[i], vdj2['unique_norm']))


###STAD

files = glob.glob('/Users/nate/Projects/EBV_interactome/stad/vdj/*vdj.tsv')
vlen = {}
vcount= {}
for i in files:
    x = pd.read_table(i, index_col=0)
    barcode = i.split('/')[-1].split('.')[0]
    x = x[x['v'].str.startswith('T')]
    vlen[barcode] = len(x.index)
    vcount[barcode] = np.sum(x.index)
    print(x.iloc[:3])
vdj= pd.DataFrame.from_dict(vlen, orient='index')
vdj.columns = ['unique']
vdj['counts'] = pd.DataFrame.from_dict(vcount,orient='index')

mir = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stadmirs_counts211neg24pos.fordeseq.tsv', index_col=0)
mir= 1000000*mir/np.sum(mir)
sebv = np.sum(mir.loc[[i for i in mir.index if 'ebv' in i]])
sebv = sebv[-24:]
scib = pd.read_csv('/Users/nate/Projects/EBV_interactome/stad/CIBERSORTx_Job190_Results.csv', index_col=0)

stad = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stad_mss_only_211neg_24pos.fordeseq.tsv', index_col=0)
scibe=scib.loc[sebv.index]
vdj.index = [i.rsplit('-',1)[0] for i in vdj.index] 
vdj2 = vdj.loc[sebv.index]
vdj2['mir'] = sebv
vdj2['mapped']=np.sum(stad) 
vdj2['unique_norm'] = 1000000*(vdj2['unique'] / vdj2['mapped'])

# for i in scib.columns:
#     print(i,pearsonr( scibe[i], sebv))
# for i in cib.columns:
#     print(i,pearsonr( cibv[i], vdj2['unique_norm']))

fig = plt.figure(figsize=(8,8))
ax = plt.subplot()

plt.scatter(vdj2['mir'], vdj2['unique'], s=50, c='orange')
#plt.scatter(vdj2['mir'], cibv['T cells CD8'],s=50,c='orange')

m, b = np.polyfit(vdj2['mir'], vdj2['unique'], 1) 
plt.plot(vdj2['mir'],vdj2['mir']*m  + b, c='k')
plt.savefig('/Users/nate/Projects/EBV_interactome/figures/stad_clonotype_v_mir.svg')
print("ebv mir vs cd3 clonotypes", spearmanr(vdj2['mir'], vdj2['unique']))




# EBV pos v neg unique clonotypes BL


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

#12 neg the rest pos
vdj = vdj.loc[b.columns]

#first sample looks like an outlier

neg_mean = np.mean(vdj.loc[b.columns[:12],'unique'])
pos_mean = np.mean(vdj.loc[b.columns[12:], 'unique'])
neg_sem = sem(vdj.loc[b.columns[:12],'unique'])
pos_sem = sem(vdj.loc[b.columns[12:], 'unique'])

print(mannwhitneyu(vdj.loc[b.columns[:12],'unique'], vdj.loc[b.columns[12:], 'unique']))

fig, ax = plt.subplots(figsize=(6,8))
plt.bar([0], neg_mean, yerr=[[0],[neg_sem]])
plt.bar([1], pos_mean, yerr=[[0],[pos_sem]])

plt.savefig('/Users/nate/Projects/EBV_interactome/bl/pos_v_neg_bl_clonotypes.svg')






# EBV pos v neg unique clonotype stad
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
vdj.index = [i[:12] for i in vdj.index] 

stad = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stad_mss_only_211neg_24pos.fordeseq.tsv', index_col=0)

vdj = vdj.loc[stad.columns]

neg_mean = np.mean(vdj.loc[stad.columns[:-24],'unique'])
pos_mean = np.mean(vdj.loc[stad.columns[-24:], 'unique'])
neg_sem = sem(vdj.loc[stad.columns[:-24],'unique'])
pos_sem = sem(vdj.loc[stad.columns[-24:], 'unique'])

fig, ax = plt.subplots(figsize=(6,8))
plt.bar([0], neg_mean, yerr=[[0],[neg_sem]])
plt.bar([1], pos_mean, yerr=[[0],[pos_sem]])

print('stad', mannwhitneyu(vdj.loc[stad.columns[:-24],'unique'], vdj.loc[stad.columns[-24:], 'unique']))

plt.savefig('/Users/nate/Projects/EBV_interactome/stad/pos_v_neg_stad_clonotypes.svg')