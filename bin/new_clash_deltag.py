import pandas as pd
import numpy as np
from scipy.stats import ttest_rel, sem, ks_2samp
import glob
import matplotlib.pyplot as plt


snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')
akata_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*annotated')
mhv68_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*annotated')
kshv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*annotated')


snu_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*shuffled_dg.tsv')
akata_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*shuffled_dg.tsv')
mhv68_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*shuffled_dg.tsv')
kshv_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*shuffled_dg.tsv')



mhv68_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/smallfrac/CL100128477_L1_WHMOUkmcSAAETAASE-545_1.fq.counts.tsv']
mhv68_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/longfrac/HE2-1.genes.tsv']

kshv_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/smallfrac/541_WT-TIVE-1.small_fraction.counts.tsv']
kshv_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/longfrac/TIVE-1.genes.tsv']

akata_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/smallfrac/*')
akata_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/akata/longfrac/Akata.genes.tsv']

snu_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/smallfrac/*')
snu_long_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/longfrac/*')



def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    return df


def calculate_average_binding_energy_by_species_unweighted(path, vir='ebv', min_counts=10):

    df = process_df(path, vir=vir)
    virus = df[df['species'] == vir].groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
    virus_filtered = virus[virus['count'] >= min_counts]
    host = df[df['species'] == 'host'].groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
    host_filtered = host[host['count'] >= min_counts]

    return np.array(virus_filtered['binding_energy']), np.array(host_filtered['binding_energy'])



colors = ['#57E2E5', '#57E2E5', '#57E2E5', '#E08DAC', '#E08DAC', '#E08DAC']
min_counts = 100
vir = 'ebv'
earr, harr, eshuff, hshuff = [], [], [], []
path = snu_mrna_clash_paths[0]

# for path in snu_shuffled_clash_paths:
path = snu_shuffled_clash_paths[0]
df = pd.read_csv(path, sep='\t',index_col=0)
blacklist = set([i for i in df['mrna'] if 'MT-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
df = df[~df['mrna'].isin(blacklist)]
df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
df['count'] = 1000000 * df['count'] / np.sum(df['count'])
df = df[df['count']>min_counts]
es = np.array(df[df['species'] == vir]['shuffle_deltaG'])
hs = np.array(df[df['species'] == 'host']['shuffle_deltaG'])
# eshuff.append(np.mean(es))
# hshuff.append(np.mean(hs))
    
e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='ebv', min_counts=min_counts)
print(ks_2samp(e,h))

    # break
    # print('ebv:', np.mean(e), '    hsa: ', np.mean(h))
    # earr.append(np.mean(e))
    # harr.append(np.mean(h))

# print(ttest_rel(earr, harr))

fig = plt.figure(figsize=(2, 6), dpi=300)
ax = plt.subplot(111)
plt.bar([0,1,2,3], [np.median(hs), np.median(es), np.median(h), np.median(e)], color=['.9','.75',colors[-1], colors[0]], yerr=[[sem(hshuff), sem(eshuff), sem(harr), sem(earr)], [0,0,0 ,0]],ec='k', lw=1)
ax.set_xticks([])
ax.grid(alpha=.3,ls='--')
ax.set_axisbelow(True)
plt.savefig('snu719_deltaG_c.svg')



colors = ['#57E2E5', '#57E2E5', '#57E2E5', '#E08DAC', '#E08DAC', '#E08DAC']
min_counts = 30
vir = 'mghv'
earr, harr, eshuff, hshuff = [], [], [], []

# for path in snu_shuffled_clash_paths:
path = mhv68_shuffled_clash_paths[0]
df = pd.read_csv(path, sep='\t',index_col=0)
blacklist = set([i for i in df['mrna'] if 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
df = df[~df['mrna'].isin(blacklist)]
df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
df['count'] = 1000000 * df['count'] / np.sum(df['count'])
df = df[df['count']>min_counts]
es = np.array(df[df['species'] == vir]['shuffle_deltaG'])
hs = np.array(df[df['species'] == 'host']['shuffle_deltaG'])
# eshuff.append(np.mean(es))
# hshuff.append(np.mean(hs))
    
e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='mghv', min_counts=min_counts)
print(ks_2samp(e,h))

    # break
    # print('ebv:', np.mean(e), '    hsa: ', np.mean(h))
    # earr.append(np.mean(e))
    # harr.append(np.mean(h))

# print(ttest_rel(earr, harr))

fig = plt.figure(figsize=(2, 6), dpi=300)
ax = plt.subplot(111)
plt.bar([0,1,2,3], [np.median(hs), np.median(es), np.median(h), np.median(e)], color=['.9','.75',colors[-1], colors[0]], yerr=[[sem(hshuff), sem(eshuff), sem(harr), sem(earr)], [0,0,0 ,0]],ec='k', lw=1)
ax.set_xticks([])
ax.grid(alpha=.3,ls='--')
ax.set_axisbelow(True)
plt.savefig('mhv68_deltaG_c.svg')



colors = ['#57E2E5', '#57E2E5', '#57E2E5', '#E08DAC', '#E08DAC', '#E08DAC']
min_counts = 30
vir = 'kshv'
earr, harr, eshuff, hshuff = [], [], [], []

# for path in snu_shuffled_clash_paths:
path = kshv_shuffled_clash_paths[0]
df = pd.read_csv(path, sep='\t',index_col=0)
blacklist = set([i for i in df['mrna'] if 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
df = df[~df['mrna'].isin(blacklist)]
df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
df['count'] = 1000000 * df['count'] / np.sum(df['count'])
df = df[df['count']>min_counts]
es = np.array(df[df['species'] == vir]['shuffle_deltaG'])
hs = np.array(df[df['species'] == 'host']['shuffle_deltaG'])
# eshuff.append(np.mean(es))
# hshuff.append(np.mean(hs))
    
e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='kshv', min_counts=min_counts)
print(ks_2samp(e,h))

    # break
    # print('ebv:', np.mean(e), '    hsa: ', np.mean(h))
    # earr.append(np.mean(e))
    # harr.append(np.mean(h))

# print(ttest_rel(earr, harr))

fig = plt.figure(figsize=(2, 6), dpi=300)
ax = plt.subplot(111)
plt.bar([0,1,2,3], [np.median(hs), np.median(es), np.median(h), np.median(e)], color=['.9','.75',colors[-1], colors[0]], yerr=[[sem(hshuff), sem(eshuff), sem(harr), sem(earr)], [0,0,0 ,0]],ec='k', lw=1)
ax.set_xticks([])
ax.grid(alpha=.3,ls='--')
ax.set_axisbelow(True)
plt.savefig('kshv_deltaG_c.svg')


colors = ['#57E2E5', '#57E2E5', '#57E2E5', '#E08DAC', '#E08DAC', '#E08DAC']
min_counts = 30
vir = 'mghv'
earr, harr, eshuff, hshuff = [], [], [], []
for path in mhv68_shuffled_clash_paths:
    df = pd.read_csv(path, sep='\t',index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    df=df[df['count']>min_counts]
    es = np.array(df[df['species'] == vir]['shuffle_deltaG'])
    hs = np.array(df[df['species'] == 'host']['shuffle_deltaG'])
    eshuff.append(np.mean(es))
    hshuff.append(np.mean(hs))
    
for path in mhv68_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='mghv', min_counts=min_counts)
    print('ebv:', np.mean(e), '    hsa: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))

print(ttest_rel(earr, harr))

fig = plt.figure(figsize=(2, 6), dpi=300)
ax = plt.subplot(111)
plt.bar([0,1,2,3], [np.mean(hshuff), np.mean(eshuff), np.mean(harr), np.mean(earr)], color=['.9','.75',colors[-1], colors[0]], yerr=[[sem(hshuff), sem(eshuff), sem(harr), sem(earr)], [0,0,0 ,0]],ec='k', lw=1)
ax.set_xticks([])
ax.grid(alpha=.3,ls='--')
ax.set_axisbelow(True)

plt.savefig('mhv68_deltaG.svg')



colors = ['#57E2E5', '#57E2E5', '#57E2E5', '#E08DAC', '#E08DAC', '#E08DAC']
min_counts = 30
vir = 'kshv'
earr, harr, eshuff, hshuff = [], [], [], []
for path in kshv_shuffled_clash_paths:
    df = pd.read_csv(path, sep='\t',index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    df=df[df['count']>min_counts]
    es = np.array(df[df['species'] == vir]['shuffle_deltaG'])
    hs = np.array(df[df['species'] == 'host']['shuffle_deltaG'])
    eshuff.append(np.mean(es))
    hshuff.append(np.mean(hs))
    
for path in kshv_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='kshv', min_counts=min_counts)
    print('ebv:', np.mean(e), '    hsa: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))

print(ttest_rel(earr, harr))

fig = plt.figure(figsize=(2, 6), dpi=300)
ax = plt.subplot(111)
plt.bar([0,1,2,3], [np.mean(hshuff), np.mean(eshuff), np.mean(harr), np.mean(earr)], color=['.9','.75',colors[-1], colors[0]], yerr=[[sem(hshuff), sem(eshuff), sem(harr), sem(earr)], [0,0,0 ,0]],ec='k', lw=1)
ax.set_xticks([])
ax.grid(alpha=.3,ls='--')
ax.set_axisbelow(True)

plt.savefig('kshv_deltaG.svg')

