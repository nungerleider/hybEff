import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr, spearmanr, ks_2samp, ttest_rel, sem
import glob


def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    return df



def integrate(mir_expr_path, mrna_expr_path, clash_path, min_counts=1000):
    
    mir = pd.read_table(mir_expr_path, index_col=0,header=None)    
    mir = mir.drop('*')
    mir[1] = 1000000 * mir[1] / np.sum(mir[1])
    mrna = pd.read_table(mrna_expr_path, index_col=0)

    df = process_df(clash_path)
    df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min, 'au_before':np.mean,'au_after':np.mean }).sort_values('count').reset_index().set_index('mrna')
    df['mrna_exp'] = mrna['tpm']
    df = df[df['count'] >= min_counts]
    df = df.reset_index().set_index('mir')
    df['mir_exp'] = mir[1]
    df['hyb'] = df['count']/ (df['mir_exp'] * df['mrna_exp'])
    df = df[(df['mrna_exp'] > 1) & (df['mir_exp'] > 1)]
    return df



mhv68_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/smallfrac/CL100128477_L1_WHMOUkmcSAAETAASE-545_1.fq.counts.tsv']
mhv68_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/longfrac/HE2-1.genes.tsv']
mhv68_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*shuffled_dg.tsv')
mhv68_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*annotated')


kshv_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/smallfrac/541_WT-TIVE-1.small_fraction.counts.tsv']
kshv_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/longfrac/TIVE-1.genes.tsv']
kshv_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*shuffled_dg.tsv')
kshv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*annotated')


akata_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*shuffled_dg.tsv')
akata_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/smallfrac/*')
akata_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/akata/longfrac/Akata.genes.tsv']
akata_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*annotated')


snu_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/smallfrac/*')
snu_long_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/longfrac/*')
snu_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*shuffled_dg.tsv')
snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')


min_counts = 30
df = integrate(mhv68_small_frac_paths[0], mhv68_long_frac_paths[0], '/Users/nate/Projects/EBV_interactome/clash_scripts/SRR8395247_1_comp_mouse_mhv68.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated', min_counts)
df['species']=df.index.map(lambda x:'host'  if 'mmu' in x else 'mghv')
dfe=df[df['species']=='mghv']
dfh=df[df['species']=='host']

human_low = dfh[round(np.sum(dfh[['au_before','au_after']],1))<=2].hyb
human_high = dfh[round(np.sum(dfh[['au_before','au_after']],1))>2].hyb
ebv_low = dfe[round(np.sum(dfe[['au_before','au_after']],1))<=2].hyb
ebv_high = dfe[round(np.sum(dfe[['au_before','au_after']],1))>2].hyb
dfe['au_sum'] = np.sum(dfe[['au_before', 'au_after']],1)
dfh['au_sum'] = np.sum(dfh[['au_before', 'au_after']],1)
print(ks_2samp(dfe['au_sum'], dfh['au_sum']))
print(np.mean(dfe['au_sum']), np.mean(dfh['au_sum']))

human_low_sem = sem(dfh[round(np.sum(dfh[['au_before','au_after']],1))<=2].hyb)
human_high_sem = sem(dfh[round(np.sum(dfh[['au_before','au_after']],1))>2].hyb)
ebv_low_sem = sem(dfe[round(np.sum(dfe[['au_before','au_after']],1))<=2].hyb)
ebv_high_sem = sem(dfe[round(np.sum(dfe[['au_before','au_after']],1))>2].hyb)

xs= [0, 1, 3, 4]
means = [human_low, human_high, ebv_low, ebv_high]
sems = [[0,0,0,0],[human_low_sem, human_high_sem,ebv_low_sem, ebv_high_sem]]

fig = plt.figure(figsize=(4,8))
ax = plt.subplot()
#fig, ax = plt.subplots()
colors = ['#C1272D', '#C1272D', '#F15A24', '#F15A24']
for n, data in enumerate([human_low, human_high, ebv_low, ebv_high]):
    xs = n + np.random.rand(len(data.index)) / 3
    median = np.median(data)
    plt.scatter(xs, data, s=2, alpha=.6, lw=0,c=colors[n])
    plt.plot([n, n+ 1/3],[median, median], lw=2,color='k')


plt.yscale('log')
ax.set_yticklabels([])
plt.xticks([])
ax.set_xticklabels([])
plt.savefig('/Users/nate/Projects/EBV_interactome/mhv_au_usage.svg')


print(np.mean(human_low), np.mean(human_high))
print(ks_2samp(human_low, human_high) )
print(np.mean(ebv_low), np.mean(ebv_high))
print (ks_2samp(ebv_low, ebv_high))


min_counts = 30
df = integrate(snu_small_frac_paths[0], snu_long_frac_paths[0], '/Users/nate/Projects/EBV_interactome/clash_scripts/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_1_1_comp_human_ebv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated', min_counts)
df['species']=df.index.map(lambda x:'host'  if 'hsa' in x else 'ebv')
dfe=df[df['species']=='ebv']
dfh=df[df['species']=='host']

human_low = dfh[round(np.sum(dfh[['au_before','au_after']],1))<=2].hyb
human_high = dfh[round(np.sum(dfh[['au_before','au_after']],1))>2].hyb
ebv_low = dfe[round(np.sum(dfe[['au_before','au_after']],1))<=2].hyb
ebv_high = dfe[round(np.sum(dfe[['au_before','au_after']],1))>2].hyb
dfe['au_sum'] = np.sum(dfe[['au_before', 'au_after']],1)
dfh['au_sum'] = np.sum(dfh[['au_before', 'au_after']],1)
print(ks_2samp(dfe['au_sum'], dfh['au_sum']))
print(np.mean(dfe['au_sum']), np.mean(dfh['au_sum']))

human_low_sem = sem(dfh[round(np.sum(dfh[['au_before','au_after']],1))<=2].hyb)
human_high_sem = sem(dfh[round(np.sum(dfh[['au_before','au_after']],1))>2].hyb)
ebv_low_sem = sem(dfe[round(np.sum(dfe[['au_before','au_after']],1))<=2].hyb)
ebv_high_sem = sem(dfe[round(np.sum(dfe[['au_before','au_after']],1))>2].hyb)

xs= [0, 1, 3, 4]
means = [human_low, human_high, ebv_low, ebv_high]
sems = [[0,0,0,0],[human_low_sem, human_high_sem,ebv_low_sem, ebv_high_sem]]

fig = plt.figure(figsize=(4,8))
ax = plt.subplot()
#fig, ax = plt.subplots()

for n, data in enumerate([human_low, human_high, ebv_low, ebv_high]):
    xs = n + np.random.rand(len(data.index)) / 3
    median = np.median(data)
    plt.scatter(xs, data, s=2, alpha=.6, lw=0)
    plt.plot([n, n+ 1/3],[median, median], lw=2,color='k')


plt.yscale('log')
plt.savefig('/Users/nate/Projects/EBV_interactome/au_usage.svg')


print(np.mean(human_low), np.mean(human_high))
print(ks_2samp(human_low, human_high) )
print(np.mean(ebv_low), np.mean(ebv_high))
print (ks_2samp(ebv_low, ebv_high))






min_counts = 30
df = integrate(kshv_small_frac_paths[0], kshv_long_frac_paths[0], '/Users/nate/Projects/EBV_interactome/clash_scripts/SRR5876952_1_comp_human_kshv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated', min_counts)
df['species']=df.index.map(lambda x:'host'  if 'hsa' in x else 'kshv')
dfe=df[df['species']=='kshv']
dfh=df[df['species']=='host']

human_low = dfh[round(np.sum(dfh[['au_before','au_after']],1))<=2].hyb
human_high = dfh[round(np.sum(dfh[['au_before','au_after']],1))>2].hyb
ebv_low = dfe[round(np.sum(dfe[['au_before','au_after']],1))<=2].hyb
ebv_high = dfe[round(np.sum(dfe[['au_before','au_after']],1))>2].hyb
dfe['au_sum'] = np.sum(dfe[['au_before', 'au_after']],1)
dfh['au_sum'] = np.sum(dfh[['au_before', 'au_after']],1)
print(ks_2samp(dfe['au_sum'], dfh['au_sum']))
print(np.mean(dfe['au_sum']), np.mean(dfh['au_sum']))

human_low_sem = sem(dfh[round(np.sum(dfh[['au_before','au_after']],1))<=2].hyb)
human_high_sem = sem(dfh[round(np.sum(dfh[['au_before','au_after']],1))>2].hyb)
ebv_low_sem = sem(dfe[round(np.sum(dfe[['au_before','au_after']],1))<=2].hyb)
ebv_high_sem = sem(dfe[round(np.sum(dfe[['au_before','au_after']],1))>2].hyb)

xs= [0, 1, 3, 4]
means = [human_low, human_high, ebv_low, ebv_high]
sems = [[0,0,0,0],[human_low_sem, human_high_sem,ebv_low_sem, ebv_high_sem]]

fig = plt.figure(figsize=(4,8))
ax = plt.subplot()
#fig, ax = plt.subplots()
colors = ['#0071BC', '#0071BC', '#662D91', '#662D91']
for n, data in enumerate([human_low, human_high, ebv_low, ebv_high]):
    xs = n + np.random.rand(len(data.index)) / 3
    median = np.median(data)
    plt.scatter(xs, data, s=2, alpha=.6, lw=0,c=colors[n])
    plt.plot([n, n+ 1/3],[median, median], lw=2,color='k')


plt.yscale('log')
ax.set_yticklabels([])
plt.xticks([])
ax.set_xticklabels([])
plt.savefig('/Users/nate/Projects/EBV_interactome/kshv_au_usage.svg')

print(np.mean(human_low), np.mean(human_high))
print(ks_2samp(human_low, human_high) )
print(np.mean(ebv_low), np.mean(ebv_high))
print (ks_2samp(ebv_low, ebv_high))

