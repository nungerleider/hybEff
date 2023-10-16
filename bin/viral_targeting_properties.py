import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr, spearmanr, ks_2samp, ttest_rel, sem
import glob
import seaborn as sns
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles


def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    df = df.groupby(['mir', 'mrna', 'mrna_feature_start'])['count'].sum().sort_values().reset_index()
    df.index = [f'{i}_{j}_{k}' for i,j,k in zip(df['mir'], df['mrna'], df['mrna_feature_start'])]
    return df

akata_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*annotated')
akata_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/smallfrac/*')
akata_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/akata/longfrac/Akata.genes.tsv']

snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')
snu_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/smallfrac/*')
snu_long_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/longfrac/*')

mghv_mrna_expression_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/longfrac/HE2-1.genes.tsv']
mghv_mir_expression_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/smallfrac/CL100128477_L1_WHMOUkmcSAAETAASE-545_1.fq.counts.tsv']
mghv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/SRR839524[5-7]*annotated')# + glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/SRR8395250*')  
mghv_mrna_clash_paths.reverse()

kshv_mir_expression_paths = ['/Users/nate/Projects/EBV_interactome/kshv/smallfrac/541_WT-TIVE-1.small_fraction.counts.tsv']
kshv_mrna_expression_paths = ['/Users/nate/Projects/EBV_interactome/kshv/longfrac/TIVE-1.genes.tsv']
kshv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*annotated')

def distplots(clash_paths=snu_mrna_clash_paths, min_counts=300, vir='ebv',cell_line='snu719'):

    l = set()
    for i in clash_paths:
        x = process_df(i)
        l = l | set(x.index)


    new = pd.DataFrame(index=l)
    for i in clash_paths:
        x = process_df(i)
        new[i] = x['count']

    new = np.mean(new,1)
    new = new.sort_values()
    new = new[new > 20]
    new = pd.DataFrame(new)

    new['mir'] = new.index.map(lambda x:x.split('_')[0])
    new['mrna'] = new.index.map(lambda x:x.split('_')[1])
    new['site'] = new.index.map(lambda x:x.split('_')[-1])
    new.species = new['species'] = new.index.map(lambda x:vir if vir in x else 'host')
    new.columns = ['counts', 'mir', 'mrna', 'site', 'species']
    new2 = new.groupby(['mir'])['counts'].sum()
    new = new.reset_index().set_index('mir')
    new['mir_sum'] = new2
    new = new.reset_index().set_index('index')
    new['frac'] = new['counts'] / new['mir_sum']
    new2 = new.groupby(['mrna'])['counts'].sum()
    new = new.reset_index().set_index('mrna')
    new['mrna_sum'] = new2
    new['mrna_frac'] = new['counts'] / new['mrna_sum']
    new = new.reset_index().set_index('index')
    virus = new[new['species']==vir].sort_values('counts').iloc[-100:]
    host = new[new['species']=='host'].sort_values('counts').iloc[-100:]

    nn2  = new[(new['counts']>min_counts) ][['mrna_frac','frac', 'species']]
    
    print('mrna_frac', ks_2samp(virus['mrna_frac'], host['mrna_frac']))
    print('mir_frac', ks_2samp(np.log2(virus['frac']), np.log2(host['frac'])))


    fig, ax = plt.subplots()
    sns.distplot(host['mrna_frac'], kde_kws={"shade": True}, hist=False, rug=True)
    sns.distplot(virus['mrna_frac'], kde_kws={"shade": True}, hist=False, rug=True)
    plt.savefig(f'mrna_exclusivity_distplot_{cell_line}.svg')

    fig, ax = plt.subplots()
    sns.distplot(np.log2(host['frac']), kde_kws={"shade": True}, hist=False, rug=True)
    sns.distplot(np.log2(virus['frac']), kde_kws={"shade": True}, hist=False, rug=True)
    plt.savefig(f'mir_exclusivity_distplot_{cell_line}.svg')
    plt.close('all')
    return virus, host

v, h = distplots(snu_mrna_clash_paths,100,'ebv','snu719')
fig, ax = plt.subplots()
venn2([set(v['mrna']), set(h['mrna'])])
plt.savefig('snu_human_v_virus_venn.svg')
plt.close('all')

v,h = distplots(kshv_mrna_clash_paths,150,'kshv','tive')
fig, ax = plt.subplots()
venn2([set(v['mrna']), set(h['mrna'])])
plt.savefig('kshv_human_v_virus_venn.svg')
plt.close('all')

v,h = distplots(mghv_mrna_clash_paths,150,'mghv','be2')
fig, ax = plt.subplots()
venn2([set(v['mrna']), set(h['mrna'])])
plt.savefig('mghv_human_v_virus_venn.svg')
plt.close('all')