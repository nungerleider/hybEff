# Figure 6A

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
    df = df.groupby(['mir', 'mrna', 'mrna_feature_start'])['count'].sum().sort_values().reset_index()
    df.index = [f'{i}_{j}_{k}' for i,j,k in zip(df['mir'], df['mrna'], df['mrna_feature_start'])]
    return df

def get_vir_clash_means(paths, min_hyb_counts=20, vir='ebv'):

    l = set()
    for i in snu_mrna_clash_paths:
        x = process_df(i)
        l = l | set(x.index)

    new = pd.DataFrame(index=l)
    for i in snu_mrna_clash_paths:
        x = process_df(i)
        new[i] = x['count']

    new = np.mean(new, 1)
    new = new.sort_values()
    new = new.loc[set(i for i in new.index if vir in i)].sort_values()
    new = new[new >= min_hyb_counts]
    return new

##SNU719
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


def bubble_plot(clash_paths, prefix='SNU719', min_hyb_counts=20, vir='ebv'):
    
    clash_means = get_vir_clash_means(clash_paths, min_hyb_counts=min_hyb_counts, vir=vir)
    mirs = list(set([i.split('_')[0] for i in clash_means.index]))
    mirs.sort(key=lambda x:int(x.split('T')[1].split('-')[0]))
    position = 0
    mirs2 = []
    fig, ax = plt.subplots(figsize=(12,8))
    for mir in mirs:
        filt = new.loc[set(i for i in new.index if mir in i)]
        if np.sum(filt) < 300:
            continue 
        filtmir = pd.DataFrame(filt.sort_values())
        filtmir['frac'] = filtmir[0]/np.sum(filtmir[0])
        color = []
        for i in filtmir.index:
            if "3'UTR" in i:
                color.append('#ed1c24')
            elif "5'UTR" in i:
                color.append('#235789')
            elif 'CDS' in i:
                color.append('#f1d302')
            else:
                color.append('0.75')
        plt.scatter([position] * len(filtmir.index), filtmir['frac'], c=color, s=filtmir[0]/15, alpha=.7, ec='k')
        position += 1
        mirs2.append(mir)

    plt.savefig(f'{prefix}_mir_clash_bubble.svg')
    with open(f'{prefix}_mir_clash_bubble_xaxis.tsv', 'w') as outfile:
        for i in mirs2:
            outfile.write(i + '\n')
    

## Figure 6A
bubble_plot(snu_mrna_clash_paths, min_hyb_counts=20)