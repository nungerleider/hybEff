## Figure 2 ##

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu as mnu
from scipy.stats import ttest_ind, pearsonr, spearmanr, ks_2samp, ttest_rel, sem
import glob


def process_mir_df(path):

    df = pd.read_table(path, index_col=0, header=None)
    if '*' in df.index:
        df = df.drop('*')  # Unmapped counts
    df = df.loc[set(i for i in df.index if 'miR' in i or 'let' in i)]
    cpm = 1000000 * df / np.sum(df)  # counts per million mapped mirs
    
    return cpm


def process_ago_bound_mirs(path):

    df = pd.read_table(path, index_col=0)
    cpm = 1000000 * df / np.sum(df)
    return cpm


def process_clash(path):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt' in i]) | set(['SLMAP','HASPIN','VMP1', 'MTRNR2L12', 'Vmp1'])  # mRNAs that have lots of small fraction coverage (not real mRNA) or mtRNA
    counts = pd.DataFrame(df[~df['mrna'].isin(blacklist)].groupby('mir')['count'].sum().sort_values())
    cpm = 1000000 * counts / np.sum(counts)
    return cpm


def vir_over_total(df, virus='ebv'):
    return float(100 * np.sum(df.loc[[i for i in df.index if virus in i]]) / np.sum(df))
     

def viral_small_fraction(paths, process_func, virus='ebv'):
    arr = []
    for path in paths:
        df = process_func(path)
        percent = vir_over_total(df, virus=virus)
        arr.append(percent)
    return np.array(arr)


def plot_viral_fraction(small_frac_paths, clash_paths, virus='ebv', output_path='ago_enrichment_summary.svg'):

    vir_frac_expressed = viral_small_fraction(paths=small_frac_paths, process_func=process_mir_df,virus=virus)
    vir_frac_mrnabound = viral_small_fraction(paths=clash_paths, process_func=process_clash,virus=virus)

    fig = plt.figure(figsize=(2,6))
    ax = plt.subplot()
    expressed_x = [.1 * i for i in range(len(vir_frac_expressed))]
    mrnabound_x = [.4, .5, .6]
    plt.scatter(expressed_x, vir_frac_expressed, c='#FF9F1C', s=150, ec='k', lw=1, marker='D')
    plt.scatter(mrnabound_x, vir_frac_mrnabound, c='#1A936F', s=150, ec='k', lw=1, marker='D')

    plt.yscale('log')
    plt.ylim([.1,100])
    plt.xlim([-.2, .8])
    ax.set_xticks([])
    ax.set_yticklabels(['' for i in ax.get_yticklabels()])
    plt.savefig(output_path, dpi=500)

    ## T-test ##
    print(ttest_ind(vir_frac_expressed, vir_frac_mrnabound))



def get_mir_means(small_frac_paths):
    '''If more than one small frac seq sample, get average c.p.m. across samples'''

    mirs = set()
    for path in small_frac_paths:
        small = process_mir_df(path)
        mirs = mirs | set(small.index)

    mirdf = pd.DataFrame(index=mirs)
    for path in small_frac_paths:
        small = process_mir_df(path)
        mirdf[path] = small[1]

    return np.mean(mirdf, 1)


def get_clash_means(clash_paths):
    '''Get average h.c.p.m. across samples'''

    clashes = set()
    for path in clash_paths:
        df = process_clash(path)
        clash = pd.DataFrame(df.groupby(['mir'])['count'].sum())
        clashes = clashes | set(clash.index)

    clashdf = pd.DataFrame(index=clashes)
    for path in clash_paths:
        df = process_clash(path)
        clash = pd.DataFrame(df.groupby(['mir'])['count'].sum())
        clashdf[path] = clash['count']

    return np.mean(clashdf.fillna(0),1)


def plot_individual_mir_ago_enrich(small_frac_paths, clash_paths, min_mir_cpm, virus, output_path):
    
    mir_means = get_mir_means(small_frac_paths)
    clash_means = get_clash_means(clash_paths)

    mir_means = mir_means[mir_means > min_mir_cpm]
    inds = list(set(clash_means.index) & set(mir_means.index))
    mir_means = mir_means[inds]
    clash_means = clash_means[inds]
    fc = pd.DataFrame(clash_means / mir_means)
    fc = fc.sort_values(0)

    fc['rank'] = range(len(fc.index))
    host = fc.loc[[i for i in fc.index if virus not in i]]
    viral = fc.loc[[i for i in fc.index if virus in i]]
    host_bubble_size = np.sqrt(mir_means.loc[host.index])
    vir_bubble_size = np.sqrt(mir_means.loc[viral.index])
    fig = plt.figure(figsize=(8,4))
    ax = plt.subplot()
    plt.scatter(host['rank'], host[0], c='#A4ADBE', s=host_bubble_size, alpha=.4, lw=1)
    plt.scatter(viral['rank'], viral[0], c='#E03654', s=vir_bubble_size, alpha=.4, lw=1)
    plt.yscale('log')
    plt.ylim([.00001, 500])
    ax.set_yticklabels([])
    plt.savefig(output_path, dpi=500)

    ## Mann Whitney U test ##
    print(mnu(viral[0], host[0], alternative='greater'))




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





# Figure 2A
# T-test (independent) for triplicate small frac seq cell lines (Akata and SNU719) - otherwise, no stats done.

plot_viral_fraction(small_frac_paths=akata_small_frac_paths, clash_paths=akata_mrna_clash_paths, virus='ebv', output_path='akata_summary.svg')
plot_viral_fraction(small_frac_paths=mhv68_small_frac_paths, clash_paths=mhv68_mrna_clash_paths, virus='mghv', output_path='mhv68_total_v_clash_summary.svg')
plot_viral_fraction(small_frac_paths=kshv_small_frac_paths, clash_paths=kshv_mrna_clash_paths, virus='kshv', output_path='kshv_total_v_clash_summary.svg')
plot_viral_fraction(small_frac_paths=snu_small_frac_paths, clash_paths=snu_mrna_clash_paths, virus='ebv', output_path='snu719_total_v_clash_summary.svg')


#################################################


# Figure 2B
# Mann Whitney U

plot_individual_mir_ago_enrich(mhv68_small_frac_paths, mhv68_mrna_clash_paths, min_mir_cpm=30, virus='mghv', output_path='mhv68_rank_plot_ago_enrichment.svg')
plot_individual_mir_ago_enrich(kshv_small_frac_paths, kshv_mrna_clash_paths, min_mir_cpm=30, virus='kshv', output_path='kshv_rank_plot_ago_enrichment.svg')
plot_individual_mir_ago_enrich(snu_small_frac_paths, snu_mrna_clash_paths, min_mir_cpm=30, virus='ebv', output_path='snu719_rank_plot_ago_enrichment.svg')
plot_individual_mir_ago_enrich(akata_small_frac_paths, akata_mrna_clash_paths, min_mir_cpm=30, virus='ebv', output_path='akata_rank_plot_ago_enrichment.svg')

#################################################
