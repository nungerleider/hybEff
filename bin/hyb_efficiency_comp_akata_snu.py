## Figure 3A, 3B 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr, spearmanr, ks_2samp, ttest_rel, sem
import glob
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles


def process_clash(paths, vir='ebv'):

    ind = set()
    for i in paths:

        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i ]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j]) for i,j in zip(df['mir'], df['mrna'])]                                                                              

        ind = ind | set(df.index)
    new = pd.DataFrame(index=ind)

    for i in paths:
        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i ]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j]) for i,j in zip(df['mir'], df['mrna'])]                                                                              
        new[i] = df['count']
    new = np.mean(new,1)
    return new


def process_mir_df(path):

    df = pd.read_table(path, index_col=0, header=None)
    if '*' in df.index:
        df = df.drop('*')  # Unmapped counts
    df = df.loc[set(i for i in df.index if 'miR' in i or 'let' in i)]
    cpm = 1000000 * df / np.sum(df)  # counts per million mapped mirs
    
    return cpm


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



def get_mrna_means(mrna_expr_paths):
    '''If more than one sample, get average t.p.m. across samples'''

    mrnas = set()
    for path in mrna_expr_paths:
        mrna = pd.read_table(path, index_col=0)
        mrnas = mrnas | set(mrna.index)

    mrna_df = pd.DataFrame(index=mrnas)
    for path in mrna_expr_paths:
        mrna = pd.read_table(path, index_col=0)
        mrna_df[path] = mrna['tpm']

    return np.mean(mrna_df, 1)



def integrate_df(mir_df, mrna_df, clash_df, min_hyb_counts=100, min_mir_counts=20, min_mrna_tpm=1):
    
    clash_df = clash_df.reset_index()
    clash_df.index = clash_df['index'].map(lambda x: x.split('_')[1])
    clash_df['mrna_exp'] = mrna_df
    clash_df = clash_df[clash_df[0] >= min_hyb_counts]

    clash_df.index = clash_df['index'].map(lambda x: x.split('_')[0])
    clash_df['mir_exp'] = mir_df
    clash_df = clash_df.set_index('index')
    clash_df['hyb'] = clash_df[0] / (clash_df['mir_exp'] * clash_df['mrna_exp'])
    clash_df = clash_df[(clash_df['mrna_exp'] > min_mrna_tpm) & (clash_df['mir_exp'] > min_mir_counts)]
    return clash_df


def correlate(ak_df, snu_df, measure='hyb'):

    rows = set(snu_df.index) & set(ak_df.index)
    new = pd.DataFrame(index=rows) 
    new['snu'] = snu_df[measure]
    new['ak'] = ak_df[measure]
    ebv = new.loc[[i for i in new.index if 'ebv' in i]]
    hsa = new.loc[[i for i in new.index if 'hsa' in i]]
    fig = plt.figure()
    ax = plt.subplot()
    plt.scatter(np.log2(hsa['ak']), np.log2(hsa['snu']), color='.7', alpha=.4, s=50)
    plt.scatter(np.log2(ebv['ak']), np.log2(ebv['snu']), color='red', alpha=.4, s=50)

    print('both', spearmanr(np.log2(new['ak']), np.log2(new['snu'])))
    #print('ebv', spearmanr(np.log2(ebv['ak']), np.log2(ebv['snu'])))
    #print('hsa', spearmanr(np.log2(hsa['ak']), np.log2(hsa['snu'])))
    # print('both', pearsonr(np.log2(new['ak']), np.log2(new['snu'])))
    # print('ebv', pearsonr(np.log2(ebv['ak']), np.log2(ebv['snu'])))
    # print('hsa', pearsonr(np.log2(hsa['ak']), np.log2(hsa['snu'])))
    ax.grid(ls='--', alpha=.4)
    ax.set_axisbelow(True)

    plt.savefig(f'{measure}_hyb_efficiency_scatter.svg')


def common_hybrids(akata_df, snu_df, min_hyb_counts=20, min_mir_counts=20, min_mrna_tpm=1):
    
    human_akata = akata_df.loc[[i for i in akata_df.index if 'hsa' in i]]
    human_snu = snu_df.loc[[i for i in snu_df.index if 'hsa' in i]]
    ebv_akata = akata_df.loc[[i for i in akata_df.index if 'ebv' in i]]
    ebv_snu = snu_df.loc[[i for i in snu_df.index if 'ebv' in i]]

    human_akata_only = set(human_akata.index) - set(human_snu.index)
    human_snu_only = set(human_snu.index) - set(human_akata.index)
    human_both = set(human_akata.index) & set(human_snu.index)

    ebv_akata_only = set(ebv_akata.index) - set(ebv_snu.index)
    ebv_snu_only = set(ebv_snu.index) - set(ebv_akata.index)
    ebv_both = set(ebv_akata.index) & set(ebv_snu.index)
    #print(human)
    print('human_akata_only', len(human_akata_only))
    print('human_snu_only', len(human_snu_only))
    print('human_both', len(human_both))
    print('ebv_akata_only', len(ebv_akata_only))
    print('ebv_snu_only', len(ebv_snu_only))
    print('ebv both', len(ebv_both))


akata_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*shuffled_dg.tsv')
akata_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/smallfrac/*')
akata_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/akata/longfrac/Akata.genes.tsv']
akata_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*annotated')


snu_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/smallfrac/*')
snu_long_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/longfrac/*')
snu_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*shuffled_dg.tsv')
snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')



## Figure 3A, 3B ##
def main(min_hyb_counts=100, min_mir_counts=20, min_mrna_tpm=5):
    akmir = get_mir_means(akata_small_frac_paths)
    akmrna = get_mrna_means(akata_long_frac_paths)
    akclash = process_clash(akata_mrna_clash_paths)
    ak_df = integrate_df(akmir, akmrna, akclash, min_hyb_counts, min_mir_counts, min_mrna_tpm)

    snumir = get_mir_means(snu_small_frac_paths)
    snumrna = get_mrna_means(snu_long_frac_paths)
    snuclash = process_clash(snu_mrna_clash_paths)
    snu_df = integrate_df(snumir, snumrna, snuclash, min_hyb_counts, min_mir_counts, min_mrna_tpm)
    print('N: ', len(ak_df), len(snu_df))
    correlate(ak_df, snu_df, measure='hyb')
    correlate(ak_df, snu_df, measure=0) # "0" column is h.c.p.m., 'hyb' columns is targeting efficiency

    common_hybrids(ak_df, snu_df, min_hyb_counts, min_mir_counts, min_mrna_tpm)

main()
