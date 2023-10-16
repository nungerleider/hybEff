# Figure 4D

import subprocess
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp, sem, spearmanr
import matplotlib.pyplot as plt
import glob
from collections import defaultdict


def m13p(mirfold, ind): 
    y = 0 
    for i in mirfold[ind-1:]: 
        if i=='(': 
            y+=1 
    return y 


def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    df['m13p'] = df['mir_fold'].map(lambda x:m13p(x,13))

    return df


def get_seed(df):

    mapping = {
        0: "No seed",
        1.5: "8mer mismatch",
        2: "6mer",
        4: "7merM8",
        6: "7merA1",
        8: "8mer",
        1: "No seed, supplemental",
        2.5: "8mer mismatch, supplemental",
        3: "6mer, supplemental",
        5: "7merM8, supplemental",
        7: "7merA1, supplemental",
        9: "8mer, supplemental"
    }

    df = df[(df['a1'] != 'unknown') & (df['a1'] != 'unannotated')]
    a1 = [int(i) for i in df['a1']]
    core = df['m2_5'] + df['m6_7']
    m8 = df['m8']
    supplemental = [True if i > 3 else False for i in df['m13_17']]
    seeds = []
    for i in range(len(df.index)):
        num = 0
        if supplemental[i] is True:
            num += 1      
        if core[i]  == 6:
            num += 2      
            if a1[i] == 1:
                num += 4         
            if m8[i] == 1:
                num += 2      
        elif core[i] < 6 and a1[i] == 1 and m8[i] == 1:
            num += 1.5
        seeds.append(num)
    df["seed"] = seeds
    return df

def map_seed(df):

    mapping = {
        0: "No seed",
        1.5: "8mer mismatch",
        2: "6mer",
        4: "7merM8",
        6: "7merA1",
        8: "8mer",
        1: "No seed, supplemental",
        2.5: "8mer mismatch, supplemental",
        3: "6mer, supplemental",
        5: "7merM8, supplemental",
        7: "7merA1, supplemental",
        9: "8mer, supplemental",
    }

    df['seed'] = df.seed.map(lambda x: mapping(x))
    return df

def load_fa(path):
    fa_dict = defaultdict(str)
    with open(path) as fa_in:
        for line in fa_in:
            line = line.strip('\n')
            if line[0] == '>':
                seq_id = line.strip('>')
            else:
                fa_dict[seq_id.split('.')[0]]+=line
    return fa_dict


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


def integrate(small_frac_paths, mrna_expr_paths, clash_path, min_hyb_counts=1000, min_mir_counts=20, min_mrna_tpm=1, virus='ebv'):
    
    mir = get_mir_means(small_frac_paths)
    mrna = get_mrna_means(mrna_expr_paths)

    df = process_df(clash_path, vir=vir)
 
    df = df.groupby(['mir','transcript']).agg({'mrna':lambda x:x[0],'count':sum, 'binding_energy':min, 'seed': max, 'species': lambda x:x[0], 'match_start': lambda x:x[0], 'mrna_feature_start':lambda x:x[0],'au_before':np.mean,'au_after':np.mean}).sort_values('count').reset_index().set_index('mrna')
    df['mrna_exp'] = mrna
    df = df.groupby(['mir','mrna']).agg({'transcript':lambda x:x[0], 'mrna_exp': sum, 'count':sum, 'binding_energy':min, 'seed': max, 'species': lambda x:x[0], 'match_start': lambda x:x[0], 'mrna_feature_start':lambda x:x[0],'au_before':np.mean,'au_after':np.mean}).sort_values('count').reset_index().set_index('mrna')

    df = df.reset_index().set_index('mir')
    df['mir_exp'] = mir
    df = df[(df['count'] >= min_hyb_counts) & (df['mir_exp'] >= min_mir_counts) & (df['mrna_exp'] >= min_mrna_tpm)] # Filter out low abundance hybs, mirs, and mrnas
    df['hyb'] = df['count']/ (df['mir_exp'] * df['mrna_exp'])
    return df


def process_clash(paths, vir='ebv'):

    ind = set()
    for i in paths:

        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j]) for i,j in zip(df['mir'], df['mrna'])]                                                                              

        ind = ind | set(df.index)
    new = pd.DataFrame(index=ind)

    for i in paths:
        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j]) for i,j in zip(df['mir'], df['mrna'])]                                                                              
        new[i] = df['count']
    new = np.mean(new,1)
    return new


def get_stabilities(clash_df, fa_path, min_hyb_counts=1000, min_mir_counts=20, min_mrna_tpm=1):
    
    fa_dict = load_fa(fa_path)
    clash_df = clash_df[(clash_df['count']>=min_hyb_counts) & (clash_df['mir_exp']>= min_mir_counts) & (clash_df['mrna_exp'] >=min_mrna_tpm)]
    stabilities = []
    n = 0
    for start, transcript in zip(clash_df['match_start'], clash_df['transcript']):
        try:
            seq = fa_dict[transcript]
        except KeyError:
            seq = ''
        folding = subprocess.check_output(f'echo {seq} | RNAplfold -u 14 -L 40 -W 80', shell=True).decode('utf8')
        x = pd.read_table('plfold_lunp', skiprows=2, header=None, index_col=0)
        start -= 7
        try:
            stability = x.loc[start, 14]
        except:
            stability = 'na'
        stabilities.append(stability)
        print(len(clash_df.index) - n)
        n+=1

    clash_df['stabilities'] = stabilities
    clash_df = clash_df[clash_df['stabilities'].map(lambda x:True if type(x) is not str else False)]
    clash_df['stabilities'] = clash_df['stabilities'].map(lambda x:float(x))
    clash_df = clash_df.dropna()
    return clash_df


def plot_stabilities(df, seed_min=7, vir='ebv', utr=True):

    if utr is True:
        df = df[df['mrna_feature_start'] == "3'UTR"]
    df = df[df['seed'] >= seed_min]
    virus = df[df['species']==vir]
    host = df[df['species']=='host']
    fig = plt.figure(figsize=(4,8))
    ax = plt.subplot()
    
    for n, data in enumerate([host, virus]):
        xs = n + np.random.rand(len(data.index)) / 2
        median = np.median(data['stabilities'])
        plt.scatter(xs, data['stabilities'], s=30, alpha=.6, lw=0)
        plt.plot([n, n+ 1/2],[median, median], lw=2,color='k')
    plt.yscale('log')

    print(ks_2samp(virus['stabilities'], host['stabilities']))

    plt.savefig(f'/Users/nate/Projects/EBV_interactome/{vir}_site_accessibility.svg')



def plot_au(df, vir):

    virus = df[df['species'] == vir]
    host = df[df['species'] != vir]

    virus_mean_au = np.mean(virus['au_sum'])
    host_mean_au = np.mean(host['au_sum'])
    host_sem = sem(host['au_sum'])
    virus_sem = sem(virus['au_sum'])

    print(ks_2samp(virus['au_sum'], host['au_sum']))
    print('virus:', virus_mean_au, 'host:', host_mean_au)

    host_low = host[host['au_sum'] <= 2].hyb
    host_high = host[host['au_sum'] > 2].hyb
    virus_low = virus[virus['au_sum'] <= 2].hyb
    virus_high = virus[virus['au_sum'] > 2].hyb

    host_low_sem = sem(host_low)
    host_high_sem = sem(host_high)
    virus_low_sem = sem(virus_low)
    virus_high_sem = sem(virus_high)
    
    print(np.mean(host_low), np.mean(host_high))
    print(ks_2samp(host_low, host_high) )
    print(np.mean(virus_low), np.mean(virus_high))
    print (ks_2samp(virus_low, virus_high))

    # Column-scatter plot
    # Hyb compared w/ binned AU content
    xs = [0, 1, 3, 4]
    sems = [[0,0,0,0], [host_low_sem, host_high_sem,virus_low_sem, virus_high_sem]]
    fig = plt.figure(figsize=(4,8))
    ax = plt.subplot()
    for n, data in enumerate([host_low, host_high, virus_low, virus_high]):
        xs = n + np.random.rand(len(data.index)) / 3
        median = np.median(data)
        plt.scatter(xs, data, s=2, alpha=.6, lw=0)
        plt.plot([n, n+ 1/3],[median, median], lw=2,color='k')
    plt.yscale('log')
    plt.savefig(f'/Users/nate/Projects/EBV_interactome/{vir}_au_usage_hyb_scatter.svg')
    
    fig2 = plt.figure(figsize=(4,8))
    ax2 = plt.subplot()
    plt.bar([0, 1], [host_mean_au, virus_mean_au], yerr=[[0,0],[host_sem, virus_sem]])
    plt.savefig(f'/Users/nate/Projects/EBV_interactome/{vir}_au_usage_mean_bar.svg')



def main(clash_paths, long_frac_paths, small_frac_paths, fa_path, vir='ebv', seed_min=7, utr=True, min_hyb_counts=100, min_mir_counts=20, min_mrna_tpm=5):
    df = process_df(clash_paths[0])
    df = get_seed(df)
    df = df.groupby(['mir','transcript']).agg({'mrna':lambda x:x[0],'count':sum, 'binding_energy':min, 'seed': max, 'species': lambda x:x[0], 'match_start': lambda x:x[0], 'mrna_feature_start':lambda x:x[0],'au_before':np.mean,'au_after':np.mean}).sort_values('count').reset_index().set_index('mrna')
    
    mrna = get_mrna_means(long_frac_paths)  
    df['mrna_exp'] = mrna
    df = df.reset_index().set_index('mir')
    mir = get_mir_means(small_frac_paths)
    df['mir_exp'] = mir
    df = df.reset_index()
    cl = process_clash(clash_paths)
    df.index = [f'{i}_{j}' for i,j in zip(df['mir'], df['mrna'])] 
    df['count'] = cl
    df['hyb'] = df['count']/ (df['mir_exp'] * df['mrna_exp'])
    df['species'] = df.index.map(lambda x: vir if vir in x else 'host')
    df['au_sum'] = np.sum(df[['au_before', 'au_after']], 1)
    df = get_stabilities(df, fa_path=fa_path, min_hyb_counts=min_hyb_counts, min_mir_counts=min_mir_counts, min_mrna_tpm=min_mrna_tpm)
    df.to_csv(f'{vir}_site_stability.tsv', sep='\t')
    plot_stabilities(df, seed_min=seed_min, vir=vir, utr=utr)
    plot_au(df, vir)






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


# KS test 
main(snu_mrna_clash_paths, snu_long_frac_paths, snu_small_frac_paths, fa_path='/Users/nate/Documents/Genomes/Fastas/gencode.v34.transcripts.fa', vir='ebv', seed_min=8, utr=True, min_hyb_counts=100, min_mir_counts=20, min_mrna_tpm=5)
main(kshv_mrna_clash_paths, kshv_long_frac_paths, kshv_small_frac_paths, fa_path='/Users/nate/Documents/Genomes/Fastas/gencode.v34.transcripts.fa', vir='kshv', seed_min=8, utr=True, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=5)
main(mhv68_mrna_clash_paths, mhv68_long_frac_paths, mhv68_small_frac_paths, fa_path='/Users/nate/Documents/Genomes/Fastas/gencode.vM25.transcripts.fa', vir='mghv', seed_min=8, utr=True, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=5)



