## Figure 3E ##

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
    

## For Violin plot
def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def violin(data, ax, position, color):

    parts = ax.violinplot(data, [position], showmeans=False, showmedians=False, showextrema=False)   
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75])
    whiskers = np.array(adjacent_values(data, quartile1, quartile3))
    whiskersMin, whiskersMax = whiskers[0], whiskers[1]

    for pc in parts['bodies']:
        pc.set_facecolor(color)  
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    ax.scatter(position, medians, marker='_', color='white', s=10, zorder=3)
    ax.vlines(position, quartile1, quartile3, color='k', linestyle='-', lw=4)
    ax.vlines(position, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)


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


def integrate(small_frac_paths, mrna_expr_paths, clash_path, min_hyb_counts=1000, min_mir_counts=20, min_mrna_tpm=1):
    
    mir = get_mir_means(small_frac_paths)
    mrna = get_mrna_means(mrna_expr_paths)

    df = process_df(clash_path)
#   df = df.groupby(['mir','mrna']).agg({'count':sum,'binding_energy':min, 'au_before':np.mean,'au_after':np.mean }).sort_values('count').reset_index().set_index('mrna')
    df = df.groupby(['mir','mrna']).agg({'count':sum,'binding_energy':min }).sort_values('count').reset_index().set_index('mrna') 
    df['mrna_exp'] = mrna
    df = df.reset_index().set_index('mir')
    df['mir_exp'] = mir
    df = df[(df['count'] >= min_hyb_counts) & (df['mir_exp'] >= min_mir_counts) & (df['mrna_exp'] >= min_mrna_tpm)] # Filter out low abundance hybs, mirs, and mrnas
    df['hyb'] = df['count']/ (df['mir_exp'] * df['mrna_exp'])
    return df


def plot_hyb_eff_violin(small_frac_paths, long_frac_paths, clash_paths, min_hyb_counts, min_mir_counts, min_mrna_tpm, virus, output_prefix):
    
    
    fig = plt.figure(figsize=(6.4, 4.8), dpi=300)
    ax = plt.subplot(111)
    colors = ['#E08DAC', '#57E2E5']
    vir_medians, host_medians = [], []
    for i in range(0, len(clash_paths)):
            
        df = integrate(small_frac_paths, long_frac_paths, clash_paths[i], min_hyb_counts, min_mir_counts, min_mrna_tpm)
        vir_hyb = df.loc[set([i for i in df.index if virus in i])]['hyb']
        host_hyb = df.loc[set([i for i in df.index if virus not in i])]['hyb']
        
        log_vir = np.log10(vir_hyb.sort_values())
        log_host = np.log10(host_hyb.sort_values())
        
        host_x = 3 * i
        vir_x = (3 * i) + 1

        host_color = colors[0]
        vir_color = colors[1]

        violin(log_host, ax, host_x, host_color)
        violin(log_vir, ax, vir_x, vir_color)

        vir_medians.append(np.median(vir_hyb))
        host_medians.append(np.median(host_hyb))

        print(ks_2samp(vir_hyb, host_hyb))
        print('medians:\n', virus, np.median(vir_hyb), 'host', np.median(host_hyb))
        print('virus N:', len(log_vir), 'host N:', len(log_host))

    print('median paired ttests (N=3):\n', ttest_rel(vir_medians, host_medians))

    ax.grid(ls='--', alpha=.2)
    ax.set_axisbelow(True)
    plt.ylim([-5.5,1])
    plt.savefig(f'{output_prefix}_hyb_efficiency_violin.svg')
  
    host_sem = sem(host_medians)
    vir_sem = sem(vir_medians)

    fig2 = plt.figure()
    ax2 = plt.subplot(111)
    plt.bar([0,1],[np.mean(host_medians), np.mean(vir_medians)], yerr=[[0,0],[host_sem, vir_sem]],color=[host_color, vir_color], ec='k', lw=2)
    ax2.set_xticks([])
    ax2.grid(alpha=.2, ls='--')
    ax2.set_axisbelow(True)
    ax2.set_yticks([ax2.get_yticks()[i] for i in range(0, len(ax2.get_yticks()),2)])
    plt.savefig(f'{output_prefix}_hyb_efficiency_median_bar.svg')
        

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


# Figure 3E
# KS test (2 sample) for targeting efficacy distributions comparing host to virus within each replicate
# Paired T-test for median targeting efficacies (N = 3) 

plot_hyb_eff_violin(snu_small_frac_paths, snu_long_frac_paths, snu_mrna_clash_paths, min_hyb_counts=100, min_mir_counts=20, min_mrna_tpm=5, virus='ebv', output_prefix='snu719')
plot_hyb_eff_violin(akata_small_frac_paths, akata_long_frac_paths, akata_mrna_clash_paths, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=5, virus='ebv', output_prefix='akata')
plot_hyb_eff_violin(mhv68_small_frac_paths, mhv68_long_frac_paths, mhv68_mrna_clash_paths, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=5, virus='mghv', output_prefix='mhv68')
plot_hyb_eff_violin(kshv_small_frac_paths, kshv_long_frac_paths, kshv_mrna_clash_paths, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=5, virus='kshv', output_prefix='kshv')

###################################################################