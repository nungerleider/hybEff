import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations 
from scipy import stats
from scipy.cluster.vq import kmeans2
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles

def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'MTRNR' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    return df

def get_counts(path):
    df = process_df(path)
    df = pd.DataFrame(df.groupby(['mrna', 'mir'])['count'].sum()).reset_index()
    df.index = [f'{i}_{j}' for i,j in zip(df['mir'], df['mrna'])]
    return df

def venn(akata_df, snu_df, min_counts=200):

    akata = set(akata_df[akata_df['mean'] >= min_counts].index)
    snu = set(snu_df[snu_df['mean'] >= min_counts].index)
    
    akata_ebv = set([i for i in akata if 'ebv' in i])
    snu_ebv = set([i for i in snu if 'ebv' in i])
    akata_human = akata - akata_ebv
    snu_human = snu - snu_ebv

    print('akata only:', len(akata - snu))
    print('snu only:', len(snu - akata))
    print('akata and snu:', len(snu & akata))

    print('akata only(ebv):', len(akata_ebv - snu_ebv))
    print('snu only(ebv):', len(snu_ebv - akata_ebv))
    print('akata and snu (ebv):', len(akata_ebv & snu_ebv))
    
    print('akata only(human):', len(akata_human - snu_human))
    print('snu only(human):', len(snu_human - akata_human))
    print('akata and snu (human):', len(akata_human & snu_human))

    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(121)
    ax2 = plt.subplot(122)
    plt.axis('scaled')
    
    h = venn2([akata_human, snu_human],  set_labels=None, ax=ax)
    h.get_label_by_id('10').set_text('')
    h.get_patch_by_id('10').set_color('#ECC30B')
    h.get_patch_by_id('10').set_alpha(.5)

    h.get_label_by_id('11').set_text('')
    h.get_patch_by_id('11').set_color('#F37748')
    h.get_patch_by_id('11').set_alpha(.8)
    h.get_label_by_id('01').set_text('')
    h.get_patch_by_id('01').set_color('#D56062')
    h.get_patch_by_id('01').set_alpha(.5)

    venn2_circles([akata_human, snu_human],ax=ax,lw=2)

    h = venn2([akata_ebv, snu_ebv],  set_labels=None, ax=ax2)
    h.get_label_by_id('10').set_text('')
    h.get_patch_by_id('10').set_color('#ECC30B')
    h.get_patch_by_id('10').set_alpha(.5)

    h.get_label_by_id('11').set_text('')
    h.get_patch_by_id('11').set_color('#F37748')
    h.get_patch_by_id('11').set_alpha(.8)
    h.get_label_by_id('01').set_text('')
    h.get_patch_by_id('01').set_color('#D56062')
    h.get_patch_by_id('01').set_alpha(.5)
    venn2_circles([akata_ebv, snu_ebv],ax=ax2,lw=2)
    plt.savefig('clash_venn.svg')

# df.groupby('mir')['arr'].apply(lambda x:np.sum(x))
def make_array(seq):
    arr = np.zeros(20)
    for i in range(20):
        if i > len(seq) - 1:
            continue
        elif seq[i] == '(':
            arr[i] = 1
            
    return arr


snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')
akata_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*annotated')
mghv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*annotated')
kshv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*annotated')

mhv68_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/smallfrac/CL100128477_L1_WHMOUkmcSAAETAASE-545_1.fq.counts.tsv']
mhv68_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/longfrac/HE2-1.genes.tsv']

kshv_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/smallfrac/541_WT-TIVE-1.small_fraction.counts.tsv']
kshv_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/longfrac/TIVE-1.genes.tsv']

akata_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/smallfrac/*')
akata_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/akata/longfrac/Akata.genes.tsv']

snu_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/smallfrac/*')
snu_long_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/longfrac/*')


mghv_mrna_expression_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/longfrac/HE2-1.genes.tsv']
mghv_mir_expression_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/smallfrac/CL100128477_L1_WHMOUkmcSAAETAASE-545_1.fq.counts.tsv']
# mghv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/SRR839524[5-7]*')# + glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/SRR8395250*')  
mghv_mrna_clash_paths.reverse()


kshv_mir_expression_paths = ['/Users/nate/Projects/EBV_interactome/kshv/smallfrac/541_WT-TIVE-1.small_fraction.counts.tsv']
kshv_mrna_expression_paths = ['/Users/nate/Projects/EBV_interactome/kshv/longfrac/TIVE-1.genes.tsv']
kshv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*annotated')



hybrids = set()
for path in akata_mrna_clash_paths:
    counts = get_counts(path)
    hybrids = hybrids | set(counts.index)

akata = pd.DataFrame(index=hybrids)
for path in akata_mrna_clash_paths:
    counts = get_counts(path)
    akata[path.split('/')[-2]] = counts['count']


akata = akata.fillna(0)
akata['mean'] = np.mean(akata, 1)
akata = akata.sort_values('mean')


hybrids = set()
for path in snu_mrna_clash_paths:
    counts = get_counts(path)
    hybrids = hybrids | set(counts.index)

snu = pd.DataFrame(index=hybrids)
for path in snu_mrna_clash_paths:
    counts = get_counts(path)
    snu[path.split('/')[-1]] = counts['count']


snu = snu.fillna(0)
snu['mean'] = np.mean(snu, 1)
snu = snu.sort_values('mean')
snu+=1
akata+=1
for x, y in combinations(snu.columns[:3], 2):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(snu[x], snu[y], s=1, c='k', alpha=.3)
    ax.set_ylim([1, 100000])
    ax.set_xlim([1, 100000])
    ax.set_xscale('log')
    ax.set_yscale('log')
    pearson, p = stats.pearsonr(snu[x], snu[y])
    m, b = np.polyfit(snu[x], snu[y], 1)
    xvals = np.linspace(np.min(snu[x]),np.max(snu[x]),500)

    #plt.plot(xvals, m*xvals + b, c='r')
    plt.savefig(f'x-{x}_y-{y}_pearson-{pearson}_p-{p}.png', dpi=500)

for x, y in combinations(akata.columns[:3], 2):
    fig = plt.figure()
    ax = plt.subplot(111)
    akata_pruned = akata[[x,y]]
    akata_max = akata_pruned[np.max(akata_pruned, 1) > 1]
    ax.scatter(akata_max[x], akata_max[y], s=1, c='k', alpha=.3)
    ax.set_ylim([1, 100000])
    ax.set_xlim([1, 100000])
    ax.set_xscale('log')
    ax.set_yscale('log')
    pearson, p = stats.pearsonr(akata_max[x], akata_max[y])
    m, b = np.polyfit(akata_max[x], akata_max[y], 1)
    xval0 = np.percentile(akata_max[x],.1)
    yval0 = xval0 * m + b
    xval1 = np.percentile(akata_max[x],.9)
    yval1 = xval1 * m + b
   # plt.plot([xval0, xval1], [yval0, yval1], c='r')
    plt.savefig(f'x-{x}_y-{y}_pearson-{pearson}_p-{p}.png', dpi=500)

snu.to_csv('snu719.mir_mrna.counts.tsv', sep='\t')
akata.to_csv('akata.mir_mrna.counts.tsv', sep='\t')


venn(akata, snu, 200)

#def plot_kmeans()
# start of kmeans clustering of seed matches
df = process_df(snu_mrna_clash_paths[0])
d = pd.DataFrame(df.groupby(['mir', 'mrna','mir_fold','binding_energy'])['count'].sum())
d = d.reset_index()
d['arr'] = d['mir_fold'].map(lambda x:make_array(x))
d = d[d['mir'].str.contains('hsa')]
f = {}


def mir_bind_arr(mirfold):
    arr = np.zeros(35)
    for i in range(len(mirfold)):
        if mirfold[i] == '(':
            arr[i] = 1
    return np.sum(arr)

counts = dict(d['mir'].value_counts())
for i,j in counts.items():
    f[i] = np.zeros([j,20])

for i,j in counts.items():
    mirs = d[d['mir']==i]
    mirs = mirs.reset_index()
    for row in mirs.index:
        f[i][row] = mirs.loc[row, 'arr']

d['arr'] = d['arr'] * abs(d['binding_energy'])

z = np.zeros([len(d.index), 20])
d = d.reset_index()
for i in d.index:
    z[i] = d.loc[i,'arr']

k_clusters = 3
centroid, label = kmeans2(z, k_clusters, minit='points', iter=1000)
w_list = [z[label == i] for i in range(k_clusters)]
# w0 = z[label == 0]
# w1 = z[label == 1]
# w2 = z[label == 2]
# w3 = z[label == 3]
# w4 = z[label == 4]
# for l in range(len(w_list)):
#     plt.matshow(w_list[l], aspect=.001, cmap='binary')
#     plt.savefig(f'{l}_unsorted.svg')
# y = np.concatenate(w_list,0)
# #y=y[:,:23]
# plt.matshow(y, aspect=.001, cmap='binary')
# plt.savefig('a.svg')
for l in range(len(w_list)):
    plt.matshow(w_list[l], aspect=.001, cmap='binary')
    plt.savefig(f'{l}hsa_unsorted.svg')
    w_list[l].sort(0)
    plt.matshow(w_list[l], aspect=.001, cmap='binary')
    plt.savefig(f'{l}hsa_sorted.svg')


y = np.concatenate(w_list, 0)
# y=y[:,:20]
plt.matshow(y, aspect=.001, cmap='binary')
plt.savefig('b.svg')
