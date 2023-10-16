import pandas as pd
import numpy as np
from random import randint
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr, spearmanr, mannwhitneyu, ks_2samp, ttest_rel, sem
import glob



def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'MTRNR' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    df = df.groupby(['mir', 'mrna', 'mrna_feature_start'])['count'].sum().sort_values().reset_index()
    df.index = [f'{i}_{j}_{k}' for i,j,k in zip(df['mir'], df['mrna'], df['mrna_feature_start'])]
    return df



def process_mrna(paths):

    ind = set()
    for i in paths:
        x = pd.read_table(i, index_col=0)
        ind = ind | set(x.index)
    new = pd.DataFrame(index=ind)
    for i in paths:
        x = pd.read_table(i, index_col=0)
        new[i] = x['tpm']
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

def process_clash(paths, vir='ebv'):


    ind = set()
    for i in paths:

        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'MTRNR' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'mrna_feature_start':lambda x:str(x[0])}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j,k]) for i,j,k in zip(df['mir'], df['mrna'],df['mrna_feature_start'])]                                                                              

        ind = ind | set(df.index)
    new = pd.DataFrame(index=ind)

    for i in paths:
        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'MTRNR' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'mrna_feature_start':lambda x:str(x[0])}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j,k]) for i,j,k in zip(df['mir'], df['mrna'], df['mrna_feature_start'])]                                                                              
        new[i] = df['count']
    new = pd.DataFrame(np.mean(new,1))
    new['mrna'] = new.index.map(lambda x: x.split('_')[1])
    new['mir'] = new.index.map(lambda x: x.split('_')[0])
    new['mrna_feature_start'] = new.index.map(lambda x: x.split('_')[2])
    new['count'] = new[0]
    return new

def integrate(small_frac_paths, mrna_expr_paths, clash_paths, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=1, virus='ebv'):
    
    mir = get_mir_means(small_frac_paths)
    mrna = get_mrna_means(mrna_expr_paths)
    df = process_df(clash_paths[0], vir=virus)
    print(df)
    #df = get_seed(df)
    #df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min, 'seed': max, 'species': lambda x:x[0], 'mrna_feature_start': lambda x:x[0]}).sort_values('count').reset_index().set_index('mrna')
    #df = df.groupby(['mir','mrna'])['count'].sum().reset_index().set_index('mrna')
    # df = df.groupby(['mir', 'mrna']).agg({'count':sum, 'mrna_feature_start': lambda x:x[0]}).reset_index()
    df = df.set_index('mrna')
    df['mrna_exp'] = mrna
    df = df.reset_index().set_index('mir')
    df['mir_exp'] = mir
    clash_means = process_clash(clash_paths)
    df = df.reset_index()
    df.index = ['_'.join([i,j,k]) for i,j,k in zip(df['mir'], df['mrna'],df['mrna_feature_start'])]    
    df['count'] = clash_means 
    df = df[(df['count'] >= min_hyb_counts) & (df['mir_exp'] >= min_mir_counts) & (df['mrna_exp'] >= min_mrna_tpm)] # Filter out low abundance hybs, mirs, and mrnas
    df['hyb'] = df['count']/ (df['mir_exp'] * df['mrna_exp'])
    return df

snu_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/smallfrac/*')
snu_long_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/longfrac/*')
snu_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*shuffled_dg.tsv')
snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')


cl = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/blood871418-suppl2.csv', sep=',')
cl = cl.set_index('Patient barcode')
mir = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/bl_mir_cpm_primary_tumors_only.tsv', index_col=0)


bl = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/bl_gene_counts_barcode.tsv', index_col=0)
discovery = set(cl[cl['Cohort']=='Discovery'].index)
ebv_pos = set(bl.loc['RPMS1'].sort_values().index[18:])
mirs = set(mir.columns)

colnames = list(discovery & ebv_pos & mirs)

bl = bl[colnames]
bl = 1000000*bl/np.sum(bl)
mir = mir[colnames]

stad = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stad_mss_only_211neg_24pos.fordeseq.tsv', index_col=0)
stad = 1000000*stad/np.sum(stad) 
stadmirs = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stad_mirs.tsv', index_col=0)

stad = stad[stad.columns[-24:]]
cols = set(stad.columns) & set(stadmirs.columns)
stad = stad[cols]
stadmirs = stadmirs[cols]



min_counts = 100
# min_eff = 3
min_mrna = 5
min_mean_tumor_mrna = 5


#clash = process_clash(snu_mrna_clash_paths)
clash = integrate(snu_small_frac_paths, snu_long_frac_paths, snu_mrna_clash_paths, min_hyb_counts=100, min_mir_counts=20, min_mrna_tpm=5, virus='ebv')
clash = clash.reset_index().set_index('mrna')
snu_exp = process_mrna(snu_long_frac_paths)     

stad = stad[np.mean(stad, 1) > min_mean_tumor_mrna]
clash['mrna_exp'] = snu_exp
clash = clash.dropna()
clash['eff'] = clash['count'] / clash['mrna_exp'] 
clash = clash.reset_index()
# clash = clash[(clash['count'] > min_counts) & (clash['mrna_exp'] > min_mrna) & (clash['eff'] > min_eff)]
clash = clash[clash['mrna'].isin(stad.index)]
clash = clash[clash['mir'].isin(stadmirs.index)]

ebv = clash[clash['mir'].str.contains('ebv')]
utr3 = ebv[ebv['mrna_feature_start']=="3'UTR"]
utr5 = ebv[ebv['mrna_feature_start']=="5'UTR"]
cds = ebv[ebv['mrna_feature_start']=="CDS"]
hsa = clash[clash['mir'].str.contains('hsa')]
utr3h = hsa[hsa['mrna_feature_start']=="3'UTR"]
utr5h = hsa[hsa['mrna_feature_start']=="5'UTR"]
cdsh = hsa[hsa['mrna_feature_start']=="CDS"]



#STAD
cor = []
for mi, mr in zip(utr3['mir'], utr3['mrna']):
    cor.append(spearmanr(np.log2(stad.loc[mr]+1), np.log2(stadmirs.loc[mi]+1))[0])
utr3['stad_corr'] = cor


n_permute = 100
stad_random_cor_e = np.zeros(n_permute * len(utr3['mir']))
i = 0
for c in range(n_permute):
    for mi in utr3['mir']:
        rand_index = randint(0, len(stad.index) - 1)
        stad_random_cor_e[i] = spearmanr(np.log2(stad.iloc[rand_index]+1), np.log2(stadmirs.loc[mi]+1))[0]
        i += 1


print(ks_2samp(utr3['stad_corr'], stad_random_cor_e, alternative="greater"), np.median(utr3['stad_corr']), np.median(stad_random_cor_e))


cor = []
for mi, mr in zip(cds['mir'], cds['mrna']):
    cor.append(spearmanr(np.log(stad.loc[mr]+1), np.log(stadmirs.loc[mi]+1))[0])
cds['stad_corr'] = cor

n_permute = 100
stad_random_cor_e_cds = np.zeros(n_permute * len(cds['mir']))
i = 0
for c in range(n_permute):
    for mi in cds['mir']:
        rand_index = randint(0, len(stad.index) - 1)
        stad_random_cor_e_cds[i] = spearmanr(np.log2(stad.iloc[rand_index]+1), np.log2(stadmirs.loc[mi]+1))[0]
        i += 1

print(ks_2samp(cds['stad_corr'], stad_random_cor_e_cds, alternative="greater"), np.median(cds['stad_corr']), np.median(stad_random_cor_e_cds))


cor = []
for mi, mr in zip(utr5['mir'], utr5['mrna']):
    cor.append(spearmanr(np.log(stad.loc[mr]+1), np.log(stadmirs.loc[mi]+1))[0])
utr5['stad_corr'] = cor

n_permute = 100
stad_random_cor_e_utr5 = np.zeros(n_permute * len(utr5['mir']))
i = 0
for c in range(n_permute):
    for mi in utr5['mir']:
        rand_index = randint(0, len(stad.index) - 1)
        stad_random_cor_e_utr5[i] = spearmanr(np.log2(stad.iloc[rand_index]+1), np.log2(stadmirs.loc[mi]+1))[0]
        i += 1

print(ks_2samp(utr5['stad_corr'], stad_random_cor_e_utr5, alternative="greater"), np.median(utr5['stad_corr']), np.median(stad_random_cor_e_utr5))



#human
cor = []
for mi, mr in zip(utr3h['mir'], utr3h['mrna']):
    cor.append(spearmanr(np.log2(stad.loc[mr]+1), np.log2(stadmirs.loc[mi]+1))[0])
utr3h['stad_corr'] = cor


n_permute = 100
stad_random_cor_h = np.zeros(n_permute * len(utr3h['mir']))
i = 0
for c in range(n_permute):
    for mi in utr3h['mir']:
        rand_index = randint(0, len(stad.index) - 1)
        stad_random_cor_h[i] = spearmanr(np.log2(stad.iloc[rand_index]+1), np.log2(stadmirs.loc[mi]+1))[0]
        i += 1

print(ks_2samp(utr3h['stad_corr'], stad_random_cor_h, alternative="greater"), np.median(utr3h['stad_corr']), np.median(stad_random_cor_h))


cor = []
for mi, mr in zip(cdsh['mir'], cdsh['mrna']):
    cor.append(spearmanr(np.log2(stad.loc[mr]+1), np.log2(stadmirs.loc[mi]+1))[0])
cdsh['stad_corr'] = cor


n_permute = 100
stad_random_cor_cdsh = np.zeros(n_permute * len(cdsh['mir']))
i = 0
for c in range(n_permute):
    for mi in cdsh['mir']:
        rand_index = randint(0, len(stad.index) - 1)
        stad_random_cor_cdsh[i] = spearmanr(np.log2(stad.iloc[rand_index]+1), np.log2(stadmirs.loc[mi]+1))[0]
        i += 1


print(ks_2samp(cdsh['stad_corr'], stad_random_cor_cdsh), np.median(cdsh['stad_corr']), np.median(stad_random_cor_cdsh))


cor = []
for mi, mr in zip(utr5h['mir'], utr5h['mrna']):
    cor.append(spearmanr(np.log2(stad.loc[mr]+1), np.log2(stadmirs.loc[mi]+1))[0])
utr5h['stad_corr'] = cor


n_permute = 100
stad_random_cor_utr5h = np.zeros(n_permute * len(utr5h['mir']))
i = 0
for c in range(n_permute):
    for mi in utr5h['mir']:
        rand_index = randint(0, len(stad.index) - 1)
        stad_random_cor_utr5h[i] = spearmanr(np.log2(stad.iloc[rand_index]+1), np.log2(stadmirs.loc[mi]+1))[0]
        i += 1


print(ks_2samp(utr5h['stad_corr'], stad_random_cor_utr5h), np.median(utr5h['stad_corr']), np.median(stad_random_cor_utr5h))


# n_permute = 1
# random_cor_e = np.zeros(n_permute * len(utr3['mir']))
# i = 0
# for c in range(n_permute):
#     for mi in utr3['mir']:
#         rand_index = randint(0, len(stad.index) - 1)
#         random_cor_e[i] = pearsonr(np.log(stad.iloc[rand_index]+1), np.log(stadmirs.loc[mi]+1))[0]
#         i += 1

# n_permute = 1
# random_cor_h = np.zeros(n_permute * len(utr3h['mir']))
# i = 0
# for c in range(n_permute):
#     for mi in utr3h['mir']:
#         rand_index = randint(0, len(bl.index) - 1)
#         random_cor_h[i] = pearsonr(np.log(stad.iloc[rand_index]+1), np.log(stadmirs.loc[mi]+1))[0]
#         i += 1


nbins = 2000
ax = plt.subplot()
n, bins, patches = ax.hist(utr3['stad_corr'], nbins, density=True, histtype='step',
                           cumulative=True, label="EBV 3'UTR", lw=1,color ='red')


#n, bins, patches = ax.hist(cds['stad_corr'], nbins, density=True, histtype='step',
#                           cumulative=True, label="EBV CDS", lw=.5)

n, bins, patches = ax.hist(utr3h['stad_corr'], nbins, density=True, histtype='step',
                           cumulative=True, label="Human 3'UTR", lw=1,color='#B2B2B2')


#n, bins, patches = ax.hist(cdsh['stad_corr'], nbins, density=True, histtype='step',
 #                          cumulative=True, label="Human CDS", lw=.5)


n, bins, patches = ax.hist(stad_random_cor_e, nbins, density=True, histtype='step',
                           cumulative=True, label="Random (EBV 3'UTR)", lw=1,ls=':',color='k')

#n, bins, patches = ax.hist(stad_random_cor_e_cds, nbins, density=True, histtype='step',
 #                          cumulative=True, label='Random (EBV CDS)', lw=.5,ls='--', color='.75')

# n, bins, patches = ax.hist(utr5['bl_corr'], nbins, density=True, histtype='step',
#                            cumulative=True, label="EBV 5'UTR", lw=2)


n, bins, patches = ax.hist(stad_random_cor_h, nbins, density=True, histtype='step',
                           cumulative=True, label="Random (Human 3'UTR)", lw=1,ls='--',color='k')

#n, bins, patches = ax.hist(stad_random_cor_cdsh, nbins, density=True, histtype='step',
 #                          cumulative=True, label='Random (Human CDS)', lw=.5,ls='--', color='.95')


# n, bins, patches = ax.hist(utr5['bl_corr'], nbins, density=True, histtype='step',
#                            cumulative=True, label="EBV 5'UTR", lw=2)



#plt.legend()
plt.savefig('/Users/nate/Projects/EBV_interactome/stad_cdf_plot.svg')





#BL

clash = process_clash(snu_mrna_clash_paths)
clash = clash.reset_index().set_index('mrna')
snu_exp = process_mrna(snu_long_frac_paths)     
bl = bl[np.mean(bl, 1) > min_mean_tumor_mrna]
clash['mrna_exp'] = snu_exp
clash = clash.dropna()
clash['eff'] = clash['count'] / clash['mrna_exp'] 
clash = clash.reset_index()
clash = clash[(clash['count'] > min_counts) & (clash['mrna_exp'] > min_mrna) & (clash['eff'] > min_eff)]
clash = clash[clash['mrna'].isin(bl.index)]
clash = clash[clash['mir'].isin(mir.index)]

ebv = clash[clash['mir'].str.contains('ebv')]
utr3 = ebv[ebv['mrna_feature_start']=="3'UTR"]
utr5 = ebv[ebv['mrna_feature_start']=="5'UTR"]
cds = ebv[ebv['mrna_feature_start']=="CDS"]
hsa = clash[clash['mir'].str.contains('hsa')]
utr3h = hsa[hsa['mrna_feature_start']=="3'UTR"]
utr5h = hsa[hsa['mrna_feature_start']=="5'UTR"]
cdsh = hsa[hsa['mrna_feature_start']=="CDS"]




cor = []
for mi, mr in zip(utr3['mir'], utr3['mrna']):
    cor.append(spearmanr(np.log2(bl.loc[mr]+1), np.log2(mir.loc[mi]+1))[0])
utr3['bl_corr'] = cor


n_permute = 100
bl_random_cor_e = np.zeros(n_permute * len(utr3['mir']))
i = 0
for c in range(n_permute):
    for mi in utr3['mir']:
        rand_index = randint(0, len(bl.index) - 1)
        bl_random_cor_e[i] = spearmanr(np.log2(bl.iloc[rand_index]+1), np.log2(mir.loc[mi]+1))[0]
        i += 1


print(ks_2samp(utr3['bl_corr'], bl_random_cor_e, alternative="greater"), np.median(utr3['bl_corr']), np.median(bl_random_cor_e))


cor = []
for mi, mr in zip(cds['mir'], cds['mrna']):
    cor.append(spearmanr(np.log(bl.loc[mr]+1), np.log(mir.loc[mi]+1))[0])
cds['bl_corr'] = cor

n_permute = 100
bl_random_cor_e_cds = np.zeros(n_permute * len(cds['mir']))
i = 0
for c in range(n_permute):
    for mi in cds['mir']:
        rand_index = randint(0, len(bl.index) - 1)
        bl_random_cor_e_cds[i] = spearmanr(np.log2(bl.iloc[rand_index]+1), np.log2(mir.loc[mi]+1))[0]
        i += 1

print(ks_2samp(cds['bl_corr'], bl_random_cor_e_cds, alternative="greater"), np.median(cds['bl_corr']), np.median(bl_random_cor_e_cds))


cor = []
for mi, mr in zip(utr5['mir'], utr5['mrna']):
    cor.append(spearmanr(np.log(bl.loc[mr]+1), np.log(mir.loc[mi]+1))[0])
utr5['bl_corr'] = cor

n_permute = 100
bl_random_cor_e_utr5 = np.zeros(n_permute * len(utr5['mir']))
i = 0
for c in range(n_permute):
    for mi in utr5['mir']:
        rand_index = randint(0, len(bl.index) - 1)
        bl_random_cor_e_utr5[i] = spearmanr(np.log2(bl.iloc[rand_index]+1), np.log2(mir.loc[mi]+1))[0]
        i += 1

print(ks_2samp(utr5['bl_corr'], bl_random_cor_e_utr5), np.median(utr5['bl_corr']), np.median(bl_random_cor_e_utr5))



#human
cor = []
for mi, mr in zip(utr3h['mir'], utr3h['mrna']):
    cor.append(spearmanr(np.log2(bl.loc[mr]+1), np.log2(mir.loc[mi]+1))[0])
utr3h['bl_corr'] = cor


n_permute = 100
bl_random_cor_h = np.zeros(n_permute * len(utr3h['mir']))
i = 0
for c in range(n_permute):
    for mi in utr3h['mir']:
        rand_index = randint(0, len(bl.index) - 1)
        bl_random_cor_h[i] = spearmanr(np.log2(bl.iloc[rand_index]+1), np.log2(mir.loc[mi]+1))[0]
        i += 1

print(ks_2samp(utr3h['bl_corr'], bl_random_cor_h, alternative="greater"), np.median(utr3h['bl_corr']), np.median(bl_random_cor_h))


cor = []
for mi, mr in zip(cdsh['mir'], cdsh['mrna']):
    cor.append(spearmanr(np.log2(bl.loc[mr]+1), np.log2(mir.loc[mi]+1))[0])
cdsh['bl_corr'] = cor


n_permute = 100
bl_random_cor_cdsh = np.zeros(n_permute * len(cdsh['mir']))
i = 0
for c in range(n_permute):
    for mi in cdsh['mir']:
        rand_index = randint(0, len(bl.index) - 1)
        bl_random_cor_cdsh[i] = spearmanr(np.log2(bl.iloc[rand_index]+1), np.log2(mir.loc[mi]+1))[0]
        i += 1


print(ks_2samp(cdsh['bl_corr'], bl_random_cor_cdsh), np.median(cdsh['bl_corr'], alternative="greater"), np.median(bl_random_cor_cdsh))


cor = []
for mi, mr in zip(utr5h['mir'], utr5h['mrna']):
    cor.append(spearmanr(np.log2(bl.loc[mr]+1), np.log2(mir.loc[mi]+1))[0])
utr5h['bl_corr'] = cor


n_permute = 100
bl_random_cor_utr5h = np.zeros(n_permute * len(utr5h['mir']))
i = 0
for c in range(n_permute):
    for mi in utr5h['mir']:
        rand_index = randint(0, len(bl.index) - 1)
        bl_random_cor_utr5h[i] = spearmanr(np.log2(bl.iloc[rand_index]+1), np.log2(mir.loc[mi]+1))[0]
        i += 1


print(ks_2samp(utr5h['bl_corr'], bl_random_cor_utr5h), np.median(utr5h['bl_corr']), np.median(bl_random_cor_utr5h))



# alpha = 0.05
#threshold = np.percentile(random_cor_e, 5)


nbins = 2000
ax = plt.subplot()
n, bins, patches = ax.hist(utr3['bl_corr'], nbins, density=True, histtype='step',
                           cumulative=True, label="EBV 3'UTR", lw=1.5, color='red')


#n, bins, patches = ax.hist(cds['bl_corr'], nbins, density=True, histtype='step',
 #                          cumulative=True, label="EBV CDS", lw=1)

n, bins, patches = ax.hist(utr3h['bl_corr'], nbins, density=True, histtype='step',
                           cumulative=True, label="Human 3'UTR", lw=1.5, color='#B2B2B2')


#n, bins, patches = ax.hist(cdsh['bl_corr'], nbins, density=True, histtype='step',
 #                          cumulative=True, label="Human CDS", lw=1)


n, bins, patches = ax.hist(bl_random_cor_e, nbins, density=True, histtype='step',
                           cumulative=True, label="Random (EBV 3'UTR)", lw=1,ls=':',color='k')

#n, bins, patches = ax.hist(bl_random_cor_e_cds, nbins, density=True, histtype='step',
 #                          cumulative=True, label='Random (EBV CDS)', lw=1)

# n, bins, patches = ax.hist(utr5['bl_corr'], nbins, density=True, histtype='step',
#                            cumulative=True, label="EBV 5'UTR", lw=2)


n, bins, patches = ax.hist(bl_random_cor_h, nbins, density=True, histtype='step',
                           cumulative=True, label="Random (Human 3'UTR)", lw=1, ls='--',color='k')

#n, bins, patches = ax.hist(bl_random_cor_cdsh, nbins, density=True, histtype='step',
 #                          cumulative=True, label='Random (Human CDS)', lw=1)


# n, bins, patches = ax.hist(utr5['bl_corr'], nbins, density=True, histtype='step',
#                            cumulative=True, label="EBV 5'UTR", lw=2)



#plt.legend()
plt.savefig('/Users/nate/Projects/EBV_interactome/bl_cdf_plot.svg')






utr3_100 = utr3[utr3['count'] > 100]
cds_100 = cds[cds['count'] > 100]
df = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/discovery_cohort_spearman.tsv', index_col=0) 
mrnas = set(utr3_100[utr3_100['corr'] < -.27]['mrna']) | set(cds_100[cds_100['corr'] < -.27]['mrna'])
mm = df[set(mrnas) & set(df.columns)]



from gprofiler import GProfiler 
gp = GProfiler(user_agent='ExampleTool', return_dataframe=True )   

a2 = gp.profile(query=list(mm['CSDE1'].sort_values().index[-500:]),sources=['GO:BP','GO:MP','REAC','KEGG','GO:CC']) 




#spearman
cor = []
for mi, mr in zip(utr3['mir'], utr3['mrna']):
    cor.append(spearmanr(bl.loc[mr], mir.loc[mi])[0])
utr3['corr'] = cor

cor = []
for mi, mr in zip(cds['mir'], cds['mrna']):
    cor.append(spearmanr(bl.loc[mr], mir.loc[mi])[0])
cds['corr'] = cor


min_counts = 25
utr3_100 = utr3[utr3['count']>min_counts]
cds_100 = cds[cds['count'] > min_counts]


# alpha = 0.05
np.percentile(random_cor_e, 5)

n_permute = 100
random_cor = np.zeros(n_permute * len(utr3_100['mir']))
i = 0
for c in range(100):
    for mi in utr3_100['mir']:
        rand_index = randint(0, len(bl.index) - 1)
        random_cor[i] = spearmanr(bl.iloc[rand_index], mir.loc[mi])[0]
        i += 1


mut = pd.read_table('blood8871418-suppl2_mutations.csv', index_col=0)
mut
pos=mut.loc[bl.columns]
mut.index = [i.rsplit('-',2)[0] for i in mut.index]
mut
pos=mut.loc[bl.columns]
neg=mut.loc[set(i for i in mut.index if i not in bl.columns)]
neg
pos
neg.iloc[3]
neg.iloc[3][:10]
neg.iloc[3][:100]
neg.iloc[3][:50]
neg.iloc[3][50:]
neg.iloc[3][50:100]
neg.iloc[3]['SYMBOL']
neg.value_counts("SYMBOL")
neg['SYMBOL'].value_counts()
neg['SYMBOL'].value_counts().sort_values()
neg = neg.reset_index()
pos = pos.reset_index()
neg.iloc[1]
neg.groupby(['index', 'SYMBOL']).unique()
neg.groupby(['index', 'SYMBOL']).count()
neg.groupby(['index', 'SYMBOL'])['SYMBOL'].count()
neg.groupby(['index', 'SYMBOL'])['SYMBOL'].count().sort_values()
neg.groupby(['index', 'SYMBOL'])['SYMBOL'].count().sort_values().reset_index()
neg.groupby(['index', 'SYMBOL'])['SYMBOL'].count().sort_values()
pd.DataFrame(neg.groupby(['index', 'SYMBOL'])['SYMBOL'].count().sort_values())
pd.DataFrame(neg.groupby(['index', 'SYMBOL'])['SYMBOL'].count().sort_values()).reset_index()
pd.DataFrame(neg.groupby(['index', 'SYMBOL']).count().sort_values())
neg.groupby(['index', 'SYMBOL']).count().reset_index()
neg.groupby(['index', 'SYMBOL']).count().reset_index()[['index', 'SYMBOL']
]
nn=neg.groupby(['index', 'SYMBOL']).count().reset_index()[['index', 'SYMBOL']]
nn
nn.value_counts('SYMBOL').count()
nn.groupby('SYMBOL').count()
nn.groupby('SYMBOL').count().sort_values()
nn.groupby('SYMBOL').count().sort_values('SYMBOL')
nn.groupby('SYMBOL').count().sort_values('index')
pp=pos.groupby(['index', 'SYMBOL']).count().reset_index()[['index', 'SYMBOL']]
nn=nn.groupby('SYMBOL').count().sort_values('index')
pp=pp.groupby('SYMBOL').count().sort_values('index')
pp
nn
pp
nn
len(set(neg.index))
neg.index
len(set(neg['index']))
len(set(pos['index']))
nn
nn/25
nn=nn/25
pp=pp/66
nn
pp
nn
pp
new = pd.DataFrame(index = set(nn.index) & set(pp.index))
new
new['pos'] = pp
new['neg'] = nn
new = new.fillna(0)
new
new.sort_values('pos')
new = new.sort_values('pos')
new.sort_values('neg')
new['fc'] = (new['neg'] +.01)/(new['pos'] + .01)
new
new = new.sort_values('fc')
new
new[new['neg'] > .15]
new[new['pos'] > .15]
