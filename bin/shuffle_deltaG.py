# Shuffle microRNA and mRNA sequences then calculate predicted binding energy 
# Used as a control for true hybrids
# Produces a new clash dataframe with additional columns
# Minimum counts set at 10 to reduce run time. Subsequent analysis always used hybrid pairs w at least 10 counts

import subprocess
import sys
import pandas as pd
from random import shuffle
import numpy as np


def process_clash(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt' in i]]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    return df


def shuffle_word(word):
    word = list(word)
    shuffle(word)
    return ''.join(word)


def run_rna_fold(mir_seq, mrna_seq):
    '''RNAcofold must be in system path'''

    folding = subprocess.check_output('echo "%s&%s" | RNAcofold --noPS' % (mir_seq, mrna_seq), shell=True).decode('utf8')
    folding_structure, free_energy = folding.split('\n')[1].split(' ', 1)
    mirna_half, mrna_half = folding_structure.split('&')
    free_energy = free_energy.lstrip('(').rstrip(')')
    return float(free_energy)


def main(clash_path, virus, min_counts):

    df = process_clash(clash_path, virus)
    df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min, 'mir_seq':lambda x:x[0], 'mrna_seq':lambda x:x[0]}).reset_index()
    df['shuffle_mir'] = df['mir_seq'].map(lambda x:shuffle_word(x))
    df['shuffle_mrna'] = df['mrna_seq'].map(lambda x:shuffle_word(x))
    df = df[df['count'] > min_counts]
    delta_gs = []
    for shuff_mir, shuff_mrna in zip(df['shuffle_mir'], df['shuffle_mrna']):
        delta_g = run_rna_fold(shuff_mir, shuff_mrna)
        delta_gs.append(delta_g)
    df['shuffle_deltaG'] = delta_gs
    df.to_csv(f'{clash_path}.shuffled_dg.tsv', sep='\t')


mhv68_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*annotated')
kshv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*annotated')
akata_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*annotated')
snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')

min_counts = 10
for path in mhv68_mrna_clash_paths:
    main(path, 'mghv', min_counts)

for path in kshv_mrna_clash_paths:
    main(path, 'kshv', min_counts)

for path in akata_mrna_clash_paths:
    main(path, 'ebv', min_counts)

for path in snu_mrna_clash_paths:
    main(path, 'ebv', min_counts)