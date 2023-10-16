import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


x = pd.read_csv('/Users/nate/Projects/EBV_interactome/Akata_Zta.summary.genes.csv', index_col=0)
x = x[x.columns[:4]]
ebv = pd.read_table('/Users/nate/Projects/EBV_interactome/ebvgenes', header=None)
ebv = set(ebv[0]) & set(x.index)
ebv = {i for i in ebv if 'miR' not in i}
ncrna = {'RPMS1', 'EBER1', 'EBER2', 'A73'}
coding = {i for i in ebv if i not in ncrna}

coding_exp = np.sum(x.loc[coding])
noncod_exp = np.sum(x.loc[ncrna])




f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
ax.bar(range(4), noncod_exp)
ax2.bar(range(4), noncod_exp)
ax.bar(range(4), coding_exp, bottom=noncod_exp)
ax2.bar(range(4), coding_exp, bottom=noncod_exp)

ax2.set_ylim([0,2000])
ax.set_ylim([2000, 1000000])

plt.savefig('/Users/nate/Projects/EBV_interactome/sinclair_lytic.svg')

latencygenes = {'LMP-2A',
 'LMP-2B',
 'LMP-1',
 'Cp-EBNA2',
 'Cp-EBNA3A',
 'Cp-EBNA3B',
 'Cp-EBNA3C',
 'Cp-EBNA1',
 'EBNA-LP',
 'Qp-EBNA1',
 'EBER1',
 'EBER2',
 'RPMS1',
 'LF3',
 'A73',
 'BHLF1'}



coding_exp = np.sum(x.loc[coding & latencygenes])
noncod_exp = np.sum(x.loc[ncrna & latencygenes])




f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
ax.bar(range(4), noncod_exp)
ax2.bar(range(4), noncod_exp)
ax.bar(range(4), coding_exp, bottom=noncod_exp)
ax2.bar(range(4), coding_exp, bottom=noncod_exp)

ax2.set_ylim([0,2000])
ax.set_ylim([2000, 1000000])

plt.savefig('/Users/nate/Projects/EBV_interactome/sinclair_lytic_latency_genes_only.svg')
