import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# STAD microRNA barplot, percent EBV (Fig 1b?)
m = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stad_mir_cpm.tsv', index_col=0)  # STAD c.p.m. microRNA

m2 = m[np.sum(m.loc[[i for i in m.index if 'ebv' in i]]).sort_values().index]
m2 = m2[m2.columns[-38:]]
sums = np.sum(m2.loc[[i for i in m2.index if 'ebv' in i]]) / 10000  # convert c.p.m. to c.p.hundred for percentage plotting

fig = plt.figure(figsize=(24,5), dpi=300)
ax = plt.subplot()
ax2 = ax.twinx()
ax2.bar([i-.2 for i in range(len(m2.columns))], sums,facecolor='#72DDF7', ec='k',linewidth=.5)
ax2.set_yticks(range(0, 101, 25))
ax.set_xlim([-.9, 37.5])
ax.set_yticks([])
ax.set_xticks([])
ax2.grid(ls='--',alpha=.5)
ax2.set_axisbelow(True)
plt.savefig('stad_sum_ebv_mirs.svg',dpi=300)