import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

hsa_color = '0.5'
ebv_color = '#149911'

s = pd.read_csv('/Users/nate/Projects/EBV_interactome/bl/blood871418-suppl2.csv', index_col=0)
m = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/bl_mir_cpm_primary_tumors_only.tsv', index_col=0)
d = {'pos':[],'neg':[]}

for i,j in zip(s['Patient barcode'], s['EBV status']):
    if 'pos' in j:
        d['pos'].append(i)
    else:
        d['neg'].append(i)

m1 = m[set(d['neg']) & set(m.columns)]
m2 = m[set(d['pos']) & set(m.columns)]
m3 = pd.DataFrame(index=m.index)

m2 = m2[np.sum(m2.loc[[i for i in m2.index if 'ebv' in i]]).sort_values().index]

m3[m1.columns] = m1
m3[m2.columns] = m2

fig = plt.figure(figsize=(24,6), dpi=300)
ax = plt.subplot()
patients = len(m3.columns)


for x,patient in enumerate(m3.columns):
    hsa = m3.loc[[i for i in m3.index if 'hsa' in i]]
    ebv = m3.loc[[i for i in m3.index if 'ebv' in i]]
    plt.scatter([x-.2]*len(hsa), hsa[patient], c=hsa_color, s=4, alpha=.8, linewidth=0)
    plt.scatter([x-.2]*len(ebv), ebv[patient], c=ebv_color, s=4, alpha=.8,linewidth=0)

means = np.sum(m3.loc[[i for i in m3.index if 'ebv' in i]])/10000

ax2 = ax.twinx()
plt.fill_between(range(patients), means, alpha=.1, color=ebv_color, linewidth=0)
ax2.plot(range(patients), means, c=ebv_color, linewidth=.3, ls='-', alpha=.3)

ax2.set_yticks(range(0, 101, 10))
ax2.set_ylim([0, 100])
ax.set_ylim([0, 500000])

ax2.grid(alpha=.7 , ls='--',lw=.2)
ax.set_xticks([])
ax.set_xlim([-1, 84])
