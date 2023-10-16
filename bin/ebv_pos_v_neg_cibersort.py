import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu as mwu
import matplotlib.pyplot as plt

v = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/stadmirs_counts211neg24pos.fordeseq.tsv', index_col=0)
s = pd.read_table('/Users/nate/Projects/EBV_interactome/stad/CIBERSORTx_Job190_Results.csv', index_col=0,sep=',')
s = s.loc[v.columns]


new = pd.DataFrame(index=s.columns[:-4])
fcs, ps, means = [], [], []
for i in new.index:
    fcs.append((np.mean(s.iloc[211:][i])+.01) / (np.mean(s.iloc[:211][i]) +.01))
    ps.append(mwu(s.iloc[211:][i], s.iloc[:211][i], alternative='two-sided')[1])
    means.append(np.mean(s[i] + .01))
    
new['fc'] = fcs
new['ps'] = ps
new['mean'] = means
new = new.sort_values('fc')

fig, ax = plt.subplots(figsize=(6,8))
plt.scatter(np.log2(new['fc']),range(len(new.index)),s= 100 * new['mean'],c=['r' if x < .05 else "0.75" for x in new['ps']],alpha=0.7, lw=0)
ax.set_yticks(range(len(new.index)))
ax.set_yticklabels(new.index)
plt.xlim([-3.5, 3.5])
plt.tight_layout()
plt.savefig('/Users/nate/Projects/EBV_interactome/stad/stad_cib_ebv_posvneg.svg')




bc = pd.read_csv('/Users/nate/Projects/EBV_interactome/bl/CIBERSORTx_Job187_Results.csv', index_col=0)
b = pd.read_table('/Users/nate/Projects/EBV_interactome/bl/revised_bl_ebvneg14ebvpos59_fordeseq.tsv', index_col=0)
b.columns = [i.replace('.', '-') for i in b.columns]
ind = set(bc.index) & set(b.columns)
b = b[ind].loc['RPMS1'].sort_values()
bc = bc.loc[b.index]


new = pd.DataFrame(index=bc.columns[:-4])
fcs, ps, means = [], [], []
for i in new.index:
    fcs.append((np.mean(bc.iloc[12:][i])+.01) / (np.mean(bc.iloc[:12][i]) +.01))
    try:
        ps.append(mwu(bc.iloc[12:][i], bc.iloc[:12][i], alternative='two-sided')[1])
        #ps.append(ttest_ind(bc.iloc[12:][i], bc.iloc[:12][i])[1])
    except:
        ps.append(1)
    means.append(np.mean(bc[i] + .01))
    
new['fc'] = fcs
new['ps'] = ps
new['mean'] = means
new = new.sort_values('fc')

fig, ax = plt.subplots(figsize=(6,8))
plt.scatter(np.log2(new['fc']),range(len(new.index)),s= 100 * new['mean'],c=['r' if x < .05 else "0.75" for x in new['ps']],alpha=0.7, lw=0)
ax.set_yticks(range(len(new.index)))
ax.set_yticklabels(new.index)
plt.xlim([-3.5, 3.5])
plt.tight_layout()
plt.savefig('/Users/nate/Projects/EBV_interactome/bl/bl_cib_ebv_posvneg.svg')
