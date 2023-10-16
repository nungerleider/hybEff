import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


x = pd.read_table('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Papers/The_EBV_microRNA_targetome/data_for_figures/figure1/C/bl_coding_transcripts_only.tpms.tsv',index_col=0)
y = 100 * x / np.sum(x)
fig = plt.figure(figsize=(4,8))
ax = plt.subplot()
parts = ax.violinplot(np.sum(y.loc[[i for i in x.index if 'ENST' not in i]]), [0], showmedians=False, showextrema=False, showmeans=False)
for pc in parts['bodies']:
    pc.set_facecolor('#F7A072')
    pc.set_edgecolor('black')
    pc.set_alpha('1')
plt.yscale('linear')
plt.ylim([0, .06])
ax.grid(lw=1, alpha=.3, ls='--')
ax.set_axisbelow(True)
ax.set_xticks([])
plt.savefig('bl_ebv_coding_percentage.svg', dpi=300)