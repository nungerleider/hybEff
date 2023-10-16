from collections import defaultdict
import pandas as pd
import sys

blastfile = sys.argv[1] # 180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_2_1_comp_human_ebv.clash.blast.natescript.hyb.mir92hg_converted.tsv
output = blastfile + '.ago_counts.tsv'
d = defaultdict(int)
with open(blastfile) as infile:
    for line in infile:
        line = line.split('\t')
        if len(line) < 9:
            continue
        for field in [1,8]:
            if 'miR' in line[field] or 'let' in line[field]:
                d[line[field]] += int(line[0].split('_')[1])

df = pd.DataFrame.from_dict(d, orient='index').sort_values(0)
df.columns = ['ago_bound_counts']
df.to_csv(output, sep='\t')