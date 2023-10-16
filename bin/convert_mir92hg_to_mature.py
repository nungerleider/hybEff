import sys

input_path = sys.argv[1] #180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_2_1_comp_human_ebv.clash.blast.natescript.hyb


hg_dict = {

'GCAAATC': 'hsa-miR-19a-3p',
'CACTTGTCCCGGC': 'hsa-miR-92a-3p',
'GCATTATGAGCACTTA': 'hsa-miR-20a-3p',
'TAAAGTGCTTATAGT': 'hsa-miR-20a-5p',
'CAAAGTGCTTACAG': 'hsa-miR-17-5p',
'CTGCAGTGAAGGCAC': 'hsa-miR-17-3p',
'ACCTATGAATTGACAGC': 'hsa-miR-192-5p',
'CAGCAACTCCA': 'hsa-miR-194-5p',
'TGAAATC': 'hsa-miR-29b-3p',
'GTGTCTTAGCTG': 'hsa-miR-34a-5p',
'AGCACAGAAATATT': 'hsa-miR-195-5p',
'GGTAGTAGGTTGTGTG': 'hsa-let-7b-5p',
'CTGCCAGTTGAAGAAC': 'hsa-miR-22-3p',
'TACAATCTACTGTCTT': 'hsa-let-7a-3p',



}

with open(input_path) as infile, open(input_path + '.mir92hg_converted.tsv', 'w') as outfile:
    for line in infile:
        if 'MIR' not in line:
            outfile.write(line)
            continue
        line = line.split('\t')
        for key in hg_dict:
            if key not in line[-1]:
                continue
            else:
                if 'MIR' in line[1]:
                    line[1] = hg_dict[key]
                    
                elif 'MIR' in line[7]:
                    line[8] = hg_dict[key]
                
                break
        outfile.write('\t'.join(line) + '\n')
