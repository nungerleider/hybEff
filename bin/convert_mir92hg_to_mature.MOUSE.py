import sys

input_path = sys.argv[1] #180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_2_1_comp_human_ebv.clash.blast.natescript.hyb


hg_dict = {

'GCAAATCTAT': 'mmu-miR-19a-3p',
'CACTTGTCCCGGC': 'mmu-miR-92a-3p',
'GCATTACGAGCA': 'mmu-miR-20a-3p',
'TAAAGTGCTTATAGT': 'mmu-miR-20a-5p',
'CAAAGTGCTTACAG': 'mmu-miR-17-5p',
'CTGCAGTGAGGG': 'mmu-miR-17-3p',
'ACCTATGAATTGACAGC': 'mmu-miR-192-5p',
'CAGCAACTCCA': 'mmu-miR-194-5p',
'AGCACCATTTGAAATCAG': 'mmu-miR-29b-3p',
'ATCTGAAATCG': 'mmu-miR-29a-3p',
'CCATTTGAAATCGG': 'mmu-miR-29c-3p',
'GTGTCTTAGCTG': 'mmu-miR-34a-5p',
'CAGCAAGTATACTGC': 'mmu-miR-34a-3p',

'AGCACAGAAATATT': 'mmu-miR-195-5p',
'GGTAGTAGGTTGTGTG': 'mmu-let-7b-5p',
'CTGCCAGTTGAAGAAC': 'mmu-miR-22-3p',
'TACAATCTACTGTCTT': 'mmu-let-7a-3p',



}

with open(input_path) as infile, open(input_path + '.mir92hg_converted.tsv', 'w') as outfile:
    for line in infile:
        if 'Mir' not in line:
            outfile.write(line)
            continue
        line = line.split('\t')
        for key in hg_dict:
            if key not in line[-1]:
                continue
            else:
                if 'Mir' in line[1]:
                    line[1] = hg_dict[key]
                    
                elif 'Mir' in line[7]:
                    line[8] = hg_dict[key]
                
                break
        outfile.write('\t'.join(line) + '\n')
