from collections import namedtuple, defaultdict
import sys

#gtf = sys.argv[1]
fasta = sys.argv[1]
blast_path = sys.argv[2]


# Map fasta sequence to sequence ID 
# fasta ='/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Papers/The_EBV_microRNA_targetome/clash/alignments/ebv/Akata_Latency1/Latency1_R1_comp_human_ebv.clash_hybrids.fasta'

fa = open(fasta)
fa_d = {}
line_n = 0
for line in fa:
    line = line.strip('\n>')
    if line_n % 2 == 0:
        ID = line
    else:
        fa_d[ID] = line
    line_n+=1

# DONT DO MTOPHITS -> DO CLASH.BLAST
# blast_path = '/home/nate/lat1.blast'

blast_handle = open(blast_path)
hits = defaultdict(list)
Blast = namedtuple('Blast', ['ID', 'gene', 'percent_aligned', 'alignment_length', 'mm0', 'mm1', 'read_start_pos', 'read_stop_pos', 'transcript_start_pos', 'transcript_stop_pos', 'p', 'score'])
for line in blast_handle:
    features = Blast(*line.strip('\n').split('\t'))
    ID, gene, read_start_pos, read_stop_pos, score = features.ID, features.gene, int(features.read_start_pos), int(features.read_stop_pos), float(features.score)
    if 'miR' in gene or 'let' in gene:
        score += 5
    hits[ID].append([ID, gene, read_start_pos, read_stop_pos, score, features.transcript_start_pos, features.transcript_stop_pos])

tophits = defaultdict(list)
for key in hits:
    alignments = hits[key]
    alignments.sort(key=lambda x:int(x[4]), reverse=True)
    top_alignment = alignments[0]
    tophits[key].append(top_alignment)
    n = 1
    while len(tophits[key]) < 2 and n < len(alignments):
        alignment = alignments[n]
        current_start = alignment[2]
        current_stop = alignment[3]
        top_start = top_alignment[2]
        top_stop = top_alignment[3]

        if current_start - top_stop >= -5 or top_start - current_stop >= -5:
            tophits[key].append(alignment)
        n += 1

with open(f'{blast_path}.natescript.hyb', 'w') as blast_out:
      
    for key in tophits:
        if key in fa_d:
            sequence = fa_d[key]
        else:
            sequence = 'Not Found'
        

        if len(tophits[key]) == 1: 
            tophits[key].append(['' for i in range(len(tophits[key][0]))])

        tophits[key].append(sequence)
        alignment1, alignment2, seq = tophits[key] 
        
        blast_out.write('\t'.join([str(i) for i in alignment1]) + '\t')
        blast_out.write('\t'.join([str(i) for i in alignment2]) + '\t')
        blast_out.write(seq + '\n')
