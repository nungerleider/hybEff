from itertools import product
from collections import defaultdict
import subprocess


def rev_comp(seq): 
    dna = {'U':'A','A':'T','C':'G','G':'C','T':'A', 'N':'N'} 
    return ''.join(reversed([dna[i.upper()] for i in seq])) 

def count_motif(): 

    fa_dict = load_fa(fa)
    print('fa loaded')
    bases = list('ACTG')
    n = [bases for i in range(kmer_length)]
    motif_counts = {rev_comp(''.join(i)) +'A': 0 for i in product(*n)}
    penalty = [''] * 8
    seq_id_motifs = defaultdict(list)

    for seq_id in fa_dict: 
        seq = fa_dict[seq_id] 
        i = 0 
        while i < len(seq) - 9: 
            motif = seq[i:i+8]
            if motif[-1] != 'A':
                i+=1
                continue
            if motif not in penalty: 
                motif_counts[motif]+= 1 
                #seq_id_motifs[motif] += seq_id
            penalty += motif
            penalty.pop(0)
            i += 1 
        print(seq_id)
    

    return motif_counts      


def load_fa(path):
    fa_dict = defaultdict(str)
    with open(path) as fa_in:
        for line in fa_in:
            line = line.strip('\n')
            if line[0] == '>':
                seq_id = line.strip('>')
            else:
                fa_dict[seq_id.split('.')[0]]+=line
    return fa_dict


def get_mir(path='/Users/nate/Downloads/mature.fa'):
    mir = {} 
    with open(path) as infile: 
        for line in infile: 
            if line[0] == '>': #and 'hsa' in line or 'ebv' in line: 
                seq_id = line.strip('>\n') 
            else: 
                mir[seq_id.split(' ')[0]] = line.strip('\n') 
    return mir
                            

kmer_length = 7
fa_path = '/Users/nate/Downloads/gencode.v34.transcripts.fa'
fa = load_fa(fa_path)
motif_counts = count_motif()


seq = 'ACTACAGAGAGAGAGAG'  # Sample sequence
folding = subprocess.check_output(f'echo {seq} | RNAplfold -u 14 -L 40 -W 80', shell=True).decode('utf8')                              
x = pd.read_table('plfold_lunp', skiprows=2, header=None, index_col=0)