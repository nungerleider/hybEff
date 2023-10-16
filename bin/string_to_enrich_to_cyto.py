import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr, spearmanr, ks_2samp, ttest_rel, sem
import glob
import requests ## python -m pip install requests
import json
import requests
from collections import Counter, defaultdict
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles


def m13p(mirfold, ind): 
    y = 0 
    for i in mirfold[ind-1:]: 
        if i=='(': 
            y+=1 
    return y 


def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    df['m13p'] = df['mir_fold'].map(lambda x:m13p(x,13))

    return df


def get_seed(df):

    mapping = {
        0: "No seed",
        1.5: "8mer mismatch",
        2: "6mer",
        4: "7merM8",
        6: "7merA1",
        8: "8mer",
        1: "No seed, supplemental",
        2.5: "8mer mismatch, supplemental",
        3: "6mer, supplemental",
        5: "7merM8, supplemental",
        7: "7merA1, supplemental",
        9: "8mer, supplemental"
    }

    df = df[(df['a1'] != 'unknown') & (df['a1'] != 'unannotated')]
    a1 = [int(i) for i in df['a1']]
    core = df['m2_5'] + df['m6_7']
    m8 = df['m8']
    supplemental = [True if i > 3 else False for i in df['m13_17']]
    seeds = []
    for i in range(len(df.index)):
        num = 0
        if supplemental[i] is True:
            num += 1      
        if core[i]  == 6:
            num += 2      
            if a1[i] == 1:
                num += 4         
            if m8[i] == 1:
                num += 2      
        elif core[i] < 6 and a1[i] == 1 and m8[i] == 1:
            num += 1.5
        seeds.append(num)
    df["seed"] = seeds
    return df


def map_seed(df):

    mapping = {
        0: "No seed",
        1.5: "8mer mismatch",
        2: "6mer",
        4: "7merM8",
        6: "7merA1",
        8: "8mer",
        1: "No seed, supplemental",
        2.5: "8mer mismatch, supplemental",
        3: "6mer, supplemental",
        5: "7merM8, supplemental",
        7: "7merA1, supplemental",
        9: "8mer, supplemental",
    }

    df['seed'] = df.seed.map(lambda x: mapping(x))
    return df


def process_mir_df(path):

    df = pd.read_table(path, index_col=0, header=None)
    if '*' in df.index:
        df = df.drop('*')  # Unmapped counts
    df = df.loc[set(i for i in df.index if 'miR' in i or 'let' in i)]
    cpm = 1000000 * df / np.sum(df)  # counts per million mapped mirs
    
    return cpm


def get_mir_means(small_frac_paths):
    '''If more than one small frac seq sample, get average c.p.m. across samples'''

    mirs = set()
    for path in small_frac_paths:
        small = process_mir_df(path)
        mirs = mirs | set(small.index)

    mirdf = pd.DataFrame(index=mirs)
    for path in small_frac_paths:
        small = process_mir_df(path)
        mirdf[path] = small[1]

    return np.mean(mirdf, 1)


def get_mrna_means(mrna_expr_paths):
    '''If more than one sample, get average t.p.m. across samples'''

    mrnas = set()
    for path in mrna_expr_paths:
        mrna = pd.read_table(path, index_col=0)
        mrnas = mrnas | set(mrna.index)

    mrna_df = pd.DataFrame(index=mrnas)
    for path in mrna_expr_paths:
        mrna = pd.read_table(path, index_col=0)
        mrna_df[path] = mrna['tpm']

    return np.mean(mrna_df, 1)



def process_clash(paths, vir='ebv'):

    ind = set()
    for i in paths:

        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j]) for i,j in zip(df['mir'], df['mrna'])]                                                                              

        ind = ind | set(df.index)
    new = pd.DataFrame(index=ind)

    for i in paths:
        df = pd.read_table(i, index_col=0)
        blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'Mt-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12', 'Vmp1'])
        df = df[~df['mrna'].isin(blacklist)]
        df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
        df['count'] = 1000000 * df['count'] / np.sum(df['count'])
        df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
        df.index =  ['_'.join([i,j]) for i,j in zip(df['mir'], df['mrna'])]                                                                              

        new[i] = df['count']
    new = np.mean(new,1)
    return new


def integrate(small_frac_paths, mrna_expr_paths, clash_paths, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=1, virus='ebv'):
    
    mir = get_mir_means(small_frac_paths)
    mrna = get_mrna_means(mrna_expr_paths)
    df = process_df(clash_paths[0], vir=virus)
    df = get_seed(df)
    df = df.groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min, 'seed': max, 'species': lambda x:x[0], 'mrna_feature_start': lambda x:x[0]}).sort_values('count').reset_index().set_index('mrna')
    #df = df.groupby(['mir','mrna'])['count'].sum().reset_index().set_index('mrna')
    df['mrna_exp'] = mrna
    df = df.reset_index().set_index('mir')
    df['mir_exp'] = mir
    clash_means = process_clash(clash_paths)
    df = df.reset_index()
    df.index = ['_'.join([i,j]) for i,j in zip(df['mir'], df['mrna'])]                                                                              
    df['count'] = clash_means 
    df = df[(df['count'] >= min_hyb_counts) & (df['mir_exp'] >= min_mir_counts) & (df['mrna_exp'] >= min_mrna_tpm)] # Filter out low abundance hybs, mirs, and mrnas
    df['hyb'] = df['count']/ (df['mir_exp'] * df['mrna_exp'])
    return df


def get_enrich(gene_list, gene_set_library):

    genes_str = '\n'.join(gene_list)

    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'

    description = 'Gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    user_list_id = data['userListId']

    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)

    header = "Rank, Term name, P-value, Z-score, Combined score, Overlapping genes, Adjusted p-value, Old p-value, Old adjusted p-value".split(', ')
    enrichment = pd.DataFrame(data[gene_set_library])
    enrichment.columns = header
    enrichment = enrichment.set_index('Rank')
    enrichment = enrichment[enrichment['Adjusted p-value']< .0001]

    return enrichment


def string_enrich(my_genes, species='9606'):
    # my_genes = ["IPO7"]

    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    request_url = "/".join([string_api_url, output_format, method])
    params = {

        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : 9606, # species NCBI identifier 
        "limit" : 20,
        "additional_network_nodes":10,
        "caller_identity" : "nate" # your app name

    }

    response = requests.post(request_url, data=params)

    interactome = pd.DataFrame([i.split('\t') for i in response.text.strip().split('\n')])
    enrich = get_enrich(set(interactome[3]), 'Reactome_2016')

    return enrich, set(interactome[3])


def get_n_without_dups(genes_sorted, n=150):
    s = set()
    count = 0
    for gene in genes_sorted[::-1]:
        
        if count >= n:
            break
        elif gene in s:
            continue

        else:
            try:
                x, y = string_enrich([gene])
            except:
                continue
                print(gene, 'not found')

            s.add(gene)
            count+=1
    return s


def get_genes(small_frac_paths, long_frac_paths, clash_paths, min_hyb_counts=30, min_mir_counts=20, min_mrna_tpm=1, vir='ebv', utr_only=False, top_n_genes=150):
    
    df = integrate(small_frac_paths, long_frac_paths, clash_paths, min_hyb_counts=min_hyb_counts, min_mir_counts=min_mir_counts, min_mrna_tpm=min_mrna_tpm, virus=vir)
    print(df)
    if utr_only is True:
        df = df[df['mrna_feature_start'] == "3'UTR"]
    virus_df = df.loc[set([i for i in df.index if vir in i])].sort_values('hyb')
    host_df = df.loc[set([i for i in df.index if vir not in i])].sort_values('hyb')
    virus_genes = get_n_without_dups(virus_df['mrna'], n=top_n_genes)
    host_genes = get_n_without_dups(host_df['mrna'], n=top_n_genes)
    return virus_genes, host_genes


def neighbor_genes(genes, species=9606):
    
    d = {}
    for gene in genes:
        try:
            x, y = string_enrich([gene])
            d[gene] = y

        except:
            print(gene, 'not found')
    return d


def assemble_pathways(gene_dict, pathway):

    d = {}
    for gene in gene_dict.keys():
        try:
            x = get_enrich(gene_dict[gene], pathway)
            names = list(x['Term name'])
            d[gene] = names
        except:
            print(gene)    

    return d


def count_pathways(assembled_pathway_dict):

    pathways = []
    for key in assembled_pathway_dict:
        for pathway in assembled_pathway_dict[key]:
            pathways.append(pathway)
    
    return Counter(pathways).most_common()


def terms(dic):
    l = []
    for i in dic:
        for j in dic[i]:
            l.append(j)
    l = set(l)
    ge = get_enrich(l, 'Reactome_2016')
    ge['Term name'] = ge['Term name'].map(lambda x:x.split(' Homo s')[0])
    ge = ge[ge['Term name'].isin([i for i in ddd if 5 < len(ddd[i]) < 100])]
    return ge


def save_csv_pathways_targets(pathways_dict, prefix):
    new = defaultdict(list)
    for i in pathways_dict:
        val = pathways_dict[i]
        for path in val:
            new[path].append(i)

    new = {i:', '.join(j) for i,j in new.items()}
    with open(f'{prefix}.pathway_gene_map.tsv', 'w') as outfile:
        outfile.write('Pathway\tSeed gene count\tSeed gene(s)\n')
        for entry in new:
            outfile.write(f'{entry}\t{len(new[entry].split(","))}\t{new[entry]}\n')

def get_pathway_set(pathway_counts, min_counts):
    return set([i[0] for i in pathway_counts if i[1] >= min_counts])


#pathways = 'KEGG_2019_Human'#'GO_Biological_Process_2018' #'GO_Molecular_Function_2018'#'Reactome_2016' #'GO_Cellular_Component_2018' #'KEGG_2019_Human' #'GO_Biological_Process_2018' #GO_Molecular_Function_2018' #

mhv68_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/smallfrac/CL100128477_L1_WHMOUkmcSAAETAASE-545_1.fq.counts.tsv']
mhv68_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/mhv68/longfrac/HE2-1.genes.tsv']
mhv68_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*shuffled_dg.tsv')
mhv68_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/mhv68/clash/*SRR839524[5-7]*annotated')

kshv_small_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/smallfrac/541_WT-TIVE-1.small_fraction.counts.tsv']
kshv_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/kshv/longfrac/TIVE-1.genes.tsv']
kshv_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*shuffled_dg.tsv')
kshv_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/kshv/clash/*annotated')

akata_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*shuffled_dg.tsv')
akata_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/smallfrac/*')
akata_long_frac_paths = ['/Users/nate/Projects/EBV_interactome/akata/longfrac/Akata.genes.tsv']
akata_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/akata/clash/*annotated')

snu_small_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/smallfrac/*')
snu_long_frac_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/longfrac/*')
snu_shuffled_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*shuffled_dg.tsv')
snu_mrna_clash_paths = glob.glob('/Users/nate/Projects/EBV_interactome/snu719/clash/*annotated')

utr_only = False
min_hyb_counts = 30
min_mir_cpm = 20
min_mrna_tpm = 5
top_n = 20 # number of "seed genes"
pathway = 'Reactome_2016'

e, eh = get_genes(snu_small_frac_paths, snu_long_frac_paths, snu_mrna_clash_paths, 100, min_mir_cpm, min_mrna_tpm, 'ebv', True, top_n)
m, mh = get_genes(mhv68_small_frac_paths, mhv68_long_frac_paths, mhv68_mrna_clash_paths, min_hyb_counts, min_mir_cpm, min_mrna_tpm, 'mghv', True, top_n)
k, kh = get_genes(kshv_small_frac_paths, kshv_long_frac_paths, kshv_mrna_clash_paths, min_hyb_counts, min_mir_cpm, min_mrna_tpm, 'kshv', True, top_n)

m = set([i.upper() for i in m])
mh = set([i.upper() for i in mh]) 

ebv = neighbor_genes(e)
mhv = neighbor_genes(m)
kshv = neighbor_genes(k)
ebvh = neighbor_genes(eh)
mhvh = neighbor_genes(mh)
kshvh = neighbor_genes(kh)

ebv_pathways = assemble_pathways(ebv, pathway)
mhv_pathways = assemble_pathways(mhv, pathway)
kshv_pathways = assemble_pathways(kshv, pathway)
ebvh_pathways = assemble_pathways(ebvh, pathway)
mhvh_pathways = assemble_pathways(mhvh, pathway)
kshvh_pathways = assemble_pathways(kshvh, pathway)

cc = count_pathways(ebv_pathways)
ccm = count_pathways(mhv_pathways)
cck = count_pathways(kshv_pathways)
cch = count_pathways(ebvh_pathways)
ccmh = count_pathways(mhvh_pathways)
cckh = count_pathways(kshvh_pathways)

min_pathway_counts = 2
ccs = get_pathway_set(cc, min_pathway_counts)
ccks = get_pathway_set(cck, min_pathway_counts)
ccms = get_pathway_set(ccm, min_pathway_counts)
cchs = get_pathway_set(cch, min_pathway_counts) 
cckhs = get_pathway_set(cckh, min_pathway_counts) 
ccmhs = get_pathway_set(ccmh, min_pathway_counts) 

fig = plt.figure()
venn3([ccs, ccks, ccms])
venn3_circles([ccs, ccks, ccms])

fig2 = plt.figure()
venn3([cchs, cckhs, ccmhs])
venn3_circles([cchs, cckhs, ccmhs])

save_csv_pathways_targets(ebv_pathways, 'ebv')
save_csv_pathways_targets(mhv_pathways, 'mhv68')
save_csv_pathways_targets(kshv_pathways, 'kshv')
save_csv_pathways_targets(ebvh_pathways, 'snuhuman')
save_csv_pathways_targets(mhvh_pathways, 'he21mouse')
save_csv_pathways_targets(kshvh_pathways, 'tivehuman')


with open(f'ebv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in e:
        outfile.write(i + '\n')
with open(f'mhv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in m:
        outfile.write(i + '\n')
with open(f'kshv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in k:
        outfile.write(i + '\n')
    
with open(f'kshv_AND_ebv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in (e & k):
        outfile.write(i + '\n')

with open(f'kshv_AND_mhv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in (k & m):
        outfile.write(i + '\n')

with open(f'ebv_AND_mhv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in (e & m):
        outfile.write(i + '\n')


with open(f'ebv_AND_mhv_and_kshv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in (e & m & k):
        outfile.write(i + '\n')


with open(f'snuhuman_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in eh:
        outfile.write(i + '\n')
with open(f'he21mouse_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in mh:
        outfile.write(i + '\n')
with open(f'tivehuman_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in kh:
        outfile.write(i + '\n')

with open(f'human_AND_ebv_venn_{top_n}_{pathway}_genes.tsv','w') as outfile:
    for i in (e & eh):
        outfile.write(i + '\n')


# Remove duplicate pathways (i.e. TP53 degredation and expression vs TP53 degradation. Leave more specific of two intact). Also remove uninformative/broad pathways (i.e. Infectious disease, or gene expression)

to_remove = ['Gene Expression Homo sapiens R-HSA-74160', 'Infectious disease Homo sapiens R-HSA-5663205','Disease Homo sapiens R-HSA-1643685' , 'Influenza Life Cycle Homo sapiens R-HSA-168255', 
'Regulation of TP53 Expression and Degradation Homo sapiens R-HSA-6806003' ,
'Nonsense-Mediated Decay (NMD) Homo sapiens R-HSA-927802', 'Regulation of Hypoxia-inducible Factor (HIF) by oxygen Homo sapiens R-HSA-1234174',
'Cellular response to hypoxia Homo sapiens R-HSA-2262749', 'Toll Like Receptor 3 (TLR3) Cascade Homo sapiens R-HSA-168164', 
'Activated TLR4 signalling Homo sapiens R-HSA-166054', 'Toll-Like Receptors Cascades Homo sapiens R-HSA-168898', 'G1 Phase Homo sapiens R-HSA-69236',
'Cyclin D associated events in G1 Homo sapiens R-HSA-69231', 'Immune System Homo sapiens R-HSA-168256', 'Toll Like Receptor 4 (TLR4) Cascade Homo sapiens R-HSA-166016','Late Phase of HIV Life Cycle Homo sapiens R-HSA-162599']

to_remove2 = ['Disease Homo sapiens R-HSA-1643685', 'Cell Cycle Homo sapiens R-HSA-1640170', 
'G1 Phase Homo sapiens R-HSA-69236', 'Infectious disease Homo sapiens R-HSA-5663205','Cell Cycle, Mitotic Homo sapiens R-HSA-69278', 'Chromatin organization Homo sapiens R-HSA-4839726','Metabolism Homo sapiens R-HSA-1430728','rRNA processing Homo sapiens R-HSA-72312','Metabolism of proteins Homo sapiens R-HSA-392499']

cc_plot = [i for i in cc[:35] if i[0] not in to_remove][:10]
cch_plot = [i for i in cch[:25] if i[0] not in to_remove2][:10]

fig = plt.figure()
ax = plt.subplot()
plt.barh(range(len(cc_plot)), [i[1] for i in cc_plot])
ax.set_xlim([0,6.5])

plt.savefig('pathways_ebv_mirs.svg')
plt.close()
with open('yaxis_pathways_ebv_mirs.txt', 'w') as outfile:
    for i in cc_plot:
        outfile.write(i[0] + '\n')

fig = plt.figure()
ax = plt.subplot()
plt.barh(range(len(cch_plot)), [i[1] for i in cch_plot])
ax.set_xlim([0,6.5])
plt.savefig('pathways_snu_cell_mirs.svg')
with open('yaxis_pathways_snu_cell_mirs.txt', 'w') as outfile:
    for i in cch_plot:
        outfile.write(i[0] + '\n')


