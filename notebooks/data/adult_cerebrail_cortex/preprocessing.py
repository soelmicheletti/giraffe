import pandas as pd
import requests
from scipy.io import mmread

def binarySearch(data, val):
    lo, hi = 0, len(data) - 1
    best_ind = lo
    while lo <= hi:
        mid = lo + (hi - lo) // 2
        if data[mid] < val:
            lo = mid + 1
        elif data[mid] > val:
            hi = mid - 1
        else:
            best_ind = mid
            break
        # check if data[mid] is closer to val than data[best_ind]
        if abs(data[mid] - val) < abs(data[best_ind] - val):
            best_ind = mid
    return best_ind

def generate_data():
    genes = pd.read_csv('data/adult_cerebrail_cortex/raw/Adult_CTX_RNA/genes.tsv', header=None, sep='\t')

    chromatine = mmread('data/adult_cerebrail_cortex/raw/Adult_CTX_DNA/matrix.mtx')
    rows = pd.read_csv('data/adult_cerebrail_cortex/raw/Adult_CTX_DNA/genes.tsv', header=None, sep='\t')

    gene_locations = pd.DataFrame(0, index = genes[0], columns = ['gene_name', 'chr', 'start', 'end'])
    gene_locations['gene_name'] = list(genes[1])
    starts = []
    ends = []
    c = []
    for g in genes[0][len(starts):]:
        print(g)
        r = requests.get("https://useast.ensembl.org/Mus_musculus/Gene/Summary?g=" + str(g))
        a = str(r.content)
        b = a.split("r=", 1)[1]
        chromosome = b.split(":", 1)[0]
        start = b.split(":", 1)[1].split("-", 1)[0]
        end = b.split(":", 1)[1].split("-", 1)[1].split("\"", 1)[0]
        if ';' in end:
            end = end.split(";", 1)[0]
        c.append(chromosome)
        starts.append(start)
        ends.append(end)
        gene_locations.loc[g, 'chr'] = chromosome
        gene_locations.loc[g, 'start'] = start
        gene_locations.loc[g, 'end'] = end

    gene_locations.to_csv("data/adult_cerebrail_cortex/gene_locations.csv")

    chr = []
    start = []
    end = []
    for entry in rows[0]:
        chr.append(entry.split(":", 1)[0])
        start.append(int(entry.split(":", 1)[1].split("-", 1)[0]))
        end.append(int(entry.split(":", 1)[1].split("-", 1)[1]))
    rows_processed = pd.DataFrame(0, index=rows.index, columns=["chr", "start", "end"])
    rows_processed['chr'] = chr
    rows_processed['start'] = start
    rows_processed['end'] = end
    access = []
    for i in range(len(genes)):
        if i % 500 == 0:
            print(i)
        gene_info = gene_locations.loc[genes[0][i]]
        chromosome = gene_info['chr']
        start = gene_info['start']
        if not start.isdigit():
            access.append(-1)
            continue
        idx = binarySearch(list(rows_processed[rows_processed['chr'] == 'chr' + str(chromosome)]['start']), int(start))
        access.append(np.sum(chromatine.getrow(idx).toarray()[0]))
    pd.DataFrame(access).to_csv("data/adult_cerebrail_cortex/chromatine.txt")

    motif = pd.read_csv("data/adult_cerebrail_cortex/raw/mm10_reseq_-750_250_1e-4.txt", sep='\t', header=None)
    motif_matrix = pd.DataFrame(0, index=genes[1], columns=set(motif[0]))
    motif = motif[motif[2] == 1]
    for i in range(motif.shape[0]):
        tf = motif[0][i]
        g = motif[1][i]
        if g in set(genes[1]):
            motif_matrix.loc[g, tf] = 1
    motif_matrix.to_csv("data/adult_cerebrail_cortex/motif.txt")
    prot_info = pd.read_csv("data/adult_cerebrail_cortex/raw/mippie_proteins_v1_0.tsv", sep='\t')
    ppi_list = pd.read_csv("data/adult_cerebrail_cortex/raw/mippie_ppi_v1_0.tsv", sep='\t')
    translate = pd.read_csv("data/adult_cerebrail_cortex/raw/Mus_musculus_motifinfo.txt", sep='\t')
    ppi = pd.DataFrame(0, index=motif_matrix.columns, columns=motif_matrix.columns)
    for i in ppi.columns:
        for j in ppi.columns:
            if i + "_1.02" in set(translate['MotifID']) and j + "_1.02" in set(translate['MotifID']):
                tf1 = translate[translate['MotifID'] == i + "_1.02"]['TF_Info'].values[0]
                tf2 = translate[translate['MotifID'] == j + "_1.02"]['TF_Info'].values[0]
                if tf1 in set(prot_info['official_symbol']) and tf2 in set(prot_info['official_symbol']):
                    entrez1 = prot_info[prot_info['official_symbol'] == tf1]['entrez'].values[0]
                    entrez2 = prot_info[prot_info['official_symbol'] == tf2]['entrez'].values[0]
                    if entrez1 in set(ppi_list['entrezA']) and entrez2 in set(ppi_list['entrezB']):
                        tmp = ppi_list[ppi_list['entrezA'] == entrez1]
                        val = tmp[tmp['entrezB'] == entrez2]['MIPPIE_score'].values
                        if val >= 0.4:
                            ppi.loc[i, j] = 1
                            ppi.loc[j, i] = 1
                    if entrez2 in set(ppi_list['entrezA']) and entrez1 in set(ppi_list['entrezB']):
                        tmp = ppi_list[ppi_list['entrezA'] == entrez2]
                        val = tmp[tmp['entrezB'] == entrez1]['MIPPIE_score'].values
                        if val >= 0.4:
                            ppi.loc[i, j] = 1
                            ppi.loc[j, i] = 1
