import pandas as pd
import requests
from scipy.io import mmread

def binarySearch(data, val):
    """
    Performs binary search to find indices of a gene promoter region in the chromatine dataset
    :param data: dataframe containing ATAC-seq data in every region of the genome
    :param val: region to be searched
    :return: index of region in data
    """
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
        if mid + 1 < len(data) and data[mid] < val and data[mid + 1] >= val:
            return mid
        if abs(data[mid] - val) < abs(data[best_ind] - val):
            best_ind = mid
    return best_ind

def get_motif(genes):
    """
    Transform list of motif edges to an adjacency matrix
    :param genes: dataframe containing gene names
    :return: motif matrix
    """
    motif = pd.read_csv("data/adult_cerebrail_cortex/raw/mm10_reseq_-750_250_1e-4.txt", sep='\t', header=None)
    motif_matrix = pd.DataFrame(0, index=genes[1], columns=set(motif[0]))
    motif = motif[motif[2] == 1]
    for i in range(motif.shape[0]):
        tf = motif[0][i]
        g = motif[1][i]
        if g in set(genes[1]):
            motif_matrix.loc[g, tf] = 1

def get_ppi(motif):
    """
    Extract the ppi network from
    :param motif: dataframe of dimension G x TF containing the prior to the regulatory network
    :return: ppi network as dataframe of dimension TF x TF
    """
    name_to_entrez = pd.read_csv("data/adult_cerebrail_cortex/raw/mippie_proteins_v1_0.tsv", sep='\t')
    ppi_list = pd.read_csv("data/adult_cerebrail_cortex/raw/mippie_ppi_v1_0.tsv", sep='\t')
    id_to_name = pd.read_csv("data/adult_cerebrail_cortex/raw/Mus_musculus_motifinfo.txt", sep='\t')
    ppi = pd.DataFrame(0, index = motif.columns, columns = motif.columns)
    for i in ppi.columns:
        for j in ppi.columns:
            if i + "_1.02" in set(id_to_name['MotifID']) and j + "_1.02" in set(id_to_name['MotifID']):
                tf1 = id_to_name[id_to_name['MotifID'] == i + "_1.02"]['TF_Info'].values[0]
                tf2 = id_to_name[id_to_name['MotifID'] == j + "_1.02"]['TF_Info'].values[0]
                if tf1 in set(name_to_entrez['official_symbol']) and tf2 in set(name_to_entrez['official_symbol']):
                    entrez1 = name_to_entrez[name_to_entrez['official_symbol'] == tf1]['entrez'].values[0]
                    entrez2 = name_to_entrez[name_to_entrez['official_symbol'] == tf2]['entrez'].values[0]
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
    return ppi

def save_gene_accessbility(genes, rows, chromatine):
    """
    Saves number of counts per cell for every gene
    :param genes: list of gene names
    :param rows: list of gene regions in the chromatine dataset
    :param chromatine: counts in the dataset
    """
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

    cells = pd.read_csv('data/adult_cerebrail_cortex/raw/Adult_CTX_DNA/barcodes.tsv', header=None)
    gene_locations = pd.read_csv("data/adult_cerebrail_cortex/raw/mart_export.txt", sep="\t", index_col=0)
    access = []
    to_drop = []

    for i in range(len(genes)):
        if i % 500 == 0:
            print(i)
        if i % 1000 == 0:
            idx = [j for j in range(i)]
            idx = list(set(idx) - set(to_drop))
            pd.DataFrame(access, index=genes[0][idx]).to_csv("data/adult_cerebrail_cortex/accessibility.txt")
        if genes[0][i] not in set(gene_locations.index):
            to_drop.append(i)
            continue
        gene_info = gene_locations.loc[genes[0][i]]
        strand = gene_info['Strand']
        chromosome = gene_info['Chromosome/scaffold name']
        if isinstance(chromosome, pd.core.series.Series):
            chromosome = chromosome.values[0]
        start_list = list(rows_processed[rows_processed['chr'] == 'chr' + str(chromosome)]['start'])
        if (isinstance(gene_locations.loc[genes[0][i]]['Strand'], np.int64) and strand == 1) or (
                not isinstance(gene_locations.loc[genes[0][i]]['Strand'], np.int64) and strand.values[0] == 1):
            start = gene_info['Gene start (bp)']
            if isinstance(start, pd.core.series.Series):
                start = start.values[0]
            start_index = binarySearch(start_list, start - 500)
            end_index = binarySearch(start_list, start)
            res = [0 for i in range(len(cells))]
            for idx in range(start_index, end_index + 1):
                if idx >= len(start_list):
                    continue
                key = "chr" + str(chromosome) + ":" + str(start_list[idx]) + "-" + str(start_list[idx] + 1000)
                offset = list(rows[0]).index(key)
                if idx != start_index and idx != end_index:
                    res += chromatine.getrow(offset).toarray()[0]
                elif idx == start_index:
                    res += chromatine.getrow(offset).toarray()[0] * (1000 - ((start - 500) % 1000)) / 1000
                elif idx == end_index:
                    res += chromatine.getrow(offset).toarray()[0] * (start % 1000) / 1000
            access.append(res)
        else:
            end = gene_info['Gene end (bp)']
            if isinstance(end, pd.core.series.Series):
                end = end.values[0]
            start_index = binarySearch(start_list, end)
            end_index = binarySearch(start_list, end + 500)
            res = [0 for i in range(len(cells))]
            for idx in range(start_index, end_index + 1):
                if idx >= len(start_list):
                    continue
                key = "chr" + str(chromosome) + ":" + str(start_list[idx]) + "-" + str(start_list[idx] + 1000)
                offset = list(rows[0]).index(key)
                if idx != start_index and idx != end_index:
                    res += chromatine.getrow(offset).toarray()[0]
                elif idx == start_index:
                    res += chromatine.getrow(offset).toarray()[0] * (1000 - (end % 1000)) / 1000
                elif idx == end_index:
                    res += chromatine.getrow(offset).toarray()[0] * ((end + 500) % 1000) / 1000
            access.append(res)
    idx = [j for j in range(len(genes))]
    idx = list(set(idx) - set(to_drop))
    pd.DataFrame(access, index=genes[0][idx]).to_csv("data/adult_cerebrail_cortex/accessibility.txt")

def generate_data():
    genes = pd.read_csv('data/adult_cerebrail_cortex/raw/Adult_CTX_RNA/genes.tsv', header=None, sep='\t')

    chromatine = mmread('data/adult_cerebrail_cortex/raw/Adult_CTX_DNA/matrix.mtx')
    rows = pd.read_csv('data/adult_cerebrail_cortex/raw/Adult_CTX_DNA/genes.tsv', header=None, sep='\t')

    save_gene_accessbility(genes, rows, chromatine)

    motif = get_motif(genes)
    motif.to_csv("data/adult_cerebrail_cortex/motif.txt")

    ppi = get_ppi(motif)
    ppi.to_csv("data/adult_cerebrail_cortex/ppi.txt")

