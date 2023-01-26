import pandas as pd

def get_preliminary():
    ppi_list = pd.read_csv("data/luad/raw/ppi_list_of_edges.csv", index_col=0)
    tf = list(set(ppi_list['V1']))
    exp = pd.read_csv("data/luad/raw/male_lung.csv", index_col=0)
    genes = list(set(exp.index))
    return tf, genes

def motif_to_matrix(motif, tf, genes):
    motif_list = pd.read_csv(motif[0], index_col = 0)
    motif_list = motif_list[motif_list['V3'] == 1]
    print(motif_list.shape)
    motif_matrix = pd.DataFrame(0, index = genes, columns = tf)
    for i in range(motif_list.shape[0]):
        if i % 10000 == 0:
            print(i)
        source = motif_list.iloc[i]['V1']
        target = motif_list.iloc[i]['V2']
        val = motif_list.iloc[i]['V3']
        if val == 1 and source in tf and target in genes:
            motif_matrix[source][target] = 1
    motif_matrix.to_csv(motif[1])


def ppi_to_matrix():
    ppi_list = pd.read_csv("data/luad/raw/ppi_list_of_edges.csv", index_col=0)
    tf = set(ppi_list['V1'])
    ppi = pd.DataFrame(0, index=list(tf), columns=list(tf))
    for i in range(ppi_list.shape[0]):
        u = ppi_list.iloc[i]['V1']
        v = ppi_list.iloc[i]['V2']
        w = ppi_list.iloc[i]['V3']
        if v not in ppi.index:
            continue
        ppi[u][v] = w
        ppi[v][u] = w
    ppi.to_csv("data/luad/ppi.txt")

def generate_data():
    tf, genes = get_preliminary()
    #ppi_to_matrix()
    #motifs = [("data/luad/raw/motif_list_of_edges_female.csv", "data/luad/motif_female.txt"),
    #          ("data/luad/raw/motif_list_of_edges_male.csv", "data/luad/motif_male.txt")
    #           ]
    #for motif in motifs:
    #    motif_to_matrix(motif, tf, genes)