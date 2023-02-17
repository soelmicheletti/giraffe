import pandas as pd
import requests

def get_neighbors(P) :
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])
    my_genes = [P]
    params = {
        "identifiers" : "%0d".join(my_genes),
        "species" : 9606,
        "network_type" : "physical",
        "caller_identity" : "www.awesome_app.org"
    }
    response = requests.post(request_url, data=params)
    res = []
    for line in response.text.strip().split("\n"):

        l = line.strip().split("\t")
        if len(l) < 3 :
            return []
        p1, p2 = l[2], l[3]
        experimental_score = float(l[10])
        if experimental_score > 0.4 and p1 == P:
            res += [p2]
    return list(set(res))

def generate_data():
    motif = pd.read_csv("data/liver/raw/motif.txt", sep="\t", header=None)
    motif_filter = motif[motif[2] == 1]
    expression = pd.read_csv("data/liver/raw/expression.gct", skiprows=2, sep='\t', index_col='Name',
                             usecols=lambda x: x != 'id')
    new_index = [idx[0:15] for idx in expression.index]
    expression.index = new_index
    motif_matrix = pd.DataFrame(0, index=expression.index, columns=list(set(motif[0])))
    for i in range(motif_filter.shape[0]):
        if i % 10000 == 0:
            print(i)
        tf = motif_filter.iloc[i][0]
        g = motif_filter.iloc[i][1]
        if g in motif_matrix.index and tf in motif_matrix.columns:
            motif_matrix[tf][g] = 1
    to_drop = []
    for i in range(len(expression.index)):
        if sum(expression.iloc[i][expression.columns[1:]]) == 0 or sum(motif_matrix.iloc[i]) == 0:
            to_drop.append(i)
    keys = [expression.index[i] for i in to_drop]
    expression = expression.drop(keys)
    motif_matrix = motif_matrix.drop(keys)
    new_index = expression['Description']
    expression.index = new_index
    expression = expression[expression.columns[1:]]
    motif_matrix.index = new_index
    motif_matrix.to_csv("data/liver/motif.txt")
    expression.to_csv("data/liver/expression.txt")
    motif = pd.read_csv("data/liver/motif.txt", index_col = 0)
    ppi = pd.DataFrame(0, index=motif.columns, columns=motif.columns)
    for i in ppi.columns:
        res = get_neighbors(i)
        for r in res:
            if r not in ppi.index:
                continue
            ppi[i][r] = 1
            ppi[r][i] = 1
    ppi.to_csv("data/liver/ppi.csv")
