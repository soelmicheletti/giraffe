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
    expression = pd.read_csv("data/liver/raw/expression.csv", index_col = 0)
    motif_matrix = pd.DataFrame(0, index=expression.index, columns=list(set(motif[0])))
    for i in range(motif_filter.shape[0]):
        if i % 10000 == 0:
            print(i)
        tf = motif_filter.iloc[i][0]
        g = motif_filter.iloc[i][1]
        if g in motif_matrix.index and tf in motif_matrix.columns:
            motif_matrix[tf][g] = 1
    motif_matrix.to_csv("data/liver/motif.csv")
    ppi = pd.DataFrame(0, index=motif_matrix.columns, columns=motif_matrix.columns)
    for i in ppi.columns:
        res = get_neighbors(i)
        for r in res:
            if r not in ppi.index:
                continue
            ppi[i][r] = 1
            ppi[r][i] = 1
    ppi.to_csv("data/liver/ppi_matrix.csv")
    motif = pd.read_csv("data/liver/motif.csv", index_col = 0)
    chip = pd.DataFrame(0, index=motif.index, columns=motif.columns)
    chip_raw = pd.read_csv("data/liver/raw/TF-Target-information.txt", on_bad_lines='skip', sep='\t')
    chip_raw = chip_raw[chip_raw['tissue'] == 'liver']
    translation = pd.read_csv("data/liver/raw/gen_v26_mapping.csv")
    for i in range(chip_raw.shape[0]):
        tf = chip_raw.iloc[i, :]['TF']
        g = chip_raw.iloc[i, :]['target']
        if tf not in chip.columns:
            continue
        translations = list(translation[translation['gene_name'] == g]['gene_id'])
        for t in translations:
            if t[0:15] in chip.index:
                chip.loc[t[0:15], tf] = 1
    chip.to_csv("data/liver/chip.csv")

