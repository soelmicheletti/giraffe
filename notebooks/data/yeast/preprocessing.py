import giraffe
import pandas as pd
import requests

def alias(gene):
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"

    params = {
        "identifiers" : "\r".join([gene]),
        "species" : 4932,
        "limit" : 1,
        "echo_query" : 1,
        "caller_identity" : "www.awesome_app.org" # your app name
    }
    request_url = "/".join([string_api_url, output_format, method])
    results = requests.post(request_url, data=params)
    alias = []
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        input_identifier, string_identifier = l[0], l[2]
        alias += [string_identifier.lstrip("4932.")]
    return set(alias)


def get_neighbors(P):
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    request_url = "/".join([string_api_url, output_format, method])
    my_genes = [P]
    params = {
        "identifiers": "%0d".join(my_genes),
        "species": 4932,
        "network_type": "physical",
        "caller_identity": "www.awesome_app.org"
    }

    response = requests.post(request_url, data=params)
    res = []
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        if len(l) < 3:
            return []
        p1, p2 = l[2], l[3]
        experimental_score = float(l[10])
        if experimental_score > 0.4:
            if P not in alias(p1) and P in alias(p2):
                res += [p1]
            if P not in alias(p2) and P in alias(p1):
                res += [p2]
    return list(set(res))

def query_ppi(TF):
    PPI = pd.DataFrame(0, index=TF, columns=TF)
    cnt = 0
    for protein in TF:
        res = get_neighbors(protein)
        cnt += 1
        for r in res:
            target = list(alias(r))
            for l in target:
                if l in PPI.index:
                    PPI[target][protein] = 1
                    PPI[protein][target] = 1
    PPI.to_csv("data/yeast/ppi.txt")

def save_motif_matrix(motif):
    giraffe.utils.transform_motif_to_matrix(motif, "data/yeast/motif_matrix.txt")

def generate_data():
    motif = pd.read_csv("data/yeast/raw/YeastCCData_Motif.txt", sep='\t', header=None)
    #  We update some names from the original file
    def replace_old_names(motif_edge_list, old_name, new_name):
        for i in range(motif_edge_list.shape[0]):
            if motif_edge_list[0][i] == old_name:
                motif_edge_list[0][i] = new_name
    replace_old_names(motif, "SIG1", "YER068W")
    replace_old_names(motif, "RCS1", "YGL071W")
    replace_old_names(motif, "RLR1", "YNL139C")
    motif.to_csv("data/yeast/motif_edge_list.txt")
    save_motif_matrix(motif)
    TF = set(motif[0])
    query_ppi(TF)