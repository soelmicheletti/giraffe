import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc


def evaluate_regulation_auroc(ground_truth, pred) :
    ground_truth[ground_truth != 0] = 1
    fpr, tpr, _ = roc_curve(ground_truth.flatten(), np.abs(pred).flatten())
    roc_aucCont1 = auc(fpr, tpr)
    return roc_aucCont1

def tansform_ppi_for_panda(ppi, path):
    ppi[ppi != 0] = 1
    ppi_panda = pd.DataFrame(0, index=[i for i in range(int(np.sum(np.sum(ppi))))], columns=["from", "to", "score"])
    idx = 0
    for i in range(ppi.shape[0]):
        for j in range(ppi.shape[1]):
            if ppi[j][i] != 0:
                ppi_panda["from"][idx] = i
                ppi_panda["to"][idx] = j
                ppi_panda["score"][idx] = 1
                idx += 1
    ppi_panda.rename(columns={'from': 0, 'to': 1, 'score': 2}, inplace=True)
    ppi_panda.to_csv(path, sep="\t", header=None, index=None)