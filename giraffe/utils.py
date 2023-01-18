import matplotlib.pyplot as plt
from netZooPy.panda import Panda
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sklearn.metrics import roc_curve, auc


def check_symmetric(ppi):
    """
    Raises an exception if ppi does not have the expected structure.

    :param ppi: matrix to be checked
    """
    if ppi.shape[0] != ppi.shape[1]:
        raise Exception(
            "PPI must be a squared matrix. "
        )
    if np.diag(ppi).any() != 1:
        raise Exception(
            "PPI matrix must have ones on the diagonal. "
        )
    if np.any(np.abs(ppi - ppi.T) > 1e-8):
        raise Exception(
            "PPI matrix must be symmetrical. "
        )


def compute_panda_locally(expression, motif, ppi, path="tmp"):
    """
    Computes PANDA using the code from netZooPy. Note that for large datasets, this be slow or run out of memory.
    If necessary, consider moving to a cluster and using the command line interface, see
    https://netzoopy.readthedocs.io/en/latest/functions/cli.html [13th Jan. 2023]

    :param expression: numpy array of dimension g x n.
    :param motif: numpy array of dimension g x tf. Prior for the regulation matrix.
    :param ppi: numpy array of timension tf x tf. Represents the adjacency matrix of the (symmetric) interactions between TFs.
    :return:regulation matrix computed by panda. Numpy matrix of dimension g x tf.
    """
    panda_obj = Panda(
        pd.DataFrame(expression),
        transform_motif_for_panda(motif, path + "/motif.txt"),
        transform_ppi_for_panda(ppi, path + "/ppi.txt"),
        save_tmp=False,
        save_memory=False,
        remove_missing=False,
        keep_expression_matrix=False
    )
    panda_obj.save_panda_results(path + "/panda_res.txt")
    return transform_panda_output(path + "/panda_res.txt", ppi.shape[0], expression.shape[0])


def evaluate_regulation_auroc(ground_truth, pred):
    """
    Computes AUROC, TPR, FPR for given predictions of a ground-truth binary array.

    :param ground_truth: binary numpy array containing the ground-truth.
    :param pred: array containing the predicted probability that corresponding entry in ground-truth is one.
    :return: AUROC score, true positive rate (tpf), false positive rate (fpr).
    """
    ground_truth[ground_truth != 0] = 1
    fpr, tpr, _ = roc_curve(ground_truth.flatten(), np.abs(pred).flatten())
    roc_aucCont = auc(fpr, tpr)
    return roc_aucCont, tpr, fpr


def plot_auroc(ground_truth, pred, model_names=None, title="AUROC"):
    """
    Plots AUROC for the various models passed in pred.

    :param ground_truth: binary numpy array containing the ground-truth.
    :param pred: list of arrays containing the predicted probability that corresponding entry in ground-truth is one.
    :param model_names: names of the models used to obtain the predictions. Must have the same length as pred.
    :param title: Title of the plotted AUROC.
    :return: list of AUROC scores
    """
    if model_names is not None:
        if len(pred) != len(model_names):
            raise Exception(
                "Model names must have the same length as the number of provided models in the prediction vector. "
            )
        scores = []
        for i, model_predictions in enumerate(pred):
            auc, tpr, fpr = evaluate_regulation_auroc(ground_truth, model_predictions)
            scores.append(auc)
            if model_names is None:
                plt.plot(fpr, tpr)
            else:
                plt.plot(fpr, tpr, label=model_names[i] + ", AUROC = " + str(auc))

        plt.legend(loc='lower right')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.title(title)
        plt.show()
        return scores


def star_plot(scores, metrics_names, model_names):
    """
    Star plot to compare score of multiple models across multiple metrics.
    Consider mod different models and met different metrics

    :param scores: numpy array of dimension mod x met. Each entry is the corresponding sccore for the given model/ metric.
    :param metrics_names: numpy array of dimension met with a description of the given metric.
    :param model_names: list of length mod with the names of the models to be compared.
    """
    if len(scores) != len(model_names):
        raise Exception(
            "The number of models provided in the score array must be equal to the length ot the list of model names!"
        )
    fig = go.Figure()

    for idx, score in enumerate(scores):
        if len(score) != len(metrics_names):
            raise Exception(
                "Number of provided scores for " + str(model_names[idx]) + " wrong. Please provide " + len(
                    metrics_names) + " scores."
            )
        fig.add_trace(go.Scatterpolar(
            r=score,
            theta=metrics_names,
            fill='toself',
            name=model_names[idx]
        ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True
            )
        ),
        template='plotly_dark',
        showlegend=True
    )
    fig.show()


def transform_motif_for_panda(motif, path):
    """
    Transforms motif matrix of dimension g x tf into a PANDA friendly .txt file, which will be saved in path.

    :param motif: motif: numpy array of dimension g x tf. Prior for the regulation matrix.
    :param path: path to store the PANDA friendly .txt file.
    :return: The path where the .txt file has been successfully stored.
    """
    motif[motif != 0] = 1  # PANDA requires a binary motif
    W = pd.DataFrame(0, index=[i for i in range(int(np.sum(np.sum(motif))))], columns=["source", "target", "weight"])
    idx = 0
    for i in range(motif.shape[0]):
        for j in range(motif.shape[1]):
            if motif[i, j] == 1:
                W["source"][idx] = i
                W["target"][idx] = j
                W["weight"][idx] = 1
                idx += 1
    W.to_csv(path, sep="\t", header=None, index=None)
    return path

def transform_motif_to_matrix(motif, path):
    """
    Transforms pandas dataframe containing edge list of motif into adjacency matrix of size

    :param motif: pandas Dataframe with three columns: tf - target - {0, 1}.
    :param path: path to store the motif in matrix format.
    :return: The path where the .txt file has been successfully stored.
    """
    W = pd.DataFrame(0, index=list(set(motif[0])), columns=list(set(motif[1])))

    for i in range(motif.shape[0]):
        tf = motif[0][i]
        target = motif[1][i]
        if motif[2][i] == 1:
            W[target][tf] = 1
    W.T.to_csv(path)
    return path

def transform_panda_output(path, tf, g):
    """
    Transforms .txt file containing regulation in a numpy matrix of dimension g x tf.

    :param path: path storing the PANDA output.
    :param tf: for convenience.
    :param g: for convenience.
    :return: numpy matrix of size g x tf.
    """
    regulation = pd.read_csv(path, header=None, sep=" ")
    R_panda = pd.DataFrame(np.zeros((g, tf)))
    for i in range(regulation.shape[0]):
        tf = regulation[0][i]
        target = regulation[1][i]
        if tf in R_panda.index and target in R_panda.columns:
            R_panda[target][tf] = regulation[3][i]
    return R_panda.to_numpy()


def transform_ppi_for_panda(ppi, path):
    """
    Transforms PPI matrix of dimension tf x tf into a PANDA friendly .txt file, which will be saved in the provided path.

    :param ppi: numpy matrix of size tf x tf. Must be symmetric.
    :param path: path to store the PANDA friendly .txt file.
    :return: The path where the .txt file has been successfully stored.
    """
    ppi[ppi != 0] = 1
    check_symmetric(ppi)
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
    return path
