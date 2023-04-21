import matplotlib.pyplot as plt
from netZooPy.panda import Panda
import numpy as np
import pandas as pd
import patsy
import plotly.express as px
import plotly.graph_objects as go
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.metrics import roc_curve, auc
from sklearn.tree import BaseDecisionTree
import sys

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

def plot_auroc(ground_truth, pred, model_names=None, title="AUROC", approx = True):
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
            if approx:
                auc = round(auc, 3)
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

def limma(pheno, exprs, covariate_formula, design_formula='1', rcond=1e-8):
    """
    Applies limma to correct, for instance, batch effects.
    :param pheno: dataframe with additional variables.
    :param exprs: expression matrix.
    :param covariate_formula: variables to be adjusted.
    :param design_formula: how to adjust for variables.
    :param rcond: tolerance for least square estimator.
    :return: corrected dataset.
    """
    design_matrix = patsy.dmatrix(design_formula, pheno)

    design_matrix = design_matrix[:, 1:]
    rowsum = design_matrix.sum(axis=1) - 1
    design_matrix = (design_matrix.T + rowsum).T

    covariate_matrix = patsy.dmatrix(covariate_formula, pheno)
    design_batch = np.hstack((covariate_matrix, design_matrix))
    coefficients, res, rank, s = np.linalg.lstsq(design_batch, exprs.T, rcond=rcond)
    beta = coefficients[-design_matrix.shape[1]:]
    return exprs - design_matrix.dot(beta).T

def get_sign_accuracy(R, R_hat):
    """
    :param R: true regulatory network.
    :param R_hat: inferred regulatory network.
    :return: sign accuracy, considering hits if corresponding entries in R and R_hat have the same sign.
    """
    hit = 0
    tot = 0
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            if R[i, j] > 0:
                tot += 1
                if R_hat[i, j] > 0:
                    hit += 1
            if R[i, j] < 0:
                tot += 1
                if R_hat[i, j] < 0:
                    hit += 1
    return hit / tot

from sklearn.tree import BaseDecisionTree
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from numpy import *
import time
from operator import itemgetter
from multiprocessing import Pool

def compute_feature_importances(estimator):
    """
    Utility function for GENIE3, code from https://github.com/vahuynh/GENIE3/tree/master/GENIE3_python.
    """
    if isinstance(estimator, BaseDecisionTree):
        return estimator.tree_.compute_feature_importances(normalize=False)
    else:
        importances = [e.tree_.compute_feature_importances(normalize=False)
                       for e in estimator.estimators_]
        importances = array(importances)
        return sum(importances, axis=0) / len(estimator)


def get_link_list(VIM, gene_names=None, regulators='all', maxcount='all', file_name=None):
    """ Code from https://github.com/vahuynh/GENIE3/tree/master/GENIE3_python.
    Gets the ranked list of (directed) regulatory links.
    Parameters
    ----------
    VIM: numpy array
        Array as returned by the function GENIE3(), in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene.

    gene_names: list of strings, optional
        List of length p, where p is the number of rows/columns in VIM, containing the names of the genes. The i-th item of gene_names must correspond to the i-th row/column of VIM. When the gene names are not provided, the i-th gene is named Gi.
        default: None
    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names), and the returned list contains only edges directed from the candidate regulators. When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
    maxcount: 'all' or positive integer, optional
        Writes only the first maxcount regulatory links of the ranked list. When maxcount is set to 'all', all the regulatory links are written.
        default: 'all'
    file_name: string, optional
        Writes the ranked list of regulatory links to the file file_name.
        default: None
    Returns
    -------
    The list of regulatory links, ordered according to the edge score. Auto-regulations do not appear in the list. Regulatory links with a score equal to zero are randomly permuted. In the ranked list of edges, each line has format:
        regulator   target gene     score of edge
    """
    # Check input arguments
    if not isinstance(VIM, ndarray):
        raise ValueError('VIM must be a square array')
    elif VIM.shape[0] != VIM.shape[1]:
        raise ValueError('VIM must be a square array')

    ngenes = VIM.shape[0]

    if gene_names is not None:
        if not isinstance(gene_names, (list, tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError(
                'input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')

    if regulators != 'all':
        if not isinstance(regulators, (list, tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')

    if maxcount != 'all' and not isinstance(maxcount, int):
        raise ValueError('input argument maxcount must be "all" or a positive integer')

    if file_name is not None and not isinstance(file_name, str):
        raise ValueError('input argument file_name must be a string')

    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]

    # Get the non-ranked list of regulatory links
    vInter = [(i, j, score) for (i, j), score in ndenumerate(VIM) if i in input_idx and i != j]

    # Rank the list according to the weights of the edges
    vInter_sort = sorted(vInter, key=itemgetter(2), reverse=True)
    nInter = len(vInter_sort)

    # Random permutation of edges with score equal to 0
    flag = 1
    i = 0
    while flag and i < nInter:
        (TF_idx, target_idx, score) = vInter_sort[i]
        if score == 0:
            flag = 0
        else:
            i += 1

    if not flag:
        items_perm = vInter_sort[i:]
        items_perm = random.permutation(items_perm)
        vInter_sort[i:] = items_perm

    # Write the ranked list of edges
    nToWrite = nInter
    if isinstance(maxcount, int) and maxcount >= 0 and maxcount < nInter:
        nToWrite = maxcount

    if file_name:

        outfile = open(file_name, 'w')

        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('%s\t%s\t%.6f\n' % (gene_names[TF_idx], gene_names[target_idx], score))
        else:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('G%d\tG%d\t%.6f\n' % (TF_idx + 1, target_idx + 1, score))

        outfile.close()

    else:

        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('%s\t%s\t%.6f' % (gene_names[TF_idx], gene_names[target_idx], score))
        else:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('G%d\tG%d\t%.6f' % (TF_idx + 1, target_idx + 1, score))


def GENIE3(expr_data, gene_names=None, regulators='all', tree_method='RF', K='sqrt', ntrees=1000, nthreads=1):
    ''' Code from https://github.com/vahuynh/GENIE3/tree/master/GENIE3_python
    Computation of tree-based scores for all putative regulatory links.
    Parameters
    ----------
    expr_data: numpy array
        Array containing gene expression values. Each row corresponds to a condition and each column corresponds to a gene.
    gene_names: list of strings, optional
        List of length p, where p is the number of columns in expr_data, containing the names of the genes. The i-th item of gene_names must correspond to the i-th column of expr_data.
        default: None
    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names). When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
    tree-method: 'RF' or 'ET', optional
        Specifies which tree-based procedure is used: either Random Forest ('RF') or Extra-Trees ('ET')
        default: 'RF'
    K: 'sqrt', 'all' or a positive integer, optional
        Specifies the number of selected attributes at each node of one tree: either the square root of the number of candidate regulators ('sqrt'), the total number of candidate regulators ('all'), or any positive integer.
        default: 'sqrt'
    ntrees: positive integer, optional
        Specifies the number of trees grown in an ensemble.
        default: 1000
    nthreads: positive integer, optional
        Number of threads used for parallel computing
        default: 1
    Returns
    -------
    An array in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene. All diagonal elements are set to zero (auto-regulations are not considered). When a list of candidate regulators is provided, the scores of all the edges directed from a gene that is not a candidate regulator are set to zero.
    '''
    def compute_feature_importances(estimator):
        if isinstance(estimator, BaseDecisionTree):
            return estimator.tree_.compute_feature_importances(normalize=False)
        else:
            importances = [e.tree_.compute_feature_importances(normalize=False)
                           for e in estimator.estimators_]
            importances = array(importances)
            return sum(importances, axis=0) / len(estimator)

    def wr_GENIE3_single(args):
        return ([args[1], GENIE3_single(args[0], args[1], args[2], args[3], args[4], args[5])])

    def GENIE3_single(expr_data, output_idx, input_idx, tree_method, K, ntrees):

        ngenes = expr_data.shape[1]

        # Expression of target gene
        output = expr_data[:, output_idx]

        # Normalize output data
        output = output / std(output)

        # Remove target gene from candidate regulators
        input_idx = input_idx[:]
        if output_idx in input_idx:
            input_idx.remove(output_idx)

        expr_data_input = expr_data[:, input_idx]

        # Parameter K of the tree-based method
        if (K == 'all') or (isinstance(K, int) and K >= len(input_idx)):
            max_features = "auto"
        else:
            max_features = K

        if tree_method == 'RF':
            treeEstimator = RandomForestRegressor(n_estimators=ntrees, max_features=max_features)
        elif tree_method == 'ET':
            treeEstimator = ExtraTreesRegressor(n_estimators=ntrees, max_features=max_features)

        # Learn ensemble of trees
        treeEstimator.fit(expr_data_input, output)

        # Compute importance scores
        feature_importances = compute_feature_importances(treeEstimator)
        vi = zeros(ngenes)
        vi[input_idx] = feature_importances

        return vi

    time_start = time.time()
    # Check input arguments
    if not isinstance(expr_data, ndarray):
        raise ValueError(
            'expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
    ngenes = expr_data.shape[1]
    if gene_names is not None:
        if not isinstance(gene_names, (list, tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError(
                'input argument gene_names must be a list of length p, where p is the number of columns/genes in the expr_data')
    if regulators != 'all':
        if not isinstance(regulators, (list, tuple)):
            raise ValueError('input argument regulators must be a list of gene names')
        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('the genes must contain at least one candidate regulator')
    if tree_method != 'RF' and tree_method != 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')
    if K != 'sqrt' and K != 'all' and not isinstance(K, int):
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
    if isinstance(K, int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
    if not isinstance(ntrees, int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')
    if not isinstance(nthreads, int):
        raise ValueError('input argument nthreads must be a stricly positive integer')
    elif nthreads <= 0:
        raise ValueError('input argument nthreads must be a stricly positive integer')
    print('Tree method: ' + str(tree_method))
    print('K: ' + str(K))
    print('Number of trees: ' + str(ntrees))
    print('\n')
    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = list(range(ngenes))
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]
    # Learn an ensemble of trees for each target gene, and compute scores for candidate regulators
    VIM = zeros((ngenes, ngenes))
    if nthreads > 1:
        print('running jobs on %d threads' % nthreads)
        input_data = list()
        for i in range(ngenes):
            input_data.append([expr_data, i, input_idx, tree_method, K, ntrees])
        pool = Pool(nthreads)
        alloutput = pool.map(wr_GENIE3_single, input_data)
        for (i, vi) in alloutput:
            VIM[i, :] = vi
    else:
        print('running single threaded jobs')
        for i in range(ngenes):
            print('Gene %d/%d...' % (i + 1, ngenes))
            vi = GENIE3_single(expr_data, i, input_idx, tree_method, K, ntrees)
            VIM[i, :] = vi
    VIM = transpose(VIM)
    time_end = time.time()
    print("Elapsed time: %.2f seconds" % (time_end - time_start))
    return VIM