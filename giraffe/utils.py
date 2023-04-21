import matplotlib.pyplot as plt
from netZooPy.panda import Panda
import numpy as np
import pandas as pd
import patsy
import plotly.express as px
import plotly.graph_objects as go
from sklearn.metrics import roc_curve, auc
import sys

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