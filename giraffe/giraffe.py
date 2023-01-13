import numpy as np
import pandas as pd
import torch
from torch.functional import F
from torch import nn

class Giraffe(object):
    """
            GIRAFFE infers regulation and transcription factor activity using

            gene expression
            TF DNA binding motif
            TF interactions (PPI)

            while (optionally) controlling for covariates of interest (e.g. age).
            There are G genes, TF transcription factors, and n samples (cells, patients, ...)

        	We minimize the following objective:
        	    min f(R, TFA)
        	    with f(R, TFA) = TODO

        	Parameters
            -----------
        	    expression : numpy array or pandas dataframe

        	    motif : numpy array or pandas dataframe
                    TF-gene regulatory network based on TF motifs as a
        	        matrix of size (t,g), g=number of genes, t=number of TFs
        	    PPI : numpy array or pandas dataframe
                    TF-TF protein interaction network as a matrix of size (TF, TF)
        	    iterations : int
                    number of iterations of the algorithm
        	    lr : float
                    the learning rate
                lam:

        	Returns
            ---------
        	    W     : array
                    Predicted TF-gene complete regulatory network as an adjacency matrix of size (t,g).

            References
            -------------
                .. TODO
        """
    def __init__(
            self,
            expression,
            motif,
            ppi,
            adjusting = None,
            iterations = 200,
            lr = 0.00001,
            lam = None
    ):
        self.process_data(
            expression,
            motif,
            ppi,
            adjusting
        )
        self._iterations = iterations
        self._lr = lr
        self._lam = lam
        self._R, self._TFA = self._compute_giraffe()

    def process_data(
            self,
            expression,
            motif,
            ppi,
            adjusting
    ):
        """ Processes data files into data matrices.
               Parameters
               ----------
                   expression_file : str
                       numpy matrix or pandas dataframe containing gene expression data
                   motif_file : str
                       numpy matrix or pandas dataframe containing motif data
                   ppi_file : str
                       Numpy or pandas dataframe containing PPI data.
                       The PPI should be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix.
        """
        self._expression = expression
        self._motif = motif
        self._ppi = ppi

        if isinstance(expression, pd.DataFrame):
            self._expression = expression.to_numpy()
        if isinstance(motif, pd.DataFrame):
            self._motif = motif.to_numpy()
        if isinstance(ppi, pd.DataFrame):
            self._ppi = ppi.to_numpy()

        if not isinstance(self._expression, np.ndarray):
            raise Exception(
                "Error in processing expression data. Please provide a numpy array or a pandas dataframe. "
            )
        if not isinstance(self._motif, np.ndarray):
            raise Exception(
                "Error in processing motif data. Please provide a numpy array or a pandas dataframe. "
            )
        if not isinstance(self._ppi, np.ndarray):
            raise Exception(
                "Error in processing PPI data. Please provide a numpy array or a pandas dataframe. "
            )

        self._adjusting = adjusting
        if isinstance(adjusting, pd.DataFrame):
                self._adjusting = adjusting.to_numpy()
        if self._adjusting is not None and not isinstance(self._adjusting, np.ndarray):
            raise Exception(
                    "Error in processing adjusting covariates. Please provide a numpy array or a pandas dataframe. "
            )
        if self._adjusting is not None :
            self._adjusting = torch.Tensor(self._adjusting)
            self._adjusting /= torch.norm(torch.Tensor(self._adjusting)) # Normalize

        # Normalize motif and PPI
        self._motif = motif / np.sqrt(np.trace(motif.T @ motif))
        self._ppi = self._ppi / np.trace(self._ppi)
        # Compute co-expression matrix
        self._C = np.corrcoef(self._expression)
        self._C /= np.trace(self._C)

    def _compute_giraffe(self):
        """
        Giraffe optimization
        :return: R : matrix G x TF, partial effects between tanscription factor and gene
                 TFA : matrix TF x n, transcrption factor activity
        """
        giraffe_model = Model(np.random.random((self._motif.shape[1], self._expression.shape[1])), self._motif)
        optim = torch.optim.Adam(giraffe_model.parameters(), lr = self._lr) # We run Adam to optimize f(R, TFA)

        for i in range(self._iterations):
            pred = giraffe_model(
                torch.Tensor(self._expression),
                torch.Tensor(self._ppi),
                torch.Tensor(self._C),
                self._lam,
                self._adjusting
            ) # Compute f(R, TFA)
            loss = F.mse_loss(pred, torch.norm(torch.Tensor(np.zeros((3, 3))))).sqrt() # Minimization problem: loss = ||f(R, TFA)||
            optim.zero_grad() # Reset gradients
            loss.backward()
            optim.step() # Adam step

        R = giraffe_model.R.detach().numpy()
        print(R)
        TFA = torch.abs(giraffe_model.TFA.detach()).numpy()
        return R, TFA

    def get_regulation(self):
        return self._R

    def get_tfa(self):
        return self._TFA

class Model(nn.Module):
    def __init__(
            self,
            tfa_prior,
            motif
    ):
        super().__init__()
        self.TFA = nn.Parameter(torch.Tensor(tfa_prior))
        self.R = nn.Parameter(torch.Tensor(motif))
        self.coefs = nn.Parameter(torch.Tensor(torch.ones(motif.shape[0])))

    def forward(
            self,
            Y,
            PPI,
            C,
            lam,
            adjusting
    ):
        if adjusting is None :
            L1 = torch.norm(Y - torch.matmul(self.R, torch.abs(self.TFA))) ** 2
            L2 = torch.norm(torch.matmul(torch.t(self.R), self.R) - PPI) ** 2
            L3 = torch.norm(torch.matmul(self.R, torch.t(self.R)) - C) ** 2
            L4 = torch.norm(torch.matmul(torch.abs(self.TFA), torch.t(torch.abs(self.TFA))) - PPI) ** 2
            L5 = torch.norm(self.R) ** 2
            weights = self._get_weights(lam, L1, L2, L3)
            return weights[0] * L1 + weights[1] * L2 + weights[2] * L3 + weights[3] * L4 + weights[4] * L5

    def _get_weights(self, lam, L1, L2, L3):
        weights = [1, 1, 1, 1, 1]
        if lam is not None:
            weights = lam
        else:
            sum = L1.item() + L2.item() + L3.item()
            weights = [1 - L1.item() / sum, 1 - L2.item() / sum, 1 - L3.item() / sum, 1, 1]
        return weights
