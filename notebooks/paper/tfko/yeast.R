library(netZooR)
library(ramify)
library(ROCR)
library(lars)
library(devtools)
install_github("jpvert/tigress")
library(tigress)

setwd("/Users/soel/Desktop/giraffe/notebooks/paper/tfko/")

BITFAM <- function(data, tf, chipseq_weight, interseted_TF = NA, scATAC_obj = NA, ncores = 1, iter = 8000, tol_rel_obj=0.005){
  
  TF_used <-tf
  gene_list <- rownames(data)
  
  if(is.na(scATAC_obj)){
  }else{
    for(i in TF_used){
      gene_list[[which(TF_used == i)]] <- gene_list[[which(TF_used == i)]][gene_list[[which(TF_used == i)]] %in% BITFAM_scATAC(scATAC_obj)]
    }
  }
  
  X <- t(as.matrix(data))
  
  
  Mask_matrix <- chipseq_weight
  X <- t(as.matrix(data))
  N <- dim(X)[1]
  D <- dim(X)[2]
  K <- length(TF_used)
  data_to_model <- list(N = N, D = D, K = K, X = X, Mask = Mask_matrix)
  
  
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = 1)
  
  set.seed(100)
  pca_beta_piror <- "
data {
int<lower=0> N; // Number of samples
int<lower=0> D; // The original dimension
int<lower=0> K; // The latent dimension
matrix[N, D] X; // The data matrix
matrix[D, K] Mask; // The binary mask of prior knowledge indicate the target of TFs
}
parameters {
matrix<lower=0, upper=1>[N, K] Z; // The latent matrix
matrix[D, K] W; // The weight matrix
real<lower=0> tau; // Noise term
vector<lower=0>[K] alpha; // ARD prior
}
transformed parameters{
matrix<lower=0>[D, K] t_alpha;
real<lower=0> t_tau;
for(wmd in 1:D){
for(wmk in 1:K){
t_alpha[wmd, wmk] = Mask[wmd, wmk] == 1 ? inv(sqrt(alpha[wmk])) : 0.01;
}
}
t_tau = inv(sqrt(tau));
}
model {
tau ~ gamma(1,1);
to_vector(Z) ~ beta(0.5, 0.5);
alpha ~ gamma(1e-3,1e-3);
for(d in 1:D){
for(k in 1:K){
W[d,k] ~ normal(0, t_alpha[d, k]);
}
}
to_vector(X) ~ normal(to_vector(Z*W'), t_tau);
} "

m_beta_prior <- stan_model(model_code = pca_beta_piror)
suppressWarnings(fit.vb <- vb(m_beta_prior, data = data_to_model, algorithm = "meanfield",
                              iter = iter, output_samples = 100, tol_rel_obj = tol_rel_obj))

result_matrix <- apply(extract(fit.vb,"W")[[1]], c(2,3), mean)
rownames(result_matrix) <- rownames(data)
colnames(result_matrix) <- TF_used
return(result_matrix)
}

E <- read.csv("E.csv")
rownames(E) <- E[,1]
E <- E[,-1]
full_E <- read.csv("full_E.csv")
rownames(full_E) <- full_E[,1]
full_E <- full_E[,-1]
M <- read.csv("M.csv")
rownames(M) <- M[,1]
M <- M[,-1]
D <- read.csv("D.csv")
rownames(D) <- D[,1]
D <- D[,-1]
gt <- D
expression <- read.csv("../yeast/expression/TFKOexpressionData.csv")
rownames(expression) <- expression[,1]
expression <- expression[,-1]

C <- cor(t(expression))
C <- C[rownames(expression) %in% rownames(M),colnames(M)]
write.csv(C, "R_corr.csv")

tiger <- netZooR::tiger(E, abs(t(M)), TFexpressed = FALSE, signed = FALSE, seed=46)
write.csv(tiger$W, "R_tiger.csv")

R <- BITFAM(E, tf=colnames(M), M)
write.csv(R, "R_bitfam.csv")


R_tigress <- tigress(
  t(expression), 
  tf=colnames(M), 
  targetlist = rownames(expression),
  allsteps = FALSE
)
R_tigress <- t(R_tigress)[rownames(E),]
hit <- sum(as.numeric(unlist(flatten(D))) * flatten(R_tigress) < 0)
miss <- sum(as.numeric(unlist(flatten(D))) * flatten(R_tigress) > 0)
print(hit / (hit + miss))
write.csv(R_tigress, "R_tigress.csv")