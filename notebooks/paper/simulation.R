library(netZooR)
library(ramify)
library(ROCR)
library(lars)
library(devtools)
install_github("jpvert/tigress")
library(tigress)

setwd("/Users/soel/Desktop/giraffe/notebooks/paper")

BITFAM <- function(data, tf, interseted_TF = NA, scATAC_obj = NA, ncores = 1, iter = 8000, tol_rel_obj=0.005){
  
  TF_used <-tf
  gene_list <- rownames(expr)
  
  if(is.na(scATAC_obj)){
  }else{
    for(i in TF_used){
      gene_list[[which(TF_used == i)]] <- gene_list[[which(TF_used == i)]][gene_list[[which(TF_used == i)]] %in% BITFAM_scATAC(scATAC_obj)]
    }
  }
  
  X <- t(as.matrix(data))
  chipseq_weight <- matrix(1, nrow = length(colnames(X)), ncol = length(TF_used))
  for(i in 1:length(TF_used)){
    chipseq_weight[, i] <- ifelse(colnames(X) %in% gene_list[[i]], 1, 0)
  }
  
  
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

shuffling = c("0.01", "0.1", "0.2", "0.3", "0.4", "0.5")
B <- 50
res <- c()
for(i in seq_len(length(shuffling))){
  S <- c()
  for(j in seq_len(B)){
    prior <- read.csv(paste0("simulation_data/prior_sim_", (j - 1), "_shuffle_", shuffling[i], ".csv"))
    prior <- prior[c(seq(1, dim(prior)[1], 1)), c(seq(2, dim(prior)[2], 1))]
    expr <- read.csv(paste0("simulation_data/gene_expression_sim_", (j - 1), "_shuffle_", shuffling[i], ".csv"))
    expr <- expr[c(seq(101, dim(expr)[1], 1)),c(seq(2, dim(expr)[2], 1))]
    R <- read.csv(paste0("simulation_data/R_sim_", (j - 1), "_shuffle_", shuffling[i], ".csv"))
    R <- R[c(seq(1, dim(R)[1], 1)), c(seq(2, dim(R)[2], 1))]
    rownames(prior) <- rownames(expr)
    
    R_tiger <- netZooR::tiger(expr, t(prior), TFexpressed = FALSE, signed = FALSE)$W
    R[R != 0] <- 1
    methodPred  <- prediction(flatten(as.matrix(R_tiger)), flatten(as.matrix(R)))
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    S <- c(S, auc.methodPred)
    print(auc.methodPred)
  }
  res <- c(res, S)
}

start <- 1
while(TRUE){
  curr <- res[seq(start, start + 49)]
  print(mean(curr))
  print(sd(curr))
  start <- start + 50
  if(start > length(res)){
    break
  }
}

res <- c()
for(j in seq_len(B)){
    expr <- read.csv(paste0("simulation_data/gene_expression_sim_", (j - 1), "_shuffle_", shuffling[i], ".csv"))
    expr <- expr[c(seq(101, dim(expr)[1], 1)),c(seq(2, dim(expr)[2], 1))]
    R <- read.csv(paste0("simulation_data/R_sim_", (j - 1), "_shuffle_", shuffling[i], ".csv"))
    R <- R[,c(seq(2, dim(R)[2], 1))]
    
    R_bitfam <- BITFAM(expr, tf=rownames(expr)[c(seq(1, 100, 1))])
    
    R[R != 0] <- 1
    methodPred  <- prediction(flatten(as.matrix(R_bitfam)), flatten(as.matrix(R)))
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    res <- c(res, auc.methodPred)
}

res <- c()
for(j in seq_len(B)){
  print(j)
  expr <- read.csv(paste0("simulation_data/gene_expression_sim_", (j - 1), "_shuffle_", shuffling[1], ".csv"))
  expr <- expr[c(seq(1, dim(expr)[1], 1)),c(seq(2, dim(expr)[2], 1))]
  R <- read.csv(paste0("simulation_data/R_sim_", j - 1, "_shuffle_", shuffling[1], ".csv"))
  R <- R[,c(seq(2, dim(R)[2], 1))]
  
  R_tigress <- tigress(
    t(expr) + matrix(rnorm(600 * 50, sd = 0.001), nrow = 50), 
    tf=rownames(expr)[c(seq(1, 100, 1))], 
    targetlist = rownames(expr)[c(seq(101, 600, 1))],
    allsteps = FALSE
  )
  R_tigress <- t(R_tigress)
  
  R[R != 0] <- 1
  methodPred  <- prediction(flatten(as.matrix(R_tigress)), flatten(as.matrix(R)))
  auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
  res <- c(res, auc.methodPred)
}
