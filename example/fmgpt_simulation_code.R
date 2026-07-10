#This code generates an example data set following the same process for Scenario 1
#when there are 50 continuous outcome variables and only 1 underlying factor
library(BoomSpikeSlab)
library(pROC)
source("R/gen_gene_by_snp.R")

#generate 100 genes, 10 cis-SNPs per gene, 10 genes per region

dat <- gen_twas_sets(n_qtl = 250, n_gwas = 5000, n_genes = 100, snps_per_gene = 10, block_sizes = rep(10, 100), rho_between = 0.33, rho_within = 0.66, response_types = c(50, 0, 0), r = 0.1, k = 1, maf = 0.24, choice = "genes", sig_rho = 0.66, pi0 = 1, pi1 = 0.1, pi2 = 1, expression_heritability = 0.1,sigma2 = 1, binflate = 0.05)
while(sum(dat$A) == 0) {
  dat <- gen_twas_sets(n_qtl = 250, n_gwas = 5000, n_genes = 100, snps_per_gene = 10, block_sizes = rep(10, 100), rho_between = 0.33, rho_within = 0.66, response_types = c(50, 0, 0), r = 0.1, k = 1, maf = 0.24, choice = "genes", sig_rho = 0.66, pi0 = 1, pi1 = 0.1, pi2 = 1, expression_heritability = 0.1,sigma2 = 1, binflate = 0.05)
}

est_coef_mat <- matrix(nrow = nrow(dat$SNP_coef), ncol = ncol(dat$SNP_coef))

for(j in 1:ncol(dat$SNP_coef)) {
  tmp_out <- lm.spike(scale(dat$E_qtl[,j], scale = FALSE) ~ scale(dat$G_qtl, scale = FALSE)+0, niter = 3000)
  est_coef_mat[,j] <- apply(tmp_out$beta[2000:3000,], 2, median)
}

E_gwas <- dat$G_gwas %*% est_coef_mat

#now split the data into the 10 distinct regions

split_idx <- rep(c(1:10), each = 10)

fit_list <- vector("list", 10)

for(j in 1:10) {
  fit_list[[j]] <- FMGPT(X = scale(E_gwas[,which(split_idx == j)], scale = FALSE), Y = scale(dat$Yo, scale = FALSE), A = 0*dat$A[which(split_idx == j)], group = seq(1, 10), response_types = dat$responses, epsilon = 0.4, B = 3000, burnin = 1000, thin = 1, k_fix = NULL)
}

#get all pips

pips <- unlist(lapply(fit_list, FUN = function(x) twas_pip(x, max(x$kstar))))

#calculate, for example, the AUC
auc(roc(dat$A, pips))

                    
#calculate power, fdr
if(all(pips == 0)) {
    fdr <- 0
    power <- 0
    true_detect[i] <- 0
    false_detect[i] <- 0
  } else {
    cutoff <- bfdr(pips, 0.1)
    signal <- pips >= cutoff
    A <- dat$A
    false_discovery <- sum(A == 0 & signal == 1)
    discoveries <- sum(signal)
    if(discoveries == 0) {
      fdr <- 0
    } else {
      fdr <- false_discovery / discoveries
    }
    
    false_negative <- sum(A == 1 & signal == 0)
    negatives <- sum(1 - signal)
    power <- sum(signal*A) / sum(A)
}
