
#'@param n_qtl The sample size of the reference qtl dataset
#'@param n_gwas The sample size of the GWAS dataset
#'@param snps_per_gene The total number of cis SNPs for each gene
#'@param block_sizes How many genes are in each block
#'@param rho_between The correlation between cis SNPs of different genes
#'@param rho_within The correlation of cis SNPs for the same gene
#'@param response_types A vector containing the number of continuous, binary and count variables (in that order)
#'@param r The parameter controlling the negative binomial distribution for count variables
#'@param k The total number of factors
#'@param maf The minor allele frequency
#'@param sigma2 (what is this again?)
#'@param choice (what is this again?)
#'@param sig_rho (what is this again?)
#'@param binflate Inflation factor for the beta coefficients for the covariates. Used to control SNP heritability on gene expression
#'@param vinfale Inflation factor for outcome variance. Used to control heritability
#'@param pi0 The underlying probability of selecting any given group of genes
#'@param pi1 The underlying probability of selecting any gene in a group
#'@param pi2 The underlying probability of selecting any factor for a gene
#'@param expression_heritability Used to control how much of an effect gene expression has on phenotype heritability
gen_twas_sets <- function(n_qtl, n_gwas, n_genes, snps_per_gene, block_sizes, rho_between, rho_within, response_types = c(1,1,1), r = NULL, k, maf, sigma2 = NULL, choice = "Z", sig_rho = 0.66, binflate = 1, vinflate = 1, pi0, pi1, pi2, expression_heritability = 0.1) {

  response_vec <- rep(c("continuous", "binary", "count"), response_types)
  p <- sum(response_types)
  r2 <- numeric(p)
  if(is.null(r)) {
    r2[which(response_vec == "count")] <- 50
  } else  {
    r2[which(response_vec == "count")] <- r
  }
  r <- r2

  n <- n_qtl + n_gwas

  xp <- sum(block_sizes) #total number of covariates
  LD_corr_matrix <- matrix(rho_between, xp, xp)
  row_start <- 1
  col_start <- 1


  for (size in block_sizes) {
    block_matrix <- matrix(rho_within, nrow = size, ncol = size)

    LD_corr_matrix[row_start:(row_start + size - 1), col_start:(col_start + size - 1)] <- block_matrix

    row_start <- row_start + size
    col_start <- col_start + size
  }

  diag(LD_corr_matrix) <-1

  mu <- numeric(xp)
  Z <- mvrnorm(n = n, mu, Sigma = LD_corr_matrix)

  prob_g0 <- dbinom(0, 2, maf)
  prob_g1 <- dbinom(1, 2, maf)
  prob_g2 <- dbinom(2, 2, maf)

  c1 <- numeric(xp)
  c2 <- numeric(xp)
  for (j in 1:xp) {

    c1[j] <- qnorm(prob_g0, mean = 0, sd = LD_corr_matrix[j,j])

    c2[j] <- qnorm(prob_g2, mean = 0, sd = LD_corr_matrix[j,j], lower.tail = FALSE)
  }


  G <- Z
  intervals <- c(-Inf, c1[1], c2[1], Inf)
  replacement_values <- c(0, 1, 2)
  for (i in seq_along(intervals)) {
    lower_bound <- intervals[i]
    upper_bound <- intervals[i + 1]

    G[(Z > lower_bound) & (Z <= upper_bound)] <- replacement_values[i]
  }

  G_qtl <- G[1:n_qtl,]
  G_gwas <- G[seq(n_qtl+1,n,1),]


  gene_mat <- matrix(nrow = n_qtl, ncol = n_genes)

  snp_sets <- seq(1, n_genes*snps_per_gene, by = snps_per_gene)
  snp_coef_mat <- matrix(0, nrow = snps_per_gene * n_genes, ncol = n_genes)
  counter <- 1
  for(i in snp_sets) {
    snp_coefs <- snp_coef_mat[seq(i,i+(snps_per_gene-1),1),counter] <- rnorm(n = snps_per_gene, sd = 1)
    eta <- as.numeric(G_qtl[,seq(i,i+(snps_per_gene-1),1)] %*% snp_coefs)
    gene_mat[,counter] <- rnorm(n_qtl, mean = eta, sd = sqrt(1-expression_heritability))
    counter <- counter + 1
  }

  group <- rep(1:length(block_sizes), times = block_sizes)

  Sig <- diag(k)
  Sig[Sig == 0] <- sig_rho


  if(choice == "G") {
    alpha <- rbinom(n = length(block_sizes), size = 1, prob = pi0)
    gamma <- rbinom(n = xp, size = 1, prob = pi1)
    omega <- rbinom(n = xp * k, size = 1, prob = pi2)

    amat <- matrix(rep(rep(alpha, times = table(group)), k), nrow = xp, ncol = k, byrow = FALSE)
    gmat <- matrix(rep(gamma, k), nrow = xp, k, byrow = FALSE)
    omat <- matrix(omega, ncol = k, byrow = FALSE)
    b <- MASS::mvrnorm(n = xp, mu = rep(0, k), Sigma = binflate * vinflate * Sig)
    beta_true <- amat * gmat * omat * b
  } else {
    gamma <- rbinom(n = n_genes, size = 1, prob = pi1)
    omega <- rbinom(n = n_genes * k, size = 1, prob = pi2)

    gmat <- matrix(rep(gamma, k), nrow = n_genes, k, byrow = FALSE)
    omat <- matrix(omega, ncol = k, byrow = FALSE)
    b <- MASS::mvrnorm(n = n_genes, mu = rep(0, k), Sigma = binflate * vinflate * Sig)
    beta_true <- gmat * omat * b
  }

  E_gwas <- G_gwas %*% snp_coef_mat

  A <- as.numeric(rowSums(abs(beta_true)) != 0)
  if(choice == "G") {
    XB <- G_gwas %*% beta_true
  } else {
    XB <- E_gwas %*% beta_true
  }


  Y_latent <- XB + MASS::mvrnorm(n = n_gwas, mu = rep(0,k), Sigma = vinflate * Sig)

  Lambda.true<- matrix(0,p,k)

  if(k == p) {
    for(h in 1:k){
      Lambda.true[sample(1:p,p),h] <- rnorm(k,0,1)
    }
  } else {
    for(h in 1:k){
      Lambda.true[sample(1:p,p)[1:(2*k-(h-1))],h] <- rnorm(2*k-(h-1),0,1)
    }
  }



  if(is.null(sigma2)) {
    Ucov <- diag(1/rgamma(p,shape=1,scale=0.25))
  } else {
    Ucov <- sigma2 * diag(p)
  }
  U <- MASS::mvrnorm(n = n_gwas, rep(0, p), Ucov)

  theta <- Y_latent %*% t(Lambda.true) + U


  Y_observed <- matrix(0, n_gwas, p)
  sigmoid <- function(x) {
    1/(1 + exp(-x))
  }
  for(q in 1:length(response_vec)) {
    if(response_vec[q] == "continuous") {
      Y_observed[,q] <- matrix(Y_latent %*% Lambda.true[q,], n_gwas, 1) + U[,q] + rnorm(n_gwas, 0, 1)
    } else if(response_vec[q] == "binary") {
      prob <- sigmoid(matrix(Y_latent %*% Lambda.true[q, ], n_gwas, 1) + U[, q])
      Y_observed[,q] <- rbinom(n_gwas, 1, prob)
    } else {
      prob <- sigmoid(matrix(Y_latent %*% Lambda.true[q, ], n_gwas, 1) + U[, q])
      prob[prob > 0.9999] <- 0.9999
      Y_observed[,q] <- stats::rnbinom(n_gwas, size = r[q], prob = 1-prob)
    }
  }


  return(list(G_qtl = G_qtl, G_gwas = G_gwas, E_qtl = gene_mat, SNP_coef = snp_coef_mat, Yl = Y_latent, Yo = Y_observed, Lambda = Lambda.true, Sigma = vinflate * Sig,
              group = group, A = A, beta = beta_true, r = r, responses = response_vec, U = U, Ucov = diag(Ucov), LD = LD_corr_matrix))

}

