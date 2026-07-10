
update_Omega5 <- function(eta, Lambda, U, n, Q, response_types, Y, r) {
  omat <- matrix(NA, n, Q)
  theta <- eta %*% Lambda + U
  for(q in 1:Q) {
    if(response_types[q,drop = FALSE] == "continuous") {
      omat[,q] <- 1
    } else if(response_types[q,drop = FALSE] == "binary") {
      omat[,q] <- pg::rpg_hybrid(rep(1,n), theta[,q,drop = FALSE])
    } else {
      omat[,q] <- pg::rpg_hybrid(Y[,q,drop = FALSE]+r[q,drop = FALSE], theta[,q,drop = FALSE])
    }
  }
  return(omat)
}

update_Z <- function(Omega, response_types, Y, r, n, Q) {
  zmat <- matrix(NA, n, Q)
  for(q in 1:Q) {
    if(response_types[q,drop = FALSE] == "continuous") {
      zmat[,q] <- Y[,q,drop = FALSE]
    } else if(response_types[q,drop = FALSE] == "binary") {
      zmat[,q] <- (Y[,q,drop = FALSE] - 1/2) / Omega[,q,drop = FALSE]
    } else {
      zmat[,q] <- (Y[,q,drop = FALSE] - (Y[,q,drop = FALSE] + r[q,drop = FALSE])/2) / Omega[,q,drop = FALSE]
    }
  }
  return(zmat)
}


update_a1 <- function(a1_cur,delta1){
  a1_new <- rtruncnorm(1,a=1,b=Inf,mean=a1_cur,sd=1)
  calc1 <- (sum(dgamma(x = delta1, shape = a1_new, rate = 1, log = TRUE))+dgamma(x = a1_new, shape = 2, rate = 1, log = TRUE)) - (sum(dgamma(x = delta1, shape = a1_cur, rate = 1, log = TRUE)) + dgamma(x = a1_cur, shape = 2, rate = 1, log = TRUE))
  calc2 <- log(dtruncnorm(x = a1_cur, a = 1, b = Inf, mean = a1_new, sd = 1)) - log(dtruncnorm(a1_new, a = 1, b = Inf, mean = a1_cur, sd = 1))
  acceptance_prob <- min(1, exp(calc1 + calc2))
  if(is.na(acceptance_prob)) {
    return(a1_cur)
  } else {
    sample(c(a1_new, a1_cur), 1, prob = c(acceptance_prob, 1 - acceptance_prob))
  }
}


update_a2 <- function(a2_cur,deltah){
  a2_new <- rtruncnorm(1,a=1,b=Inf,mean=a2_cur,sd=1)
  calc1 <- (sum(dgamma(x = deltah, shape = a2_new, rate = 1, log = TRUE))+dgamma(x = a2_new, shape = 2, rate = 1, log = TRUE)) - (sum(dgamma(x = deltah, shape = a2_cur, rate = 1, log = TRUE)) + dgamma(x = a2_cur, shape = 2, rate = 1, log = TRUE))
  calc2 <- log(dtruncnorm(x = a2_cur, a = 1, b = Inf, mean = a2_new, sd = 1)) - log(dtruncnorm(a2_new, a = 1, b = Inf, mean = a2_cur, sd = 1))
  acceptance_prob <- min(1, exp(calc1 + calc2))
  if(is.na(acceptance_prob)) {
    return(a2_cur)
  } else {
    sample(c(a2_new, a2_cur), 1, prob = c(acceptance_prob, 1 - acceptance_prob))
  }
}

update_phi <- function(v, tauh, Lambda, p, Q) {
  pmat <- matrix(NA, p, Q)
  for(k in 1:p) {
    for(q in 1:Q) {
      pmat[k,q] <- rgamma(n = 1, shape = (v+1)/2, rate = (v+tauh[k]*Lambda[k,q,drop = FALSE]^2)/2)
    }
  }
  return(pmat)
}

update_phijh <- function(v,tauh,lambdajh){
  return(rgamma(n=1,
                shape=(v+1)/2,
                rate=(v+tauh*lambdajh^2)/2))
}

update_delta1 <- function(a1,Q,p,tau1l,phi,Lambda){
  return(rgamma(n=1,
                shape=a1+Q*p/2,
                rate=1+0.5*sum(tau1l*rowSums(phi*Lambda^2))))
}


update_deltah <- function(a2, Q, p, h, phi, Lambda, delta) {
  tauvec <- numeric(length(h:p))
  counter <- 1
  for(l in h:p) {
    tauvec[counter] <- taulh(delta, l, h)
    counter <- counter + 1
  }
  return(rgamma(n = 1,
                shape = a2 + Q/2*(p-h+1), rate = 1 + 1/2 * sum(tauvec * rowSums(phi[h:p,,drop = FALSE] * Lambda[h:p,,drop = FALSE]^2))))
}

update_lambda <- function(Omega, eta, Z, U, phi, delta1, deltah, p, Q) {
  lmat <- matrix(NA, p, Q)
  for(q in 1:Q) {
    if(p == 1) {
      invDq <- (diag(phi[1:p,q,drop = FALSE]*cumprod(c(delta1))))
    } else {
      invDq <- (diag(c(phi[1:p,q,drop = FALSE]*cumprod(c(delta1,deltah[1:(p-1),drop = FALSE])))))
    }
    V <- solve(invDq + eigenMapMatMult(t(eta), (eta * Omega[,q]))) #p x p
    M <- eigenMapMatMult(V, eigenMapMatMult(t(eta * Omega[,q]), (Z[,q,drop = FALSE] - U[,q,drop = FALSE])))
    lmat[,q] <- mvnfast::rmvn(n = 1, mu = M, sigma = V)
  }
  return(lmat)
}

update_sigma2 <- function(a_sigma, b_sigma, n, U, Q) {
  s2mat <- numeric(Q)
  for(q in 1:Q) {
    s2mat[q] <- LaplacesDemon::rinvgamma(n = 1, a_sigma + n / 2, b_sigma + sum(U[,q,drop = FALSE]^2)/2)
  }
  return(s2mat)
}

taulh <- function(delta, l, h) {
  tau <- prod(delta[1:l])
  if(h <= l) {
    tau <- tau / delta[h]
  }
  return(tau)
}


b_update <- function(X, Y, Z, beta, bsigma, Sigma, XtX, XtY) {
  tmp_b <- matrix(nrow = ncol(X), ncol = ncol(Y))
  for(colx in 1:ncol(X)) {
    if(all(Z[colx,,drop = FALSE] == 0)) {
      m_gj <- rep(0, ncol(Y))
      sig_gj <- (bsigma * Sigma)
    } else {
      m_gj <- as.numeric(solve(XtX[colx,colx,drop = FALSE] + 1/bsigma)) * (as.numeric(X[,colx,drop = FALSE]) %*% (Y - eigenMapMatMult(X[,-colx,drop = FALSE], beta[-colx,,drop = FALSE])))
      m_gj <- Z[colx,,drop = FALSE] * m_gj
      ZZt <- Z[colx,] %*% t(Z[colx,])
      sig_gj <- (ZZt) * (as.numeric(solve(XtX[colx,colx] + (1/bsigma)))) * Sigma + (1 - ZZt) * (bsigma * Sigma)
      sig_gj[which(Z[colx,,drop = FALSE]==1),which(Z[colx,,drop = FALSE]!=1)] <- sqrt((1/bsigma)  * as.numeric(solve(XtX[colx,colx,drop = FALSE] + 1/bsigma))) * (bsigma * Sigma)[which(Z[colx,,drop = FALSE]==1),which(Z[colx,,drop = FALSE]!=1),drop = FALSE]
      sig_gj[which(Z[colx,,drop = FALSE]!=1),which(Z[colx,,drop = FALSE]==1)] <- sqrt((1/bsigma)  * as.numeric(solve(XtX[colx,colx,drop = FALSE] + 1/bsigma))) * (bsigma * Sigma)[which(Z[colx,,drop = FALSE]!=1),which(Z[colx,,drop = FALSE]==1),drop = FALSE]
    }
    tmp_b[colx,] <- mvnfast::rmvn(n = 1, mu = m_gj, sigma = sig_gj)
  }
  return(tmp_b)
}

bsigma_update <- function(X, Y, b, beta, Sigma, Z, t1, t2) {
  return(LaplacesDemon::rinvgamma(n = 1, shape = t1 + 0.5 * (sum(rowSums(Z))), scale = t2 + 0.5 * sum(eigenMapMatMult(beta, solve(Sigma)) * beta)))
}

gpi_update <- function(alpha) {
  return(rbeta(n = 1, shape1 = 1 + sum(alpha), shape2 = 1 + sum(1-alpha)))
}

chi_update <- function(gamma, AMAT, a) {
  mu <- c(AMAT %*% a)
  potential_chis <- rnorm(length(gamma), mu, sd = 1)
  return(ifelse(gamma == 0, -abs(potential_chis), abs(potential_chis)))
}

a_update <- function(AMAT, chi) {
  return(c(MASS::mvrnorm(n = 1, mu = solve(eigenMapMatMult(t(AMAT), AMAT) + diag(2) * (1/10)) %*% (eigenMapMatMult(t(AMAT), chi) + (diag(2) * (1/10)) %*% c(0,0)), Sigma = solve(eigenMapMatMult(t(AMAT), AMAT) + diag(2) * (1/10)))))
}

spi_update <- function(AMAT, a) {
  return(pnorm(0, mean = eigenMapMatMult(AMAT, a), sd = 1, lower.tail = FALSE))
}

tpi_update <- function(omega) {
  return(rbeta(n = 1, shape1 = 1 + sum(omega), shape2 = 1 + sum(1-omega)))
}


alpha_update <- function(X, Y, b, group, pi0, alpha, gamma, omega, Sigma) {
  gunique <- unique(group)
  gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
  omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
  SSigma <- solve(Sigma)
  for(g in 1:length(gunique)) {
    if(g %in% which(table(group)==1)) {
      alpha[g] <- 1
    } else {
      amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
      gg <- which(group == g)
      amatg <- 0 * amat
      amatg[gg,] <- 1
      amatng <- amat
      amatng[gg,] <- 0
      Zg <- amatg * gmat * omat
      Zng <- amatng * gmat * omat
      Yng <- Y - eigenMapMatMult(X, (Zng * b))
      Bg <- Zg * b
      XB <- eigenMapMatMult(X, Bg)
      calc1 <- (-0.5) * sum(eigenMapMatMult(Yng - XB, SSigma) * (Yng - XB))
      calc2 <- (-0.5) * sum(eigenMapMatMult(Yng, SSigma) * (Yng))
      prob <- pi0 / (pi0 + (1-pi0) * exp(calc2 - calc1))
      alpha[g] <- rbinom(1, 1, prob)
    }
  }
  return(alpha)
}

gamma_update <- function(X, Y, b, group, pi1, alpha, gamma, omega, Sigma) {
  amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
  omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
  SSigma <- solve(Sigma)
  for(xcol in 1:ncol(X)) {
    gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
    gmats <- gmat * 0
    gmats[xcol,] <- 1
    gmatns <- gmat
    gmatns[xcol,] <- 0
    Zs <- amat * gmats * omat
    Zns <- amat * gmatns * omat
    Yns <- Y - eigenMapMatMult(X, (Zns * b))
    Bs <- Zs * b
    XB <- eigenMapMatMult(X, Bs)
    calc1 <- (-0.5) * sum(eigenMapMatMult(Yns - XB, SSigma) * (Yns - XB))
    calc2 <- (-0.5) * sum(eigenMapMatMult(Yns, SSigma) * (Yns))
    prob <- pi1[xcol,drop = FALSE] / (pi1[xcol,drop = FALSE] + (1-pi1[xcol,drop = FALSE]) * exp(calc2 - calc1))
    gamma[xcol] <- rbinom(1, 1, prob)
  }
  return(gamma)
}


omega_update <- function(X, Y, b, group, pi2, alpha, gamma, omega, Sigma) {
  amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
  gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
  SSigma <- solve(Sigma)
  tracker <- 1
  for(ycol in 1:ncol(Y)) {
    for(xcol in 1:ncol(X)) {
      omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
      omatt <-  0 * omat
      omatt[xcol,ycol] <- 1
      omatnt <- omat
      omatnt[xcol,ycol] <- 0
      Zt <- amat * gmat * omatt
      Znt <- amat * gmat * omatnt
      Ynt <- Y - eigenMapMatMult(X, (Znt * b))
      Bt <- Zt * b
      XB <- eigenMapMatMult(X, Bt)
      calc1 <- (-0.5) * sum(eigenMapMatMult(Ynt - XB, SSigma) * (Ynt - XB))
      calc2 <- (-0.5) * sum(eigenMapMatMult(Ynt, SSigma) * Ynt)
      prob <- pi2 / (pi2 + (1-pi2) * exp(calc2 - calc1))
      omega[tracker] <- rbinom(1,1, prob)
      tracker <- tracker + 1
    }
  }
  return(omega)
}

Z_update <- function(X, Y, group, alpha, gamma, omega) {
  amat <- matrix(rep(rep(alpha, times = table(group)), ncol(Y)), nrow = ncol(X), ncol = ncol(Y), byrow = FALSE)
  gmat <- matrix(rep(gamma, ncol(Y)), nrow = ncol(X), ncol(Y), byrow = FALSE)
  omat <- matrix(omega, ncol = ncol(Y), byrow = FALSE)
  return(amat * gmat * omat)
}

beta_update <- function(b, Z) {
  return(Z * b)
}


Sigma_update <- function(X, Y, b, beta, Sigma, Z, bsigma) {
  idx <- which(rowSums(Z)!=0)
  return(LaplacesDemon::rinvwishart(nu = length(idx) + ncol(Y) + nrow(Y), S = diag(ncol(Y)) + eigenMapMatMult(t(Y - eigenMapMatMult(X, beta)), (Y - eigenMapMatMult(X, beta))) + (1/bsigma) * eigenMapMatMult(t(beta), beta)))
}

update_t1_bak <- function(s2, tcur, t2) {
  tstar <- rtruncnorm(n = 1, a = 0, b = Inf, mean = tcur, sd = 1)
  num <- sum(LaplacesDemon::dinvgamma(x=s2,shape=tstar,scale=t2, log = TRUE))+LaplacesDemon::dinvgamma(x=tstar,shape=2,scale=1, log = TRUE)
  denom <- sum(LaplacesDemon::dinvgamma(x=s2,shape=tcur,scale=t2, log = TRUE))+LaplacesDemon::dinvgamma(x=tcur,shape=2,scale=1, log = TRUE)
  r <- min(1,exp(num-denom))
  if(r>runif(1)) {
    tnew <- tstar
  } else {
    tnew <- tcur
  }
  return(tnew)
}

update_t1 <- function(s2, tcur, t2){
  tstar <- rtruncnorm(n = 1, a = 0, b = Inf, mean = tcur, sd = 1)
  calc1 <- (sum(LaplacesDemon::dinvgamma(x = s2, shape = tstar, scale = t2, log = TRUE)) + LaplacesDemon::dinvgamma(x = tstar, shape = 2, scale = 1, log = TRUE)) - (sum(LaplacesDemon::dinvgamma(x = s2, shape = tcur, scale = t2, log = TRUE)) + LaplacesDemon::dinvgamma(x = tcur, shape = 2, scale = 1, log = TRUE))
  calc2 <- log(dtruncnorm(x = tcur, a = 1, b = Inf, mean = tstar, sd = 1)) - log(dtruncnorm(tstar, a = 1, b = Inf, mean = tcur, sd = 1))
  acceptance_prob <- min(1, exp(calc1 + calc2))
  if(is.na(acceptance_prob)) {
    return(tcur)
  } else {
    sample(c(tstar, tcur), 1, prob = c(acceptance_prob, 1 - acceptance_prob))
  }
}

update_t2_bak <- function(s2, tcur, t1) {
  tstar <- rtruncnorm(n = 1, a = 0, b = Inf, mean = tcur, sd = 1)
  num <- sum(LaplacesDemon::dinvgamma(x=s2,shape=t1,scale=tstar, log = TRUE))+LaplacesDemon::dinvgamma(x=tstar,shape=2,scale=1, log = TRUE)
  denom <- sum(LaplacesDemon::dinvgamma(x=s2,shape=t1,scale=tcur, log = TRUE))+LaplacesDemon::dinvgamma(x=tcur,shape=2,scale=1, log = TRUE)
  r <- min(1,exp(num-denom))
  if(r>runif(1)) {
    tnew <- tstar
  } else {
    tnew <- tcur
  }
  return(tnew)
}

update_t2 <- function(s2, tcur, t1){
  tstar <- rtruncnorm(n = 1, a = 0, b = Inf, mean = tcur, sd = 1)
  calc1 <- (sum(LaplacesDemon::dinvgamma(x = s2, shape = t1, scale = tstar, log = TRUE)) + LaplacesDemon::dinvgamma(x = tstar, shape = 2, scale = 1, log = TRUE)) - (sum(LaplacesDemon::dinvgamma(x = s2, shape = t1, scale = tcur, log = TRUE)) + LaplacesDemon::dinvgamma(x = tcur, shape = 2, scale = 1, log = TRUE))
  calc2 <- log(dtruncnorm(x = tcur, a = 1, b = Inf, mean = tstar, sd = 1)) - log(dtruncnorm(tstar, a = 1, b = Inf, mean = tcur, sd = 1))
  acceptance_prob <- min(1, exp(calc1 + calc2))
  if(is.na(acceptance_prob)) {
    return(tcur)
  } else {
    sample(c(tstar, tcur), 1, prob = c(acceptance_prob, 1 - acceptance_prob))
  }
}

#'
#' FM-GPT Fine Mapping Function
#'
#'@param X Covariate matrix X
#'@param Y Output matrix Y, possible containing mixed data types
#'@param A Indicator for eQTL status of a SNP
#'@param group Grouping indicator
#'@param respones_types A vector indicating the variable types in Y
#'@param epsilon A numeric giving the cutoff for loadings before a factor is dropped
#'@param B The total number of iterations
#'@param burnin The total number of burnin iterations
#'@param thin The thinning parameter
#'@param k_fix Fixed number of factors to use
#'@param s2.fit Logical indicating whether s^2 should be estimated
#'@param Sigma Logical indicating whether Sigma should be estimated
#'@param scale.eta A logical indicating whether the factor scores should be mean centered
#'
#'@return A list containing the kept iterations for the loading matrix, number of factors,
#'beta, a1 and a2, sigma2, sigma_beta, count data parameter r, adjustment effect for eQTL SNPs
#'@export
FMGPT <- function(X, Y, A, group, response_types, epsilon, B, burnin, thin = 1, k_fix = NULL, s2.fit = TRUE, Sigma.fit = TRUE, scale.eta = FALSE) {
  X <- scale(X, scale = TRUE)
  Q <- length(response_types)
  n <- nrow(Y)
  xp <- ncol(X)
  AMAT <- cbind(1, A)
  ngroup <- length(unique(group))
  post_inds <- seq(1+burnin,B,by=thin)
  store_length <- length(post_inds)

  v <- 3
  a_sigma <- 1
  b_sigma <- 0.3
  if(is.null(k_fix)) {
    kstore <- kstar <- Q
  } else {
    kstore <- kstar <- k_fix
  }
  r.store <- matrix(NA, store_length, Q)
  sigma2.store <- matrix(NA,store_length,Q)
  a1.store <- a2.store <- rep(NA,store_length)
  Lambda.store <- array(NA,dim=c(store_length,kstore,Q))
  bsigma.store <- rep(NA,store_length)
  pi0.store <- rep(NA,store_length)
  pi1.store <- matrix(NA,store_length,xp)
  pi2.store <- rep(NA,store_length)
  Z2.store <- array(NA,dim=c(store_length,xp,kstore))
  beta.store <- array(NA, dim = c(store_length, xp, kstore))
  a.store <- matrix(NA,store_length,2)
  t1.store <- t2.store <- rep(NA, store_length)
  mt.store <- rep(NA, store_length)
  kstar.store <- rep(NA, store_length)
  ut.store <- rep(NA, store_length)
  pt.store <- rep(NA, store_length)
  mt.store <- rep(NA, store_length)

  Omega.init <- matrix(0, n, Q)
  U.init <- matrix(0, n, Q)
  Z.init <- matrix(0, n, Q)
  r.init <- rep(10,Q)
  sigma2.init <- rep(1,Q)
  eta.init <- MASS::mvrnorm(n = n, mu = rep(0, kstore), Sigma = diag(kstore))
  a1.init <- rgamma(1,shape=2,scale=1)
  a2.init <- 3
  delta1.init <- rgamma(1,a1.init,1)
  deltah.init <- rgamma(kstore-1,a2.init,1)
  phi.init <- matrix(rgamma(Q*kstore,shape=v/2,rate=v/2),kstore,Q)
  Lambda.init <- matrix(1,kstore,Q)
  Sigma.init <- LaplacesDemon::rinvwishart(nu = kstore, S = diag(kstore))
  bsigma.init <- runif(1,0,1)
  effb.init <- MASS::mvrnorm(n = xp, mu = rep(0,kstore), Sigma = Sigma.init)
  pi0.init <- runif(1,0,1)
  pi1.init <- runif(xp,0,1)
  pi2.init <- runif(1,0,1)
  alpha.init <- rbinom(ngroup,1,pi0.init)
  gamma.init <- rbinom(xp,1,pi1.init)
  omega.init <- rbinom(xp*kstore,1,pi2.init)
  Z2.init <- matrix(rep(rep(alpha.init, times = table(group)), kstore), nrow = xp, ncol = kstore, byrow = FALSE) *
    matrix(rep(gamma.init, kstore), nrow = xp, ncol = kstore, byrow = FALSE) *
    matrix(omega.init, ncol = kstore, byrow = FALSE)
  beta.init <- effb.init * Z2.init
  chi.init <- rnorm(xp, 0, 1); chi.init <- ifelse(A == 1, abs(chi.init), -abs(chi.init))
  a.init <- rnorm(2,0,1)
  t1.init <- rgamma(n = 1, 2, 1)
  t2.init <- rgamma(n = 1, 2, 1)

  Omega <- Omega.init
  U <- U.init
  Z <- Z.init
  r <- r.init
  sigma2 <- sigma2.init
  eta <- eta.init
  a1 <- a1.init
  a2 <- a2.init
  delta1 <- delta1.init
  deltah <- deltah.init
  phi <- phi.init
  Lambda <- Lambda.init
  Sigma <- Sigma.init
  bsigma <- bsigma.init
  effb <- effb.init
  pi0 <- pi0.init
  pi1 <- pi1.init
  pi2 <- pi2.init
  alpha <- alpha.init
  gamma <- gamma.init
  omega <- omega.init
  Z2 <- Z2.init
  beta <- beta.init
  chi <- chi.init
  a <- a.init
  t1 <- t1.init
  t2 <- t2.init
  XtX <- eigenMapMatMult(t(X), X)

  if(is.null(k_fix)) {
    kstore <- kstar <- floor(log(Q)*3)
  }

  counter <- 1
  for(b in 2:B) {
    print(b)
    start <- proc.time()

    U <- update_U_cpp(Z = Z, Omega = Omega, Sigma = diag(sigma2), eta = eta[,1:kstar,drop = FALSE], Lambda = Lambda[1:kstar,,drop = FALSE], n = n, Q = Q)
    Omega <- update_Omega5(eta = eta[,1:kstar,drop = FALSE], Lambda = Lambda[1:kstar,,drop = FALSE], U = U, n = n, Q = Q, response_types = response_types, Y = Y, r = r)
    Z <- update_Z(Omega = Omega, response_types = response_types, Y = Y, r = r, n = n,  Q = Q)
    r <- update_R_cpp(response_types = as.numeric(factor(response_types, c("continuous", "binary", "count"))), Y = Y, eta = eta[,1:kstar,drop = FALSE], Lambda = Lambda[1:kstar,,drop = FALSE], U = U, Q = Q, n = n, r = r)
    sigma2 <- update_sigma2(a_sigma = a_sigma, b_sigma = b_sigma, n = n, U = U, Q = Q)
    eta[,1:kstar] <- scale(update_eta_cpp2(Omega = Omega, Lambda = Lambda[1:kstar,,drop = FALSE], Z = Z, p = kstar, U = U, Sigma = Sigma[1:kstar,1:kstar,drop = FALSE], B = beta[,1:kstar,drop = FALSE], n = n, X = X, rand_mat = MASS::mvrnorm(n, mu = rep(0, kstar), Sigma = diag(kstar))), scale = scale.eta)
    a1 <- update_a1(a1_cur=a1,delta1=delta1)
    if(kstar > 1) {
      a2 <- update_a2(a2_cur=a2,deltah=deltah[1:(kstar-1),drop = FALSE])

    }
    for(j in 1:Q){
      for(h in 1:kstar){
        phi[h,j] <-
          update_phijh(v=v,tauh=cumprod(c(delta1,deltah[1:(kstar-1)]))[h],
                       lambdajh=Lambda[h,j])
      }
    }

    if(kstar == 1) {
      delta1 <-
        update_delta1(a1=a1,Q = Q, p = kstar,
                      tau1l=cumprod(c(delta1))/delta1,
                      phi=phi[1:kstar,,drop = FALSE],Lambda=Lambda[1:kstar,,drop = FALSE])
    } else {
      delta1 <-
        update_delta1(a1=a1,Q = Q, p = kstar,
                      tau1l=cumprod(c(delta1,deltah[1:(kstar-1),drop = FALSE]))/delta1,
                      phi=phi[1:kstar,,drop = FALSE],Lambda=Lambda[1:kstar,,drop = FALSE])
    }

    if(kstar > 1) {
      deltatmp <- numeric(kstar-1)
      for(h in 2:kstar) {
        deltatmp[h-1] <-
          update_deltah(a2 = a2, Q = Q, p = kstar, h = h, delta = c(delta1, deltah),
                        phi=phi[1:kstar,,drop = FALSE],Lambda=Lambda[1:kstar,,drop = FALSE])
      }
      deltah[1:(kstar-1)] <- deltatmp
    }

    Lambda[1:kstar,] <- update_lambda(Omega = Omega, eta = eta[,1:kstar,drop = FALSE], Z = Z, U = U, phi = phi, delta1 = delta1, deltah = deltah, p = kstar, Q = Q)

    effb[,1:kstar] <- b_update(X = X, Y = eta[,1:kstar,drop = FALSE], Z = Z2[,1:kstar,drop = FALSE], beta = beta[,1:kstar,drop = FALSE], bsigma = bsigma, Sigma = Sigma[1:kstar,1:kstar,drop = FALSE], XtX = XtX, XtY = eigenMapMatMult(t(X), eta[,1:kstar,drop = FALSE]))
    if(s2.fit) {
      bsigma <- bsigma_update(X = X, Y = eta[,1:kstar,drop = FALSE], b = effb[,1:kstar,drop = FALSE], Sigma = Sigma[1:kstar,1:kstar,drop = FALSE], t1 = 2, t2 = 0.1, Z = Z2, beta = beta[,1:kstar,drop = FALSE])
    } else {
      bsigma <- 1
    }

    pi0 <- gpi_update(alpha)
    chi <- chi_update(gamma = gamma, AMAT = AMAT, a = a)
    a <- a_update(AMAT, chi)
    pi1 <- spi_update(AMAT, a)
    pi2 <- tpi_update(omega[1:(xp*kstar),drop = FALSE])
    alpha <- alpha_update(X = X, Y = eta[,1:kstar,drop = FALSE], b = effb[,1:kstar,drop = FALSE], group = group, pi0 = pi0, alpha = alpha, gamma = gamma, omega = omega[1:(xp*kstar),drop = FALSE], Sigma = Sigma[1:kstar,1:kstar,drop = FALSE])
    gamma <- gamma_updateC2(X = X, Y = eta[,1:kstar,drop = FALSE], b = effb[,1:kstar,drop = FALSE],
                                               group = group, pi1 = pi1,
                                               amat = matrix(rep(rep(alpha, times = table(group)), kstar), nrow = xp, ncol = kstar, byrow = FALSE),
                                               gmat = matrix(rep(gamma, kstar), nrow = xp, kstar, byrow = FALSE),
                                               omat = matrix(omega[1:(xp*kstar),drop = FALSE], ncol = kstar, byrow = FALSE),
                                               Sigma = solve(Sigma[1:kstar,1:kstar,drop = FALSE]), ncolY = kstar, ncolX = xp)
    omega[1:(xp*kstar)] <- omega_updateC2(X = X, Y = eta[,1:kstar,drop = FALSE], b = effb[,1:kstar,drop = FALSE],
                                                                         group = group, pi2 = pi2,
                                                                         amat = matrix(rep(rep(alpha, times = table(group)), kstar), nrow = xp, ncol = kstar, byrow = FALSE),
                                                                         gmat = matrix(rep(gamma, kstar), nrow = xp, kstar, byrow = FALSE),
                                                                         omat =  matrix(omega[1:(xp*kstar),drop = FALSE], ncol = kstar, byrow = FALSE), Sigma = solve(Sigma[1:kstar,1:kstar,drop = FALSE]), ncolY = kstar, ncolX = xp)
    Z2[,1:kstar] <- Z_update(X = X, Y = eta[,1:kstar,drop = FALSE], group = group, alpha = alpha, gamma = gamma, omega = omega[1:(xp*kstar),drop = FALSE])
    beta[,1:kstar] <- beta_update(b = effb[,1:kstar,drop = FALSE], Z = Z2[,1:kstar,drop = FALSE])
    if(Sigma.fit) {
      Sigma[1:kstar,1:kstar] <- Sigma_update(X = X, Y = eta[,1:kstar,drop = FALSE], beta = beta[,1:kstar,drop = FALSE], b = effb[,1:kstar,drop = FALSE], bsigma = bsigma, Z = Z2[,1:kstar,drop = FALSE])
    } else {
      Sigma[1:kstar,1:kstar] <-  diag(kstar)
    }



    if(b %in% post_inds) {
      r.store[counter,] <- r
      sigma2.store[counter,] <- sigma2
      a1.store[counter] <- a1
      if(kstar > 1) {
        a2.store[counter] <- a2
      }
      for(j in 1:Q){
        for(h in 1:kstar){
        }
      }


      if(kstar > 1) {
      }

      Lambda.store[counter,1:kstar,] <- Lambda[1:kstar,]

      if(s2.fit) {
        bsigma.store[counter] <- bsigma
      } else {
        bsigma.store[counter] <- 1
      }
      pi0.store[counter] <- pi0
      a.store[counter,] <- a
      pi1.store[counter,] <- pi1
      pi2.store[counter] <- pi2
      beta.store[counter,,1:kstar] <- beta[,1:kstar]
      if(Sigma.fit) {
      } else {
      }
    }

    lambda.max <- apply(Lambda,1,function(x) max(abs(x)))

    ut <- runif(1,0,1)
    pt <- exp(-1-0.0005*b)
    mt <- sum(lambda.max[1:kstar,drop = FALSE] <= epsilon)
    if(is.null(k_fix)) {
      if(ut <= pt & b > 20) {
        if(mt == 0) {
          kstar <- kstar + 1
          if(kstar > Q) {
            kstar <- Q
            deltah[is.na(deltah)] <- rgamma(sum(is.na(deltah)), a2.init, 1)
          } else {
            deltah[is.na(deltah)] <- rgamma(sum(is.na(deltah)), a2.init, 1)
            Sigma2 <- matrix(NA, nrow = kstar, ncol = kstar)
            Sigma2[1:(kstar-1),1:(kstar-1)] <- Sigma[1:(kstar-1),1:(kstar-1),drop=FALSE]
            Sigma <- Sigma2
            Sigma[kstar,] <- 0
            Sigma[,kstar] <- 0
            Sigma[kstar,kstar] <- 1
            eta[,kstar] <- rnorm(n = nrow(eta), 0, Sigma[kstar,kstar,drop = FALSE])
            phi[kstar,] <- rgamma(n = ncol(phi), shape = v/2, rate = v/2)
            Lambda[kstar,] <- rnorm(n = Q, mean = 0, sd = sqrt(1/(phi[kstar,]*prod(c(delta1,deltah[1:(kstar-1)])))))
          }

        } else {
          kstar <- kstar - mt
          deltah[is.na(deltah)] <- rgamma(sum(is.na(deltah)), a2.init, 1)
        }
        if(kstar < 1) {
          kstar <- 1
          deltah[is.na(deltah)] <- rgamma(sum(is.na(deltah)), a2.init, 1)
        }
      } else {
        deltah[is.na(deltah)] <- rgamma(sum(is.na(deltah)), a2.init, 1)
      }
    } else {
      kstar <- k_fix
    }

    if(b %in% post_inds) {
      kstar.store[counter] <- kstar
      counter <- counter + 1
    }

    print(kstar)
    print(mean(kstar.store, na.rm = TRUE))
    end <- proc.time()
    print(end - start)
  }
  post_inds <- seq(1+burnin,B,by=thin)
  return(list(lambda = Lambda.store[,,,drop = FALSE], kstar = kstar.store[,drop = FALSE],
              beta = beta.store[,,,drop = FALSE], a = a.store[,,drop = FALSE],
              sigma2 = sigma2.store[,,drop = FALSE], bsigma = bsigma.store[,drop = FALSE], r = r.store[,,drop = FALSE],
              a2 = a2.store, a1 = a1.store))
}

twas_pip <- function(samps, kstar) {
  if(kstar == 1) {
    pip <- apply(abs(samps$beta[,,1]), 2, FUN = function(x) mean(x>0))
  } else {
    indsums <- apply(abs(samps$beta[,,1:kstar]), 1, rowSums, na.rm = TRUE)
    pip <- apply(indsums, 1, FUN = function(x) mean(x>0, na.rm = TRUE))
  }
  return(pip)
}
