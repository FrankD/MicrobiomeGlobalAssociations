#############################################################
###           Author: Jun Chen												    ###
###           Date: 2011_09_15                            ###
###Variable selection for Dirichlet-Multinomial regression###
#############################################################

require(dirmult)
require(MASS)

#############################################################
Loglik <- function(Y, X, b,  model) {
  # Compute the log likelihood, constant part discarded
  # The likelihood is scaled. Be careful when computing AIC BIC
  
  g <- exp(X %*% t(b))	# n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  if (model == "dirmult") {
    res <- 	sum(lgamma(gs) - lgamma(ys+gs) +
                  rowSums(lgamma(Y+g) - lgamma(g))) 
  }
  if (model == "mult") {
    res <- 	sum(rowSums(Y * log(g)) - ys * log(gs))
  }
  if (model == "dir") {
    res <- 	sum(lgamma(gs) + rowSums((g-1) * log(Y) - lgamma(g)))
  }
  
  res / nrow(X)
}

ObjFunc <- function(bi, i, b0, Y, X, c1=0, c2=0, model){
  # Objective function: -Loglik + group lasso + lasso for ith group
  
  b0[, i] <- bi
  -Loglik(Y, X, b0, model) + c1*sqrt(sum(bi^2)) + c2*sum(abs(bi))
}


S <- function(Y, X, b, model) {
  # Compute the Score function at b
  
  S <- 0
  g <- exp(X %*% t(b))	# n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  
  if (model == "dirmult") {
    S <-  t((digamma(gs) - digamma(ys+gs) +
               digamma(Y+g) - digamma(g)) * g) %*% X
  }
  
  if (model == "mult") {
    S <- t((Y / g - ys / gs) * g) %*% X
  }
  
  if (model == "dir") {
    S <-  t((digamma(gs)  - digamma(g) + log(Y)) * g) %*% X
  }
  
  S / nrow(X)
}

H <- function(Y, X, b, model){
  #	Compute the diagonal of the hessian matrix at b
  
  H <- 0
  g <- exp(X %*% t(b))	# n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  if (model == "dirmult") {
    H <- t((trigamma(gs) - trigamma(ys+gs) +
              trigamma(Y+g) - trigamma(g)) * g^2 +
             (digamma(gs) - digamma(ys+gs) +
                digamma(Y+g) - digamma(g)) * g) %*% X^2
  }
  if (model == "mult") {
    H <- t((-Y / g^2 + ys / gs^2) * g^2  +
             (Y / g - ys / gs) * g) %*% X^2
  }
  if (model == "dir") {
    H <- t((trigamma(gs) - trigamma(g)) * g^2  +
             (digamma(gs)  - digamma(g) + log(Y)) * g) %*% X^2
  }
  #	Divided by the sample size
  H / nrow(X)
}
#############################################################
#	Major function

DirmultGrp <- function(Y, X, b0, lamda=0, cc=0, model="dirmult",
                       alpha0   = 1,
                       delta    = 0.5,
                       sigma    = 0.1,
                       cstar    = 0.001,
                       tol     = 1e-4,
                       iter.lim = 666) {
  # Minimize objective function: -Loglik +  group lasso + lasso
  # 	
  # Args:
  #		Y: species count matrix, row - n samples, column - q species
  #		X: design matrix with intercept, n * p
  #		b0: the initial coeffcient matrix, eg, MLE without covariates
  #		lamda: tuning parameter that controls the overall sparsity;
  #		cc: tunign parameter that controls the proportion of group lasso
  #				in the sparse group lasso penalty
  #		model: 
  #				"dirmult" - Dirichlet-Multinomial model
  #				"dir" - Dirichlet model
  #       "mult" - Multinomial model
  #   alpha0, delta, sigma, cstar: parameters used in Armijo rule
  #   tol: the numeric precision for b0;
  #
  # Returns:
  # 		List: b1, iter.no, Loglik, lamda, cc
  #					The coefficients are shrinked due to double sparse penalties. 
  #					Re-estimation is recommonded for tuning parameter selection
  
  b0 <- as.matrix(b0)		
  n <- nrow(X)
  p <- ncol(b0)
  q <- nrow(b0)
  c1 <- sqrt(q) * lamda * cc # so c1, c2 are on the same level
  c2 <- lamda * (1-cc)
  if (ncol(Y) != q || ncol(X) != p) stop("Dimension does not match!\n")
  if (sum(X[, 1] != 1) != 0) {
    warning("No intercepts in X! Intercepts added\n")
    X <- cbind(rep(1, nrow(X)), X)
  }
  if (nrow(Y) != nrow(X)) stop("Y and X have different sample sizes!\n")	
  if (model == "dir") {	
    # Y are not proportions
    if (max(Y) > 1) {
      Y <- (Y + 0.5) / rowSums(Y)
    } else {
      if (min(Y) == 0) {
        Y[Y==0] <- min(Y[Y!=0])/10
      }
    }
  }
  
  b1 <- b0
  iter <- 0
  
  repeat {
    b0 <- b1
    h <- H(Y, X, b0, model)
    s <- S(Y, X, b0, model)
    
    #	Minimize with respect to the intercepts
    h1 <- -max(-h[, 1], cstar)
    d1 <- -s[, 1] / h1
    if (sum(d1 != 0) != 0){
      # inexact search using the Armijo rule
      Delta <- - sum(d1 * s[, 1])
      alpha <- alpha0
      alpha.d1 <- alpha * d1
      
      while (ObjFunc(b1[, 1] + alpha.d1, 1, b1, Y, X, model=model) - 
             ObjFunc(b1[, 1], 1, b1, Y, X, model=model) > alpha*sigma*Delta) {
        alpha <- alpha * delta
        alpha.d1 <- alpha * d1
        # This step is used to eliminate those small non-zero coefficients
        if (max(abs(alpha.d1)) < 0.1*tol){alpha.d1 <- 0; break}
      }
      b1[, 1] <- b1[, 1] + alpha.d1
    }
    
    # Minimize with respect to the penalized groups
    di <- numeric(q)
    
    if (p != 1){
      for (i in 2:p) {
        hi <- -max(-h[, i], cstar)
        # First determine which elements of ith group should be zeros 
        # (Lasso in action)
        diff.v <- s[, i]- hi * b1[, i]
        ind <- (abs(diff.v) <= c2)
        di[ind] <-  -b1[ind, i]	
        if (sum(ind) != q) {
          ind <- !ind				
          si.temp <- s[ind, i] - c2*(sign(s[, i] - hi*b1[, i])[ind])				
          # Group Lasso in action for those not shrunk to zeros by Lasso
          diff.v <- si.temp - hi*b1[ind, i] #	dimension changed
          diff.n <- sqrt(sum(diff.v^2))
          if (diff.n <= c1){				
            di[ind]  <- -b1[ind, i]			
          } else {	
            di[ind] <- -(si.temp - c1*diff.v/diff.n) / hi
          }
        }
        
        if (sum(di != 0) != 0){
          # When adding the lasso penalty, this delta should be changed accordingly
          Delta <- - sum(di * s[, i]) + 
            c1 * (sqrt(sum((b1[, i] + di)^2)) - sqrt(sum(b1[, i]^2))) +
            c2 * sum(abs(b1[, i] + di) - abs(b1[, i]))
          alpha <- alpha0
          alpha.di <- alpha * di
          while (ObjFunc(b1[, i] + alpha.di, i, b1, Y, X, c1, c2, model) - 
                 ObjFunc(b1[, i], i, b1, Y, X, c1, c2, model) > alpha*sigma*Delta) {
            alpha <- alpha * delta
            alpha.di<- alpha * di
            # This step is used to eliminate those small non-zero coefficients
            if (max(abs(alpha.di)) < 0.1*tol){ alpha.di <- 0; break }
          }
          b1[, i] <- b1[, i] + alpha.di
        }
      }
    }
    
    iter <- iter + 1
    if (iter > iter.lim) {
      stop("Iter exceeds limit!")
    }
    if (max(abs(b1-b0)) <= tol) break
  }
  
  lik <- n * Loglik(Y, X, b1, model)
  return(list(b=b1, lamda=lamda, cc=cc, iter=iter, loglik=lik))
}

TuningMax <- function(Y, X, incpt, cc=0.00, model, st=NULL) {
  # Determine the value of tuning parameter lamda to achieve maximum sparsity
  # 	
  # Args:
  #		Y: species count matrix, row - n samples, column - q species
  #		X: design matrix with intercept, n * p
  #   incpt: MLE of the intercept without any covariate
  #		cc: the proportion of group lasso in the sparse group lasso penalty
  #		model: 
  #				"dirmult" - Dirichlet-Multinomial model
  #				"dir" - Dirichlet model
  #       "mult" - Multinomial model	
  #		st: the starting maximum value of lambda. If it is NULL,
  #				the default value is used
  #
  # Returns:
  # 		lamda.max	
  if (is.null(st)) {
    if (model == "dirmult" | model == "dir") {
      st <- 4
    } else {
      st <- 50
    }
  }
  
  cat("Determine max lamda...\n")
  f <- function(lamda) {
    flag <- 0
    c1 <- lamda * cc * sqrt(q)
    c2 <- lamda * (1-cc)
    for (i in 2:p) {
      diff.v <- s[, i]
      ind <- (abs(diff.v) <= c2)
      if (sum(ind) == q) next
      ind <- !ind
      si.temp <- s[ind, i] - c2*(sign(s[, i])[ind])
      diff.v <- si.temp
      diff.n <- sqrt(sum(diff.v^2))
      if (diff.n > c1) {
        flag <- 1
        break
      }
    }
    flag
  }
  p <- ncol(X)
  q <- ncol(Y)
  b0 <- matrix(0, q, p)
  b0[, 1] <- incpt
  s <- S(Y, X, b0, model)
  lamda <- st
  flag <- 0
  while (flag == 0) {
    lamda <- lamda * 0.96
    flag <- f(lamda)
  }
  lamda / 0.96	
}

DirmultGrpGrid <- function(Y, X, n.grid=50, lamda.max0=NULL, cc=0.20, nz=50,
                           model="dirmult", initscalar=10, ...) {
  # Perform DirmultGrp over a grid of parameter values
  # 	
  # Args:
  #		Y: species count matrix, row - n samples, column - q species
  #		X: design matrix with intercept, n * p
  #   n.grid: the number of grid points
  #		lamda.max0: the maximum value of lamda that controls the overall sparsity;
  #			 	it could be a vector, with different value for each different 'cc'.
  #			 	If it is not supplied, the maximum value will be computed automatically
  #		cc: the proportion of group lasso in the sparse group lasso penalty
  #		nz: stop if the number of nonzeros reaches this number 
  #		model: 
  #				"dirmult" - Dirichlet-Multinomial model
  #				"dir" - Dirichlet model
  #       "mult" - Multinomial model
  #		initiscalar: the parameter for "dirmult" function. 
  #				In case the "dirmult" function fails, try to use a smaller value
  #   ...: parameters pass to DirmultGrp
  #
  # Returns:
  # 		res: list of objects from DirmultGrp
  #					The coefficients are shrinked due to double sparse penalties. 
  #					Re-estimation is recommonded for tuning parameter selection
  res <- list()
  p <- ncol(X)
  q <- ncol(Y)
  
  if (ncol(Y) != q || ncol(X) != p) stop("Dimension not match!\n")
  if (sum(X[, 1] != 1) != 0) {
    warning("No intercepts in X! Intercepts added\n")
    X <- cbind(rep(1, nrow(X)), X)
  }
  if (nrow(Y) != nrow(X)) stop("Y and X have different sample size!\n")
  if (sum(round(Y) != Y)) stop("Y should contain counts!\n")
  if (model == "dir") {	
    Y <- (Y + 0.5) / rowSums(Y)
  }	
  if (model == "dirmult") {
    cat("Sparse Dirichlet-Multinomial Regression\n")
  }
  if (model == "dir") {
    cat("Sparse Dirichlet Regression\n")
  }
  if (model == "mult") {
    cat("Sparse Multinomial Regression\n")
  }
  cat("Initial MLE of the coefficient matrix ...\n")
  b00 <- matrix(0, q, p)
  if (model == "dirmult" ) {
    b00[, 1] <- log(dirmult(Y, initscalar=initscalar, trace=F)$gamma)
  } 
  if (model == "mult" | model == "dir") {
    f1 <- function(b) {
      b00[, 1] <- b
      -Loglik(Y, X, b00, model)
    }	
    b00[, 1] <- nlm(f1, rep(0, q))$estimate
  } 
  
  cat("Grid search:\n")
  ct <- 0
  for (j in 1:length(cc)) {
    b0 <- b00
    if(sum(is.null(lamda.max0)) != 0) {
      lamda.max <- TuningMax(Y, X, incpt=b00[, 1], cc=cc[j], model=model)
      lamda.range <- lamda.max * 0.96^(0:n.grid) 
    } else if(length(lamda.max0) == 1) {
      lamda.range <- lamda.max0 * 0.96^(0:n.grid) 
    } else {
      lamda.range <- lamda.max0[j] * 0.96^(0:n.grid) 
    }
    
    for (i in 1:length(lamda.range)) {
      lamda <- lamda.range[i]
      e <- try(obj <- DirmultGrp(Y, X, b0, lamda, cc[j], model=model, ...))
      if (class(e) == "try-error") stop("DirmultGrp error!")
      ct <- ct+1
      res[[ct]] <- obj
      b0 <- obj$b
      cat(".")
      # To restrict the number of paramter 
      if(sum(b0[, -1] != 0) >= nz) break
    }
    cat("\n")
  }
  cat("Finished!\n")
  res
}
#############################################################
# Example

DirmultSim <- function(n=100, nr=1000, s=2, rho=0.4,
                       p=100, q=40, p.g=4, q.g=4, f=0.80, theta0=0.025){
  # Simulation strategy 
  # 	
  # Args:
  #		n: the number of samples
  #		nr: number of reads for each sample
  #   s: scenario, 2 - exponential growth, 3 - linear growth
  #   rho: the correlation between covariates
  #   p: number of covariates excluding the intercept
  #		q: number of species
  #   p.g: number of relevant nutrients
  #		q.g: number of relevant species
  #   f: controls the signal
  #		theta0: the dispersion parameter
  #
  # Returns:
  #		Y: count matrix, row - n samples, column - q species
  #		X: design matrix, with intercepts, n * (p+1)
  #		b: simulated coefficients q * (p+1)
  #		pi, theta: used by dirmult
  theta <- theta0
  cat("Generating data ...\n")
  
  ct <- sample(nr:(2*nr), n, rep=T)
  
  Sigma <- matrix(1, p, p)
  Sigma <- rho^abs(row(Sigma)-col(Sigma))
  
  X <- cbind(rep(1, n), scale(mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)))
  Y <- matrix(0, n, q)
  b <- matrix(0, q, p)
  pi <- matrix(0, n, q)
  
  # coefficient with signs alternating
  st <- 0
  if (q.g != 1) {
    coef <- seq(0.6, 0.9, len=q.g) * c(1, -1)	
  } else {
    coef <- (0.6 + 0.9) / 2
  }
  
  coef.g <- seq(1.0, 1.0, len=p.g)
  
  for (i in 1:p.g) {
    # overlap two species
    # No overlaps for simplicity
    # q may be small so enable wrap-up
    b[(st:(st+q.g-1))%%q+1, 3*i-2] <- coef.g[i] * coef[((i-1):(i+q.g-2))%%q.g+1]
    st <- st+1  
  }
  
  if (s==1) {
    gs <- (1-theta) / theta
    # base proportion max diff 100 fold
    icpt <- runif(q, 0.02, 2)
    # The theta for each sample is different, but let them close to supplied value
    icpt <- gs * icpt / sum(icpt) 
    icpt <- log(icpt)		
  } 
  if (s==2) {
    # exponential growth
    # base proportion max diff 100 fold
    icpt <- runif(q, -2.3, 2.3)
  }
  if (s==3) {
    # linear growth
    # base proportion max diff 100 fold
    icpt <- runif(q, 0.02, 2)
    g.m <- X[, -1] %*% t(f*b) 
    adj <- apply(g.m, 2, min)
    icpt[adj < 0] <- icpt[adj < 0] - adj[adj < 0] + 0.02
  }
  
  b <- cbind(icpt, f*b)
  
  for (i in 1:n){
    if (s==1 || s==2){
      # Exponential function
      g <- as.vector(exp(b %*% X[i, ]))
    } else if (s==3) {
      # Linear function
      g <- as.vector(b %*% X[i, ])
    } else if (s==4) {
      # atan function
      X.temp <- c(1, 0.25 * X[i, -1])
      g <- as.vector(atan(b %*% X.temp))
    }
    pi[i, ] <- g / sum(g)
    if (s==1) {
      # Exactly the same as we model the data
      theta <- 1 / (sum(g)+1)
    }
    if (theta == 0){
      Y[i, ] <- rmultinom(1, ct[i], pi[i, ])[, 1]
    } else {
      Y[i, ] <- simPop(J=1, n=ct[i], pi=pi[i, ], theta=theta)$data[1, ]
    }
  }
  cat("Finished!\n")
  return(list(X=X, Y=Y, b=b, theta=theta0, pi=pi, s=s, p=p, q=q, nr=nr,
              p.g=p.g, q.g=q.g, rho=rho, f=f))
}

sim.obj <- DirmultSim()
Y <- sim.obj$Y
X <- sim.obj$X
dm.obj <- DirmultGrpGrid(Y, X, model="dirmult")
dir.obj <- DirmultGrpGrid(Y, X, model="dir")
mult.obj <- DirmultGrpGrid(Y, X, model="mult")