# Log-Likelihood for missing data -----------------------------------------

obs.llik.binom.std <- 
  function(par, # initial parameters vector
           J, # number of control units
           y, # response to question
           treat, # treatment status
           x, # covariates matrix
           const = FALSE) # include constant term
  {
    k <- ncol(x) # k - number of covariates
    if (const) {
      coef.h0 <- coef.h1 <- par[1:k] 
      coef.g <- par[(k+1):(2*k)]
    } else {
      coef.h0 <- par[1:k]
      coef.h1 <- par[(k+1):(2*k)]
      coef.g <- par[(2*k+1):(3*k)]
    }
    
    h0X <- logistic(x %*% coef.h0)
    h1X <- logistic(x %*% coef.h1)
    gX <- logistic(x %*% coef.g)
    
    ind10 <- ((treat == 1) & (y == 0))
    ind1J1 <- ((treat == 1) & (y == (J+1)))
    ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
    
    if (sum(ind10) > 0) {
      p10 <- sum(log(1-gX[ind10]) + 
                   dbinom(x = 0, size = J, prob = h0X[ind10], log = TRUE))
    } else {
      p10 <- 0
    }
    if (sum(ind1J1) > 0) {
      p1J1 <- sum(log(gX[ind1J1]) + dbinom(J, size = J, prob = h1X[ind1J1], log = TRUE))
    } else {
      p1J1 <- 0
    }
    if (sum(ind1y) > 0) {
      p1y <- sum(log(gX[ind1y]*dbinom(y[ind1y]-1, size = J, prob = h1X[ind1y],
                                      log = FALSE) + (1-gX[ind1y])*
                       dbinom(y[ind1y], size = J, prob = h0X[ind1y], log = FALSE)))
    } else {
      p1y <- 0
    }
    
    if (sum(treat == 0) > 0) {
      p0y <- sum(log(gX[!treat]*dbinom(y[!treat], size = J, prob = h1X[!treat], log = FALSE) +
                       (1-gX[!treat]) * 
                       dbinom(y[!treat], size = J, prob = h0X[!treat], log = FALSE)))
    } else {
      p0y <- 0
    }
    
    return(p10+p1J1+p1y+p0y)
  }

# E-step: Conditional expectation for binomial ----------------------------

Estep.binom.std <- 
  function(par, # initial parameters vector
           J, # number of control units
           y, # response to question
           treat, # treatment status
           x, # covariates matrix
           const = FALSE) # include constant term
  {
    k <- ncol(x)
    if (const) {
      coef.h0 <- coef.h1 <- par[1:k]
      coef.g <- par[(k+1):(2*k)]
    } else {
      coef.h0 <- par[1:k]
      coef.h1 <- par[(k+1):(2*k)]
      coef.g <- par[(2*k+1):(3*k)]
    }
    
    h0X <- logistic(x %*% coef.h0)
    h1X <- logistic(x %*% coef.h1)
    gX <- logistic(x %*% coef.g)
    
    ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
    w <- rep(NA, length(y))
    w[ind] <- gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) /
      (gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) +
         (1-gX[ind])*dbinom(y[ind], size = J, prob = h0X[ind], log = FALSE))
    w[(treat == 1) & (y == 0)] <- 0
    w[(treat == 1) & (y == (J+1))] <- 1
    
    return(w)
  }

# M-step 1: Weighted MLE for logistic regression --------------------------

wlogit.fit.std <- 
  function(y, # response to question
           treat, # treatment status
           x, # covariates matrix
           w, # weights from E-step
           par = NULL) # initial parameters vector for linear prediction in glm()
  {
    yrep <- rep(c(1,0), each = length(y))
    xrep <- rbind(x, x)
    wrep <- c(w, 1-w)
    return(glm(cbind(yrep, 1-yrep) ~ xrep - 1, weights = wrep, family = binomial(logit), 
               start = par))
  }

# M-step 2: Weighted MLE for binomial regression --------------------------

wbinom.fit.std <- 
  function(J, # number of control units
           y, # response to question
           treat, # treatment status
           x, # covariates matrix
           w, # weights from E-step
           par0 = NULL, # initial parameters vector
           par1 = NULL) # new parameters vector
  {
    Not0 <- ((treat == 1) & (y == (J+1)))
    y0 <- y[!Not0]
    x0 <- x[!Not0,]
    w0 <- 1-w[!Not0]
    fit0 <- glm(cbind(y0, J-y0) ~ x0, family = binomial(logit), weights = w0, start = par0)
    
    Not1 <- ((treat == 1) & (y == 0))
    y1 <- y
    y1[treat == 1] <- y1[treat == 1] - 1
    y1 <- y1[!Not1]
    x1 <- x[!Not1,]
    w1 <- w[!Not1]
    fit1 <- glm(cbind(y1, J-y1) ~ x1, family = binomial(logit), weights = w1, start = par1)
    
    return(list(fit1 = fit1, fit0 = fit0))
  }


# EM Algorithm ------------------------------------------------------------

## input parameters of the function:
## coef.control.start, coef.treat.start, 
## formula, maxIter, verbose,
## J - Number of non-sensitive (control) survey items.
## nPar <- ncol(x.all)
## x.all <- model.matrix.default(attr(mf, "terms"), mf)
## y.all <- model.response(mf)
## t: 
## if (design == "standard") 
##    t <- data[na.x == 0 & na.y == 0 & na.w == 0, paste(treat)]
## if (design == "modified") 
##    t <- as.numeric(is.na(y.all[na.x == 0 & na.y == 0 & na.w == 0, J + 1]) == FALSE)
## n:
## n <- nrow(x.treatment) + nrow(x.control)



## set some starting values
par <- c(coef.control.start, coef.treat.start)

## start off with an infinitely negative log likelihood
pllik.const <- -Inf

## calculate the log likelihood at the starting values
llik.const <- 
  obs.llik.binom.std(par, 
                     J = J, 
                     y = y.all, 
                     treat = t,
                     x = x.all, 
                     const = TRUE)

## set a counter to zero, which you will iterate
## (this just allows to you stop after a maximum # of iterations if you want
## e.g. 10k

counter <- 0

## begin the while loop, which goes until the log likelihood it calculates
## is very very close to the log likelihood at the last iteration (the
## difference is less than 10^(-8) different

while (((llik.const - pllik.const) > 10^(-8)) & (counter < maxIter)) {
  
  ## first the E step is run each time, from the parameters at the end of the
  ## last iteration, or the starting values in the first iteration
  w <- Estep.binom.std(par, J, y.all, t, x.all, const = TRUE)

  ## then the M step is run. there are two parts to that in our case
  ## this is the first one - the fit for the treatment
  ## this logistic is weighted by the weights from the E step
  lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar+1):(nPar*2)])

  ## a bunch of stuff for the package you can ignore
  y.var <- as.character(formula)[[2]]
  
  x.vars <- strsplit(gsub(" ", "",as.character(formula)[[3]]),
                     split="+", fixed=T)[[1]]

  ## now create the fake dataset that replicates the data a few times
  ## so that we can set up the weights from the E step to weight
  ## at each missing data possibility
  data.all <- as.data.frame(cbind(y.all, x.all))
  names(data.all)[1] <- y.var
  
  dtmp <- rbind(cbind(data.all, w, t), cbind(data.all, w, t)[t==1, ])
  
  dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] <-
    dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] - 1
  
  dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))] <-
    1 - dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))]
  
  dtmp$w[dtmp$t == 0] <- 1
  
  dtmp <- dtmp[dtmp$w > 0, ]

  ## now run the second part of the M step in our case, again weighted
  ## using the multiplied dataset
  fit <- glm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ",
                              paste(x.vars, collapse=" + "))),
             family = binomial(logit), weights = dtmp$w,
             start = par[1:(nPar)], data = dtmp)
  
  par <- c(coef(fit), coef(lfit))
  
  pllik.const <- llik.const
  
  if(verbose==T)
    cat(paste(counter, round(llik.const, 4), "\n"))

  ## calculate the new log likelihood with the parameters you just
  ## estimated in the M step
  llik.const <- obs.llik.binom.std(par, J = J, y = y.all, treat = t,
                                   x = x.all, const = TRUE)
  
  counter <- counter + 1

  ## stop for obvious reasons - this means your code is wrong
  ## or the data is really fucked up
  if (llik.const < pllik.const)
    stop("log-likelihood is not monotonically increasing.")
  
  if(counter == (maxIter-1))
    warning("number of iterations exceeded maximum in ML")
  
} ## end of while loop

## you also want standard errors, so do that by running optim
## once on the final parameter estimates from EM
MLEfit <- optim(par, obs.llik.binom.std, method = "BFGS", J = J,
                y = y.all, treat = t,  x = x.all, const = TRUE,
                hessian = TRUE, control = list(maxit = 0))

vcov.mle <- solve(-MLEfit$hessian)
se.mle <- sqrt(diag(vcov.mle))


