#' @title bord.prior Function 
#'
#' @description This function is a code to sample from the prior distribution for the 
#' cumulative logistic distribution. The important assumption is proportional odds.
#'
#' @param formula as in \code{\link{lm}}
#' @param data as in \code{\link{lm}}
#' @param family as in \code{\link{glm}}. Only "multinomial" is supported
#' @param nlev number of categories in the toxicity profile. If NULL this will be counted
#' from the data
#' @param nsample number of samples from prior distribution. This is only required if 
#' Prior.Only = TRUE
#' @param prior.mean prior mean of the regression parameter
#' @param prior.var prior variance of the regression parameter.
#' @param newx new covariates for prediction. If NULL the observed data will be used
#' @param alpha.min the lowest point of the intercept parameter space
#' @param alpha.max the highest point of the intercept parameter space
#' @param ci significance level to calculate Bayesian credible interval. Default is 0.95
#' @param Pcutoffs target toxicity interval -- a vector containg two probability elements -- first 
#' is the boundary point between underdose and target toxicity region, second is the boundary 
#' between target toxicity and overdose region
#' @param Pcutoffs3 same as above for grade 3
#' @param Pcutoffs4 same as above for grade 4
#' @param DLT.cutoff adverse event grade cut off for DLT, must be 3 or 4.
#' If 3 then DLT >= grade 3, if 4 then DLT >= grade 4
#' @param seed Seed for R for reproducibility
#'
#'
#' @details For single-agent data, the base model, is a standard cumulative logistic 
#' proportional odds model in log-dose. We place the Normal(prior.mean, prior.var) distribution on the 
#' regression paramter \eqn{\beta}
#' 
#' @references Albert, James A. and Chib, Siddhartha. (1993) Bayesian Analysis of Binary
#' and Polychotomous Response Data. Journal of the American Statistical Association,
#' 88(422), 669-679.
#' 
#' Sorensen, D. A., Anderson, S., Gianola, D., and Korsgaard, I. (1995) Bayesian 
#' inference in thereshold models using Gibbs sampling. Genetics Selection Evolution,
#' 27(3), 229-249.
#' 
#' Polson, Nicholas G., Scott, James G., and Windle, Jesse. (2013) Bayesian Inference 
#' for Logistic Models Using Polya-Gamma Latent Variables, Journal of the American
#' Statistical Association, 108:504, 1339-1349.
#' 
#' Montesinos-Lopez, O. A., Montesinos-Lopez, A., Crossa, J., Burgueno, J., 
#' & Eskridge, K. (2015). Genomic-enabled prediction of ordinal data with Bayesian 
#' logistic ordinal regression. G3: Genes, Genomes, Genetics, 5(10), 2113-2126.
#' 
#' 
#' @seealso \code{\link{bordCRM}}
#' 
#' @return \item{beta.post}{Posterior mean of \eqn{\beta}}
#' \item{beta.sd}{Posterior standard deviation of \eqn{\beta}}
#' \item{beta.median}{Posterior meadian of \eqn{\beta}}
#' \item{beta.summary}{Posterior summary of \eqn{\beta}}
#' \item{left.ci}{Lower point of the Bayesian credible interval for \eqn{\beta}}
#' \item{right.ci}{Upper point of the Bayesian credible interval for \eqn{\beta}}
#' \item{beta.samples}{Posterior samples of \eqn{\beta}}
#' \item{pTox.post}{Posterior predictive summaries of the dose levels}
#' \item{pCat.post}{Posterior predictive probabilities of doses}
#' \item{predTox.post}{Posterior predictive toxicity probilities of the future doses}
#' \item{pred.dose}{Doses for compound at which the prediction is computed. This is 
#' same as newx}
#' \item{alpha.post}{Posterior means of the intercepts}
#' \item{grade.prob.post}{Predictive probabilities of doses for each adverse event grade}
#' 
#' @examples 
#' data(aegrade)
#' 
#' fit <- bord.prior(formula    = as.factor(AEgrade) ~ Dose, 
#'                   data       = aegrade, 
#'                   nsample    = 100, 
#'                   prior.mean = 0,
#'                   prior.var  = 100, 
#'                   Pcutoffs   = c(0.16, 0.33),
#'                   Pcutoffs3  = c(0.10, 0.25),
#'                   Pcutoffs4  = c(0.05, 0.10),
#'                   DLT.cutoff = 3,
#'                   newx       = 3, 
#'                   nlev       = 5,
#'                   alpha.min  = -1000,
#'                   alpha.max  = 1000
#'                   )
#'               
#'  fit$beta.post
#'  fit$alpha.post
#'
#'
#'
#' @export


bord.prior <- function(formula, 
                       data       = NULL, 
                       family     = "multinomial",
                       nlev       = NULL,
                       nsample    = 1000,
                       prior.mean = 0,
                       prior.var  = 1,
                       alpha.min  = -1000,
                       alpha.max  = 1000,
                       newx       = NULL,
                       ci         = 0.95,
                       Pcutoffs   = NULL,
                       Pcutoffs3  = NULL,
                       Pcutoffs4  = NULL,
                       DLT.cutoff = 3,
                       seed       = 1
)
{
  call <- match.call()
  if (family != "multinomial")
    stop("families other than 'multinomial' are not currently implemented")
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "any")
  if (!is.factor(y)) 
    stop("response must be a factor")
  if(is.null(nlev))
  {
    nlev <- nlevels(y)  # number of categories 
  } 
  if (nlev <= 2L) 
    stop("response must have 3 or more levels")
  y <- unclass(y)
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) {
    X <- NULL
  } else {
    X <- model.matrix(mt, mf)
  }
  n <- length(y)
  p <- ncol(X)
  x <- X[, -1]
  DoseAdm <- sort(unique(x))
  
  
  if(is.null(newx))
  {
    d <- DoseAdm
  } else {
    d <- newx   # we want to find the toxicity probability at this dose
  }
  
  Pcutoffs1  <- c(0, Pcutoffs, 1)
  Pcutoffs31 <- c(0, Pcutoffs3, 1)
  Pcutoffs41 <- c(0, Pcutoffs4, 1)
  Ncat       <- length(Pcutoffs) + 1
  pCat       <- matrix(NA, nrow = length(d), ncol = Ncat)
  pCat3      <- matrix(NA, nrow = length(d), ncol = Ncat)  # For Grade 3 
  pCat4      <- matrix(NA, nrow = length(d), ncol = Ncat)  # For Grade 4
  grade.prob <- matrix(NA, nrow = length(d), ncol = nlev)
  
  # MCMC parameters
  effsamp <- nsample
  
  id1 <- nid1 <- id2 <- NULL  # to pass the package check
  
  z  <- rep(0, n)
  for(j in 1:nlev)
  {
    yj <- which(y == j)
    assign(paste("id", j, sep = ""), yj)
    assign(paste("nid", j, sep = ""), length(yj))
  }
  
  alpha.out      <- matrix(NA, nrow = nlev - 1, ncol = effsamp)
  beta.out       <- rep(NA, effsamp)
  pTox.out       <- matrix(NA, nrow = length(d), ncol = effsamp)
  pCat.out       <- array(NA, dim = c(length(d), Ncat, effsamp))
  pCat3.out      <- array(NA, dim = c(length(d), Ncat, effsamp))
  pCat4.out      <- array(NA, dim = c(length(d), Ncat, effsamp))
  predTox.out    <- matrix(NA, nrow = length(d), ncol = effsamp)
  grade.prob.out <- array(NA, dim = c(length(d), nlev, effsamp))
  
  beta.mean <- prior.mean
  beta.var  <- prior.var
  
  alpha.min <- alpha.min
  alpha.max <- alpha.max
  
  #######################################################
  ################ Gibbs Sampling #######################
  #######################################################
  
  # Initialization to start the chain
  alpha <- seq(from = -2, to = 2, length.out = nlev - 1)
  beta  <- 0
  
  cat("Sampling from prior \n")
  
  set.seed(seed)
  for(iter in 1:nsample)
  {
    # Sample beta
    var.beta  <- beta.var
    mean.beta <- beta.mean
    beta      <- rnorm(n = 1, mean = mean.beta, sd = sqrt(var.beta))
    
    # Sample alpha
    alpha <- alpha.min + (alpha.max - alpha.min) * sort(orderstats::order_probs(draw_size = nlev - 1, k = nlev - 1, n = nlev - 1))
    
    if(DLT.cutoff == 3)
    {
      # Probability of toxicity assuming toxicity is >= 3
      # Use brob() replacing exp() to avoid infinity and NaN
      pTox <- as.numeric(Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d)) -
                           Brobdingnag::brob(alpha[2] - beta*d)/(1 + Brobdingnag::brob(alpha[2] - beta*d))
                         + Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)) -
                           Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d))
                         + 1 - Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)))
      


      # Predictive toxicity assuming toxicity is >= 3
      # predTox <- (exp(alpha[3] - beta*d)/(1 + exp(alpha[3] - beta*d)) - exp(alpha[2] - beta*d)/(1 + exp(alpha[2] - beta*d))
      #          + exp(alpha[4] - beta*d)/(1 + exp(alpha[4] - beta*d)) - exp(alpha[3] - beta*d)/(1 + exp(alpha[3] - beta*d))
      #          + 1 - exp(alpha[4] - beta*d)/(1 + exp(alpha[4] - beta*d)))
      predTox <- as.numeric(Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d)) -
                                Brobdingnag::brob(alpha[2] - beta*d)/(1 + Brobdingnag::brob(alpha[2] - beta*d))
                              + Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)) -
                                Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d))
                              + 1 - Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)))
        
      
    } else if (DLT.cutoff == 4) {
      # Probability of toxicity assuming toxicity is >= 4
      
      pTox <- as.numeric(Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)) - 
                                 Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d))
                               + 1 - Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)))
      
      
      # Predictive toxicity assuming toxicity is >= 4
      # predTox <- (exp(alpha[4] - beta*d)/(1 + exp(alpha[4] - beta*d)) - exp(alpha[3] - beta*d)/(1 + exp(alpha[3] - beta*d))
      #          + 1 - exp(alpha[4] - beta*d)/(1 + exp(alpha[4] - beta*d)))
      predTox <- as.numeric(Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)) - 
                              Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d))
                  + 1 - Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)))
      
    }
    
    
    # Predictive toxicity assuming toxicity is = 3
    predTox3 <- as.numeric(Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d)) - 
                             Brobdingnag::brob(alpha[2] - beta*d)/(1 + Brobdingnag::brob(alpha[2] - beta*d)))
    
    
    # Predictive toxicity assuming toxicity is = 4
    predTox4 <- as.numeric(Brobdingnag::brob(alpha[4] - beta*d)/(1 + Brobdingnag::brob(alpha[4] - beta*d)) - 
                   Brobdingnag::brob(alpha[3] - beta*d)/(1 + Brobdingnag::brob(alpha[3] - beta*d)))
    
    
    
    # Interval probabilities
    for(jCat in 1:Ncat)
    {
      for(k in 1:length(d))
      {
        pCat[k, jCat] <- (Pcutoffs1[jCat] <= predTox[k]) && (predTox[k] <= Pcutoffs1[jCat + 1])
        pCat3[k, jCat] <- (Pcutoffs31[jCat] <= predTox3[k]) && (predTox3[k] <= Pcutoffs31[jCat + 1])
        pCat4[k, jCat] <- (Pcutoffs41[jCat] <= predTox4[k]) && (predTox4[k] <= Pcutoffs41[jCat + 1])
      }
      
    }
    
    # Grade probabilities
    for(ilev in 1:nlev)
    {
      if(ilev == 1)
      {
        grade.prob[, ilev] <- as.numeric(Brobdingnag::brob(alpha[ilev] - beta*d)/
                                           (1 + Brobdingnag::brob(alpha[ilev] - beta*d)))
      } else if (ilev == nlev) {
        grade.prob[, ilev] <- as.numeric(1 - Brobdingnag::brob(alpha[ilev - 1] - beta*d)/
          (1 + Brobdingnag::brob(alpha[ilev - 1] - beta*d)))
      } else {
        grade.prob[, ilev] <- as.numeric(Brobdingnag::brob(alpha[ilev] - beta*d)/(1 + Brobdingnag::brob(alpha[ilev] - beta*d)) -
          Brobdingnag::brob(alpha[ilev - 1] - beta*d)/(1 + Brobdingnag::brob(alpha[ilev - 1] - beta*d)))
      }
    }

    
    
    if(sum(is.na(alpha) > 0))
      stop("The posterior space for some of the threshold parameters becomes degenerate 
           Consider providing alpha.min and alpha.max towards -Inf and Inf respectively")
    
    if (iter %% 1000 == 0)
    {
      cat("iteration = ", iter, "\n")
    }
    
    alpha.out[, iter]        <- alpha
    beta.out[iter]           <- beta
    pTox.out[, iter]         <- pTox
    pCat.out[, , iter]       <- pCat
    pCat3.out[, , iter]      <- pCat3
    pCat4.out[, , iter]      <- pCat4
    predTox.out[, iter]      <- predTox
    grade.prob.out[, , iter] <- grade.prob
    
  }
  
  alpha.post      <- apply(alpha.out, 1, mean)
  beta.post       <- mean(beta.out)
  pCat.post       <- apply(pCat.out, c(1, 2), mean, na.rm = TRUE)
  pCat3.post      <- apply(pCat3.out, c(1, 2), mean, na.rm = TRUE)
  pCat4.post      <- apply(pCat4.out, c(1, 2), mean, na.rm = TRUE)
  grade.prob.post <- apply(grade.prob.out, c(1, 2), mean, na.rm = TRUE) 
  
  if(length(d) == 1)
  {
    pTox.post       <- summary(as.vector(pTox.out), na.rm = TRUE)
    predTox.post    <- mean(as.vector(predTox.out), na.rm = TRUE)
  } else {
    pTox.post       <- apply(pTox.out, 1, summary, na.rm = TRUE) 
    predTox.post    <- apply(predTox.out, 1, mean, na.rm = TRUE)
  }
  
  
  left        <- floor(ci * effsamp/2)
  right       <- ceiling((1 - ci/2) * effsamp)
  beta.sort   <- sort(beta.out, decreasing = F)
  left.point  <- beta.sort[left]
  right.point <- beta.sort[right]
  
  return(list(beta.post       = beta.post,
              beta.sd         = sd(beta.out),
              beta.median     = median(beta.out),
              beta.summary    = summary(beta.out),
              left.ci         = left.point,
              right.ci        = right.point,
              beta.samples    = beta.out,
              pTox.post       = pTox.post,
              pCat.post       = pCat.post,
              pCat3.post      = pCat3.post,
              pCat4.post      = pCat4.post,
              predTox.post    = predTox.post,
              pred.dose       = d,
              alpha.post      = alpha.post,
              grade.prob.post = grade.prob.post,
              grade.prob      = grade.prob.out
  ))
}