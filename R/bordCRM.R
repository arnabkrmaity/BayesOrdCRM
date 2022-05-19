#' @title bordCRM Function 
#'
#' @description This function is a code to the continuous assesment method (CRM) with ordinal 
#' endpoints in the oncology participants data. The important assumption is proportional odds.
#'
#' @param Data Dose DLT data -- a data frame with the first vector as the dose limiting toxicity 
#' response which is ordinal and the second vector is the dose
#' @param Agent Compound name
#' @param Doses Doses for compound at which the prediction is required
#' @param DoseRef Dose reference for compound
#' @param nlev number of categories in the toxicity profile. If NULL this will be counted
#' from the data
#' @param prior.mean prior mean of the regression parameter
#' @param prior.var prior variance of the regression parameter.
#' @param Prior.Only either TRUE or FALSE. If TRUE, Data must be equal to NULL meaning that there 
#' is no data and posterior summaries are basically the summaries from the prior distribution. 
#' If FALSE, then Data must be provided and posterior will be computed. 
#' @param Pcutoffs target toxicity interval -- a vector containg two probability elements -- first 
#' is the boundary point between underdose and target toxicity region, second is the boundary 
#' between target toxicity and overdose region
#' @param Pcutoffs3 same as above for grade 3
#' @param Pcutoffs4 same as above for grade 4
#' @param DLT.cutoff adverse event grade cut off for DLT, must be 3 or 4.
#' If 3 then DLT >= grade 3, if 4 then DLT >= grade 4
#' @param Outfile output file
#' @param Plot Indicator if plot will be generated
#' @param nburnin number of burnin samples
#' @param nmc number of markov chain monte carlo samples after burnin
#' @param nthin thinning interval, default is 1 (no thinning)
#' @param nchain number of markov chains 
#' @param nsample number of samples from prior distribution. This is only required if 
#' Prior.Only = TRUE
#' @param int.crit Critical value for interval to change color
#' @param alpha.min the lowest point of the intercept parameter space
#' @param alpha.max the highest point of the intercept parameter space
#' @param DIC DIC to save
#' @param seed Seed for R for reproducibility
#'
#'
#' @details For single-agent data, the base model, is a standard cumulative logistic 
#' proportional odds model in log-dose. We place the Normal(prior.mean, prior.var) distribution on the 
#' regression paramter \eqn{\beta}
#' 
#' 
#' @references Meter, E. M. V., Garrett-Mayer, E., & Bandyopadhyay, D. (2010).
#' Proportional odds model for dose-finding clinical trial designs with ordinal
#' toxicity grading. Statistics in Medicine, 30(17), 2070-2080.
#' 
#' Darssan, D., Thompson, M. H., & Pettitt, A. N. (2014). 
#' Incorporating adverse event relatedness into dose-finding clinical trial designs.
#' Statistics in Medicine, 33(7), 1146-1161.
#' 
#' @return \item{Data}{The dataset}
#' \item{Agent}{The agent's name}
#' \item{prior.mean}{Prior mean used for the analysis}
#' \item{prior.var}{Prior variance used for the analysis}
#' \item{beta}{Posterior summary of \eqn{\beta}}
#' \item{grade.prob}{Predictive probabilities of doses for each adverse event grade}
#' \item{P}{Posterior predictive summaries of the dose levels}
#' \item{Pcat}{Posterior predictive probabilities of doses}
#' \item{pred.Tox}{Posterior predictive toxicity probilities of the future doses}
#' \item{Doses}{Doses for compound at which the prediction is computed}
#' \item{nburnin}{number of burnin samples}
#' \item{nmc}{number of markov chain monte carlo samples after burnin}
#' \item{seed}{Seed for R used in the computation}
#' In addition, a text output containing all of the above and two plots -- posterior DLT rates and 
#' interval probabilities of DLT's, are produced. 
#' 
#' 
#' @examples 
#' \dontrun{
#' data(crs.pk)
#' 
#' y        <- as.factor(crs.pk$DLT)
#' x        <- log(as.vector(crs.pk$Dose))
#' DoseAdm  <- x  
#' Toxicity <- y
#' Data     <- data.frame(Toxicity, DoseAdm)
#' Doses    <- sort(unique(DoseAdm))  # we want to find the toxicity probability at this dose
#' 
#' 
#' fit <- bordCRM(Data       = Data,
#'                Agent      = "Example",
#'                Doses      = Doses,
#'                DoseRef    = NULL,
#'                prior.mean = 0,
#'                prior.var  = 1,
#'                Prior.Only = FALSE,
#'                Pcutoffs   = c(0.16, 0.33),
#'                Pcutoffs3  = c(0.10, 0.25),
#'                Pcutoffs4  = c(0.05, 0.10),
#'                DLT.cutoff = 3,
#'                Outfile    = "out_ordinal",
#'                Plot       = TRUE,
#'                nburnin    = 500,
#'                nmc        = 1000,
#'                nlev       = 5,
#'                alpha.min  = -1000,
#'                alpha.max  = 1000,
#'                seed       = 1452671
#'                )
#'                
#' fit$Pcat 
#' fit$Pcat3
#' fit$Pcat4
#' fit$grade.prob
#'}
#'
#'
#'
#' @export
#' 
#' 
#' @importFrom stats is.empty.model median model.matrix model.response rnorm na.omit
#' runif sd



bordCRM <- function(Data       = NULL,
                    Agent      = NULL,
                    Doses      = NULL,
                    DoseRef    = NULL,
                    nlev       = NULL,
                    prior.mean = 0,
                    prior.var  = 1,
                    Prior.Only = FALSE,
                    Pcutoffs   = c(0.16, 0.33),
                    Pcutoffs3  = c(0.10, 0.25),
                    Pcutoffs4  = c(0.05, 0.10),
                    DLT.cutoff = 3,
                    Outfile    = "out_ordinal",
                    Plot       = TRUE,
                    nburnin    = 1000,
                    nmc        = 5000,
                    nthin      = 1,
                    nchain     = 1,
                    nsample    = 1000,
                    int.crit   = c(1, 1, 0.25),
                    alpha.min  = -1000,
                    alpha.max  = 1000,
                    DIC        = NULL,
                    seed       = 1452671
)
{
  
  
  version <- "16-July-2019"
  
  if(is.null(Data) && (Prior.Only == FALSE))
    stop("Either provide data OR set Prior.Only = TRUE")
  
  set.seed(seed)
  
  if(Prior.Only){DIC=FALSE}
  
  
  
  # Reading Data
  
  DoseAdm <- Data$DoseAdm
  Tox     <- Data$Toxicity
  Npat    <- nrow(Data)
  if(is.numeric(DoseRef))
  {
    DoseAdm/DoseRef
  }
  
  if(length(DoseAdm) != length(Tox))
    stop("All the number of doses is not equal to the number of toxicities")
  # DosesAdm[DosesAdm ==0] = 0.00001
  
  Ndoses <- length(Doses)
  
  Ncat <- length(Pcutoffs) + 1
  
  Pcutoffs1 <- c(0, Pcutoffs, 1)
  Pcutoffs31 <- c(0, Pcutoffs3, 1)
  Pcutoffs41 <- c(0, Pcutoffs4, 1)
  
  intLabels  <- paste(Pcutoffs1[1:Ncat], "-", Pcutoffs1[2:(Ncat + 1)], sep = "")
  intLabels3 <- paste(Pcutoffs31[1:Ncat], "-", Pcutoffs31[2:(Ncat + 1)], sep = "")
  intLabels4 <- paste(Pcutoffs41[1:Ncat], "-", Pcutoffs41[2:(Ncat + 1)], sep = "")
  
  Doses.Label = paste(Agent, "=", exp(Doses), sep = '')
  
  
  if(!Prior.Only){
    
    fit <- bord(formula    = as.factor(Tox) ~ DoseAdm, 
                data       = NULL, 
                prior.mean = prior.mean,
                prior.var  = prior.var,
                nburnin    = nburnin,
                nmc        = nmc,
                newx       = Doses,
                seed       = seed,
                Pcutoffs   = Pcutoffs,
                Pcutoffs3  = Pcutoffs3,
                Pcutoffs4  = Pcutoffs4,
                DLT.cutoff = DLT.cutoff,
                nlev       = nlev
    )
    
  } else {
    
    fit <- bord.prior(formula    = as.factor(Tox) ~ DoseAdm, 
                      data       = NULL, 
                      prior.mean = prior.mean,
                      prior.var  = prior.var,
                      nsample    = nsample,
                      newx       = Doses,
                      seed       = seed,
                      Pcutoffs   = Pcutoffs,
                      Pcutoffs3  = Pcutoffs3,
                      Pcutoffs4  = Pcutoffs4,
                      DLT.cutoff = DLT.cutoff,
                      nlev       = nlev
    )
  }
  
  #Processing output
  
  if(length(unique(Doses)) == 1)
  {
    P <- fit$pTox.post[c(4, 1, 6)]
  } else {
    P <- fit$pTox.post[c(4, 1, 6), ] 
  }
  
  Pcat           <- fit$pCat.post 
  Pcat3          <- fit$pCat3.post
  Pcat4          <- fit$pCat4.post
  pred.Tox       <- fit$predTox.post
  pred.dose      <- exp(fit$pred.dose)
  beta           <- fit$beta.summary
  grade.prob.hat <- fit$grade.prob.post
  grade.prob     <- fit$grade.prob
  
  Pcat <- BUGS.Table2Array.bord(Pcat, Labels = list(Doses.Label, intLabels), 
                           Dim = c(Ndoses, Ncat))[,, 1]
  Pcat3 <- BUGS.Table2Array.bord(Pcat3, Labels = list(Doses.Label, intLabels3), 
                                Dim = c(Ndoses, Ncat))[,, 1]
  Pcat4 <- BUGS.Table2Array.bord(Pcat4, Labels = list(Doses.Label, intLabels4), 
                                Dim = c(Ndoses, Ncat))[,, 1]
  
  # Collecting outputs
  
  outlist = list(Data       = Data, 
                 Agent      = Agent,
                 prior.mean = prior.mean,
                 prior.var  = prior.var,
                 beta       = beta,
                 grade.prob = grade.prob.hat,
                 P          = P,
                 Pcat       = Pcat,
                 Pcat3      = Pcat3,
                 Pcat4      = Pcat4,
                 pred.Tox   = pred.Tox,
                 pred.dose  = pred.dose,
                 nburnin    = nburnin,
                 nmc        = nmc, 
                 seed       = seed
  )
  
  if (nlevels(as.factor(Doses)) < 2L)
  {
    Plot = FALSE
    warning("Some plots will be removed as the prediction doses have only 1 level")
  }
  
  
  if(Plot == TRUE){
    
    P1 <- t(P)
    colnames(P1) <- c("mean", "2.5%", "97.5%")
    
    Plot.P(Est = P1, RowNames = exp(sort(unique(Doses))))  # this will plot locally
    Plot.Pcat(Pcat = Pcat, crit = int.crit, RowNames = pred.dose)  # this will plot locally
    
    grDevices::pdf(file = paste(Outfile,"_P",'.pdf',sep=''), onefile = F)
    Plot.P(Est = P1, RowNames = exp(sort(unique(Doses))))
    grDevices::dev.off()
    
    grDevices::pdf(file = paste(Outfile,"_Pcat",'.pdf',sep=''), onefile = F)
    Plot.Pcat(Pcat = Pcat, crit = int.crit, RowNames = pred.dose) 
    grDevices::dev.off()
    
  }
  
  # Grade probabilities plot by predictive doses
  box.data <- NULL
  for(idose in 1:length(Doses))
  {
    for(ilev in 1:nlev)
    {
      box.data <- rbind(box.data, cbind(exp(Doses[idose]), ilev, fit$grade.prob[idose, ilev, ]))
    }
  }
  colnames(box.data) <- c("Dose", "Grade", "Probability")
  box.data <- na.omit(data.frame(box.data))
  box.data$Dose <- as.factor(box.data$Dose)
  box.data$Grade <- as.factor(box.data$Grade)
  
  if(0 %in% Tox)
  {
    boxplot.label <- as.character(1:nlev - 1)
  } else {
    boxplot.label <- as.character(1:nlev)
  }
  
  Probability <- NULL  # to avoid NOTE of cmd_check()
  Grade <- NULL        # to avoid NOTE of cmd_check()
  
  plot <- ggplot2::ggplot(box.data, 
                          ggplot2::aes(y = Probability, fill = Grade)) + 
    ggplot2::geom_boxplot() + ggplot2::facet_wrap(~ Dose, ncol = 3) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +  # remove x axis label
    ggplot2::scale_fill_discrete(name = "Grade", labels = boxplot.label)
  print(plot)
  
  grDevices::pdf(file = paste(Outfile, "_Pgrade", '.pdf', sep=''), onefile = F)
  print(plot)
  grDevices::dev.off()
  
  
  #Creating Output file with necessary Outputs
  
  sink(paste(Outfile, ".txt", sep = ""), append = F)
  cat("\n --------------------------------------------------------------------------------\n")
  cat(" Dose Escalation: Bayesian Inference")
  cat("\n", date())
  cat("\n Working directory: ", getwd())
  cat("\n --------------------------------------------------------------------------------\n")
  
  cat("\n\n")
  cat(" Data")
  cat("\n\n")
  print(outlist$Data)
  cat("\n")
  
  cat("\n --------------------------------------------------------------------------------\n")
  cat(" Posterior Results")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("\n")
  cat(" Parameters")
  cat("\n")
  print(outlist$beta)
  cat("\n")
  cat("\n")
  cat("\n")
  
  cat(" Grade Probabilities")
  cat("\n")
  out.matrix <- cbind(pred.dose, outlist$grade.prob)
  colnames(out.matrix)[1] <- c("Doses")
  print(out.matrix)
  cat("\n")
  cat("\n")
  cat("\n")
  
  cat(" DLT rates")
  cat("\n")
  if(length(unique(Doses)) == 1)
  {
    out.matrix <- matrix(c(exp(sort(unique(Doses))), outlist$P), nrow = 4, ncol = 1)
  } else {
    out.matrix <- rbind(exp(sort(unique(Doses))), outlist$P)
  }
  rownames(out.matrix) <- c("DoseAdm", "Minimum", "Mean", "Maximum")
  print(out.matrix)
  cat("\n")
  
  cat("\n")
  cat(" Interval Probabilities")
  cat("\n")
  if(length(unique(Doses)) == 1)
  {
    out.matrix <- matrix(c(exp(sort(unique(Doses))), outlist$Pcat), nrow = 1, ncol = 4)
  } else {
    out.matrix <- cbind(pred.dose, outlist$Pcat)
  }
  colnames(out.matrix) <- c("Doses", "Underdose Probability", "Target Toxocity Probability", 
                            "Overdose Probability")
  print(out.matrix)
  cat("\n")
  
  cat("\n")
  cat(" Predictive Toxicities (Posterior Means)")
  cat("\n")
  out.matrix <- rbind(pred.dose, outlist$pred.Tox)
  rownames(out.matrix) <- c("Doses", "Toxocity probability")
  print(out.matrix)
  cat("\n")
  
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat(" Model Specification")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("\n")
  cat(" Prior distributions")
  cat("\n")
  print(c(outlist$prior.mean, outlist$prior.var))
  cat("\n")
  
  cat("\n")
  cat(" MCMC parameters")
  cat("\n")
  print(c(outlist$nburnin, outlist$nmc))
  cat("\n")
  sink()
  
  return(outlist)
  
}
