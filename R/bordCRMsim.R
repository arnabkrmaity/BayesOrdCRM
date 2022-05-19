#' @title bordCRMsim Function
#'
#' @description This function performs Phase I Oncology trial simulation for single agent using Bayesian 
#' Ordinal Regression Model (bordCRM) and Escalation with Overdose Criteria (EWOC)
#'
#' @param Design   Design specification: contains required parameters for \code{\link{SingleBLRM}}. This incudes provisional
#' dose levels, reference dose, and prior specifications. See \code{\link{SingleBLRM}} for more details.
#' @param nTrials Number of simulated trials.
#' @param TrueProbs a vector of true toxicity  probabilities for different dose levels.
#' @param StartDose Start dose for the escalation trial
#' @param CohortRnd Cohort size (Random). A list with the first element a vector of the available
#'                  cohort size(s), and the second element a vector of the same length with the 
#'                  respective selection probability of the cohort size. For example, when set to  
#'                  list(c(3, 2), c(0.7, 0.3)), then the cohortsize is chosen as 3 or 2 at random 
#'                  with probability 0.7 or 0.3 respectively. If a fixed cohort size is desired use 
#'                  e.g. CohortRnd = list(c(3), c(1))
#' @param Pcutoffs cutoffs to define the intervals for the true DLT probability. The default is
#' (0.16, 0.33) defining three categories: underdosing ({[}0--0.16{]}), target ((0.16--0.33)), and 
#' overdosing ({[}0.33--1.00{]}).
#' @param Nsim four component vector of MCMC parameters n.burnin, n.iterations, n.chains, n.thin.
#' @param Decision Function which decides after each simulated cohort how to proceed
#' @param Outfile A R object saves all simulations output objects
#' @param Out.summary A txt file with simulation summary
#' @param RndSeed Seed used for the random number generator. Required parameter
#'
#' @details See vignettes. browseVignettes("OncoPh1BLRM")
#' 
#' @references 
#' #'Neuenschwander B, Branson M, Gsponer T. Critical aspects of the Bayesian approach to phase I
#'cancer trials. Stat Med. 2008;27(13):2420-39.
#' 
#' Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. A Bayesian Industry
#' Approach to Phase I Combination Trials in Oncology. Statistical Methods for Drug Combination
#' Studies. Taylor & Francis, 2015.
#' 
#' @return Input arguments are 
#' \item{Design}{The input data and prior information}
#' \item{nTrials}{Number of simulated trials}
#' \item{TrueProbs}{True toxicity vector}
#' \item{Pcutoffs}{Target toxicity interval}
#' \item{CohortRnd}{A list of size two with first quantity is the vector of cohortsize and second quantity is 
#'                  the vector of probabilities assigned to each cohortsize. For example, when set to 
#'                  list(c(3, 2), c(0.7, 0.3)), then the cohortsize is chosen as 3 or 2 at random 
#'                  with probability 0.7 or 0.3 respectively.}
#' \item{Decision}{}
#' \item{Outfile}{Name of the output file}
#' \item{Out.summary}{Name of the summary output file}
#' \item{RndSeed}{The seed which was fixed to reproduce the result}
#' 
#' Output arguments are 
#' \item{DLTMATRX}{the occurrences of DLT's per dose per trial}
#' \item{MTDSIM}{a binary matrix which is only TRUE when MTD has been observed in the trial}
#' \item{NNMATRX}{the cohort sizes per dose per trial}
#' \item{res.sim$simul}{descrioption of each trial -- Trial no. (Trial), Cohort no. (Cohort), 
#' Number of toxicity (Ntox), Number of patients (Npat), Weight, Current dose (DoseAdm), 
#' underdose probability (under), probability of target toxicity (target), 
#' overdose probability (over), probability of the current dose (true.prob), 
#' suggestion for the next dose, if any, (nextDoseAdm), underdose probability of the next 
#' dose (nextunder), target toxicity probability of the next dose (nexttarget), 
#' probability of the next dose (nextprob), a binary variable which is only 1 (TRUE) 
#' if the current dose is declared as the MTD, if the first dose is so toxic that 
#' \eqn{\Pr(\pi_d > 0.33) \ge EWOC} then the trial is stooped and too.toxic variable is set to 1, 
#' the Rhat.max information, crash = 1 if the code is crashed, the seed number for R (Rseed) 
#' and JAGS (WBseed)}
#' \item{res.sim$simsummary}{description of those cohorts within those trials 
#' in which the corresponding dose is the MTD}
#' \item{sim.sum}{simulation summary -- Proportion of patients in target dose, 
#' Proportion of patients in over dose, Proportion of patients in under dose,
#' Proportion of trials with MTD within target dose region, 
#' Proportion of trails with MTD within over dose region, 
#' Proportion of trails with MTD within under dose region, 
#' Proportion of trials stopped due to too toxicity at the start dose, 
#' Average sample size}
#' 
#' In addition, an R object output and a text output are generated containing all of the above. 
#' 
#' @seealso \code{\link{SingleBLRM}}, \code{\link{SingleBLRMSimDecisionRule}}
#' 
#' @examples 
#' \dontrun{
#' 
#' data(crs.pk)
#' 
#' Doses <- sort(unique(crs.pk$Dose))  # we want to find the toxicity probability at this dose
#' 
#' 
#' # Design Parameters
#' Design <- list(Agent1     = "Agent 1",
#'                Doses      = Doses,
#'                DoseRef    = NULL,
#'                prior.mean = 0,
#'                prior.var  = 1
#' )
#' 
#' # Simulation Scenario 1
#' TrueProbs.prior <- c(0.05, 0.11, 0.18, 0.25, 0.28, 0.31, 0.36)
#' 
#' #Running Simulation
#' sim.args <- list(Design      = Design,
#'                  nTrials     = 10,
#'                  TrueProbs   = TrueProbs.prior,
#'                  StartDose   = Doses[1],
#'                  Pcutoffs    = c(0.16, 0.33),
#'                  CohortRnd   = list(c(3), c(1)),
#'                  Decision    = bordCRMSimDecisionRule(MaxN = 60, mnMTD = 6, 
#'                                                       mnTOT = 15, mnTAR = 0.5, 
#'                                                       EWOC = 0.25, 
#'                                                       Escalation  = "maximal", 
#'                                                       # maxiumum dose which satisfies EWOC criteria
#'                                                       doseMaxEsc  = c(2, 3, 4, 3, 4, 2)
#'                                                       ),
#'                  Outfile     = "Output_Simulation",
#'                  Out.summary = "Output_Summary_Simulation.txt",
#'                  RndSeed     = 123
#' )
#' 
#' # Setting up back-end for parallel computing (DO NOT CHANGE)
#' library(doParallel)
#' require(doRNG)
#' nCores <- detectCores()  # number of cores
#' cl <- makeCluster(nCores) 
#' registerDoParallel(cl)
#' registerDoRNG(seed = 123, once = FALSE) #SEED1 : To initiate enviorment. 
#' # THIS IS THE IMPORTANT STEP FOR PARALLELISM AND REPRODUCIBILITY
#' getDoParWorkers()
#' 
#' simulation.bord <- do.call(bordCRMSim, sim.args)
#' 
#' # Stop all clusters (DO NOT CHANGE)
#' OncoPh1BLRM:::unregister()
#' stopImplicitCluster() 
#' 
#' simulation.bord$sim.sum
#' }
#' 
#' 
#' @export


bordCRMSim <- function (Design      = NULL,
                        nTrials     = 1000,
                        TrueProbs   = NULL,
                        StartDose   = NULL,
                        CohortRnd   = list(c(  3,   4,    5,    6),
                                         c(0.8, 0.1, 0.05, 0.05)),
                        Pcutoffs    = c(0.16, 0.33),
                        Nsim        = c(1000, 2000, 2, 1),
                        Decision    = bordCRMSimDecisionRule(MaxN  = 60, mnMTD = 6, 
                                                             mnTOT = 15, mnTAR = 0.5, 
                                                             EWOC  = 0.25, 
                                                             Escalation = "maximal", 
                                                             doseMaxEsc = 2),
                        Outfile     = "Sims",
                        Out.summary = "sims_summary.txt",
                        RndSeed     = 12
)
{

  set.seed(RndSeed)
  SeedParam <- RndSeed

  #SeedParam <- sample.int(10, nTrials)

  # suppressPackageStartupMessages(require(foreach, warn.conflicts=FALSE, quietly = TRUE))
  # suppressPackageStartupMessages(require(doRNG, warn.conflicts=FALSE, quietly = TRUE))

  if(missing(SeedParam)) {
    stop("ERROR: Random seed must be set for this operation!\nPlease specify SeedParam.")
  }
  

  forArgs <- list(i = 1:nTrials, .options.RNG=SeedParam)

  if(foreach::getDoParRegistered()) {
    ## for the MPI backend set the chunksize to be at most 10
    if(getDoParName() == "doMPI") {
      forArgs$.options.mpi <- list(chunkSize=min(10, ceiling(nTrials/getDoParWorkers())))
      cat("Setting doMPI backend chunkSize to", forArgs$.options.mpi$chunkSize, "\n")
    }
  }


  #set.seed(RndSeed)

  Results <- do.call(foreach::foreach, forArgs) %dorng% {
    cat("Running trial", i, "\n")

    Ntot <- 0
    is.MTD <- FALSE
    stop.trial <- FALSE
    too.toxic <- FALSE
    SimulationTrial <- matrix(nrow = 50, ncol = 21)
    colnames(SimulationTrial) = c("Trial", "Cohort", "Ntox",
                                  "Npat", "Weight", "DoseAdm","under",
                                  "target", "over", "true.prob", "nextDoseAdm",
                                  "nextunder", "nexttarget", "nextover", "nextprob",
                                  "MTD", "too.toxic", "Rhat.max", "crash", "Rseed",
                                  "WBseed")


    Agent   <- Design$Agent
    Doses   <- Design$Doses
    DoseRef <- Design$DoseRef

    prior.mean <- Design$prior.mean
    prior.var  <- Design$prior.var

    DosesAdm.obj <- NULL
    Toxicity.obj <- NULL

    currentDoseAdm <- StartDose

    cohort  <- 1
    crashed <- FALSE
    
    
    while (!stop.trial & !crashed) {
      if (cohort > nrow(SimulationTrial)) {
        SimulationTrialL <- matrix(nrow = nrow(SimulationTrial) + 15, ncol = 21)
        colnames(SimulationTrialL) <- colnames(SimulationTrial)
        SimulationTrialL[1:nrow(SimulationTrial), ] <- SimulationTrial
        SimulationTrial <- SimulationTrialL
      }
      SimulationTrial[cohort, "Trial"]  <- i
      SimulationTrial[cohort, "Cohort"] <- cohort
      SimulationTrial[cohort, "Rseed"]  <- i

      if (length(CohortRnd[[1]]) == 1) {
        CohortSize <- CohortRnd[[1]]
      }
      else {
        CohortSize <- sample(CohortRnd[[1]], 1, prob = CohortRnd[[2]])
      }
      
      
      tox.prob <- TrueProbs[which(Doses == currentDoseAdm)]
      
      if(tox.prob <= 0.1) {
        tox.vec  <- c(0.8, 0.1, 0.05, 0.04, 0.01)
      } else if ((tox.prob > 0.1) & (tox.prob <= 0.2)) {
        tox.vec  <- c(0.6, 0.2, 0.1, 0.08, 0.02)
      } else if ((tox.prob > 0.2) & (tox.prob <= 0.3)) {
        tox.vec  <- c(0.4, 0.3, 0.15, 0.1, 0.05)
      } else if ((tox.prob > 0.3) & (tox.prob <= 0.35)) {
        tox.vec  <- c(0.2, 0.2, 0.3, 0.2, 0.1)
      } else {
        tox.vec  <- c(0.1, 0.1, 0.3, 0.3, 0.2)
      }
      
      
      Ntox.vec <- apply(X = as.matrix(c(1:CohortSize)), MARGIN = 1, FUN = function (x) {
        sample(c(1:5), size = 1, prob = tox.vec)})
      
      
      DosesAdm.obj <- c(DosesAdm.obj, rep(currentDoseAdm, CohortSize))
      Toxicity.obj <- c(Toxicity.obj, Ntox.vec)
      
      
      # Dose Limiting Toxicity assuming toxicity is >= 4
      temp.DLT  <- ifelse(Toxicity.obj >= 4, 1, 0) 
      temp.data <- cbind(NTOX = temp.DLT, DOSEADM = DosesAdm.obj, NPAT = 1)
      

      SimulationTrial[cohort, "Npat"] <- CohortSize
      Ntot <- Ntot + CohortSize
      SimulationTrial[cohort, "DoseAdm"] <- currentDoseAdm
      SimulationTrial[cohort, "Ntox"] <- sum(temp.DLT)  #paste0(Ntox.vec, collapse = "-")
      SimulationTrial[cohort, "true.prob"] <- TrueProbs[which(Doses == currentDoseAdm)]
      SimulationTrial[cohort, "Weight"] <- 1

      
      current.data <- data.frame(Toxicity = as.factor(Toxicity.obj), DoseAdm = log(DosesAdm.obj))

      n.att <- 1
      did.run <- FALSE


      sink(paste(Outfile, "_log.txt", sep=""))

      while((n.att <= 10) & !did.run) {
        simSeeds <- sample.int(.Machine$integer.max, 2)
        currentRun <- try(bordCRM(Data       = current.data,
                                  Agent      = Agent,
                                  Doses      = log(Doses),
                                  DoseRef    = DoseRef,
                                  prior.mean = prior.mean,
                                  prior.var  = prior.var,
                                  Prior.Only = FALSE,
                                  Pcutoffs   = Pcutoffs,
                                  Outfile    = "out_ordinal",
                                  Plot       = FALSE,
                                  nburnin    = Nsim[1],
                                  nmc        = Nsim[2],
                                  nlev       = 5,
                                  alpha.min  = -1000,
                                  alpha.max  = 1000,
                                  seed       = simSeeds[1]
                                  ))
        did.run <- is.null(attr(currentRun, "class"))
        n.att <- n.att + 1
        SimulationTrial[cohort, "WBseed"] <- simSeeds[2]
      }
      sink()

      #print(currentRun$Pcat)

      if (!did.run) {
        crashed <- TRUE
        SimulationTrial[, "crash"] <- TRUE
      } 
      else {
        SimulationTrial[, "crash"] <- FALSE
        
        SimulationTrial[cohort, "under"]  <- currentRun$Pcat[, 1][which(currentDoseAdm == Doses)]
        SimulationTrial[cohort, "target"] <- currentRun$Pcat[, 2][which(currentDoseAdm == Doses)]
        SimulationTrial[cohort, "over"]   <- currentRun$Pcat[, 3][which(currentDoseAdm == Doses)]

        #print(SimulationTrial)
        
        # prepare data for decision rule
        currentRun$Data$data <- temp.data
        
        
        

        decisionVars <- names(formals(Decision))
        decision <- do.call(Decision, mget(decisionVars))

        too.toxic <- decision$too.toxic
        stop.trial <- decision$stop.trial
        is.MTD <- decision$is.MTD
        nextDoseAdm <- decision$DoseAdm

        if (!too.toxic) {

          SimulationTrial[cohort, "nextunder"] <- currentRun$Pcat[, 1][which(nextDoseAdm == Doses)]
          SimulationTrial[cohort, "nexttarget"] <- currentRun$Pcat[, 2][which(nextDoseAdm == Doses)]
          SimulationTrial[cohort, "nextover"] <- currentRun$Pcat[, 3][which(nextDoseAdm == Doses)]
          SimulationTrial[cohort, "nextprob"] <- TrueProbs[which(Doses == nextDoseAdm)]


          SimulationTrial[cohort, "MTD"] <- is.MTD
          SimulationTrial[cohort, "too.toxic"] <- too.toxic
          SimulationTrial[cohort, "nextDoseAdm"] <- nextDoseAdm


        }

        else {
          SimulationTrial[cohort, "nextunder"] <- NA
          SimulationTrial[cohort, "nexttarget"] <- NA
          SimulationTrial[cohort, "nextover"] <- NA
          SimulationTrial[cohort, "MTD"] <- is.MTD
          SimulationTrial[cohort, "too.toxic"] <- too.toxic
          SimulationTrial[cohort, "nextprob"] <- NA
          SimulationTrial[cohort, "nextDoseAdm"] <- nextDoseAdm
            }
        currentDoseAdm <- nextDoseAdm
        cat(paste("Cohort: ", cohort, "\n", sep = ""))
        cohort <- cohort + 1
        #print(SimulationTrial)
      }
    }

    mtd.done <- subset(SimulationTrial, SimulationTrial[, "MTD"] == 1)
    simsummaryTrial <- mtd.done
    #print(SimulationTrial[1:(cohort - 1), ])
    ncohort <- matrix(NA, nrow = NROW(Doses), ncol = 1)
    dltcohort <- matrix(NA, nrow = NROW(Doses), ncol = 1)
    for (l in 1:NROW(Doses)) {
        ncohort1 <- subset(SimulationTrial, SimulationTrial[, "DoseAdm"] == Doses[l])
        ncohort[l,1] <- sum(ncohort1[, "Npat"])
        dltcohort[l,1] <- sum(ncohort1[, "Ntox"])
                 }

    NMATRXTrial <- ncohort/sum(ncohort)
    NNMATRXTrial <- ncohort
    DLTMATRXTrial <- dltcohort
    cat(paste("Trial ", i, ": finalized.\n", sep = ""))

    res <- list(Simulation = SimulationTrial, simsummary = simsummaryTrial,
                NMATRX = NMATRXTrial, NNMATRX = NNMATRXTrial, DLTMATRX = DLTMATRXTrial)
    res
  }


  Doses = Design$Doses
  Agent =Design$Agent

  extract <- function(struct, elem) lapply(struct, "[[", elem)
  Simulation <- extract(Results, "Simulation")
  simsummary <- extract(Results, "simsummary")
  NMATRX <- extract(Results, "NMATRX")
  NNMATRX <- extract(Results, "NNMATRX")
  DLTMATRX <- extract(Results, "DLTMATRX")
  RNGSEEDS <- attr(Results, "rng")



  MTDsim <- list()
  simsummarys <- do.call("rbind", simsummary)
  if (NROW(simsummarys) > 0) {
    for (q in 1:nrow(simsummarys)) {
      MTDmatrix <- matrix(NA, nrow = NROW(Doses), ncol = 1)
      for (a in 1:NROW(Doses)) {
          MTDmatrix[a, 1] <- (simsummarys[q, "DoseAdm"] == Doses[a])
                               }
      MTDsim[[q]] <- MTDmatrix
    }
  }

  if (NROW(simsummarys) == 0) {
    MTDmatrix <- matrix(NA, nrow = NROW(Doses), ncol = 1)

    MTDsim[[1]] <- MTDmatrix

  }

  res.sim=list(simul = Simulation,
               simsummary = simsummarys,
               NMATRX = NMATRX,
               NNMATRX = NNMATRX,
               DLTMATRX = DLTMATRX,
               MTDsim = MTDsim,
               RNGSEEDS = RNGSEEDS)


  sim.sum.Vars <- names(formals(bordCRMSimSummary))

  sim.sum <- do.call(bordCRMSimSummary, mget(sim.sum.Vars))

  save(list = names(res.sim), file = paste(Outfile, ".Rdata", sep = ""), envir = list2env(res.sim))

  res.sim.f <- list(res.sim=res.sim, sim.sum= sim.sum)

  return(res.sim.f)

}
