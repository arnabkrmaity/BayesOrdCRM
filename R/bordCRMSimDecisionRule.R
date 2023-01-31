#' @title bordCRMSimDecisionRule Function
#'
#' @description This function performs escalation and decision rule trial simulation with 
#' single agent
#' 
#' @param MaxN Maximum number of patients in each trial. The default is 60. See Details
#' @param mnMTD At least these many patients have been treated at MTD. The default is 6. See Details
#' @param mnTOT A minimum number of these many patients should have already treated in the trial. 
#' The default is 15. See Details
#' @param mnTAR Probability thresold of target toxicity. The default is 0.5. See details
#' @param EWOC Escalation with overdose control criterion. The default is 0.25. See details
#' @param Escalation Escalation The escalation method to be used. Possible methods -- "maximal":
#' recomended, maxiumum dose which satisfies EWOC criteria; "target": dose which has 
#' highest target toxicity probability; "dose_max": maximum dose, doses will not be skipped  
#' @param doseMaxEsc Either a scaler or a vector of length of number of provisional doses. If a scaler
#' \eqn{k} then it used assumed that for all doses escaltion upto \eqn{k}-fold is allowed. If a vector
#' then for each dose the corresponding fold is allowed to escalate upto that dose. The first element is set to 2.
#' 
#' @details At any trial dose escalation is achieved satisfying the escalation with overdose control 
#' (EWOC) principle. This is implemented by setting the target toxicity interval to Pcutoffs, with
#' an overdose control criterion \eqn{P(\pi_d > Pcutoffs[2]) < EWOC}. 
#' 
#' Dose escalations continue until declaration of the maximum tolerated dose (MTD) \eqn{d}. 
#' This dose must meet the following conditions:
#' 
#' 1. At least mnMTD patients have been treated at dose \eqn{d}.
#' 
#' 2. This dose satisfies one of the following conditions:
#' 
#' i. The probability of targeted toxicity at \eqn{d} exceeds mnTAR.
#' 
#' ii. A minimum number of mnTOT patients should have already been treated in the trial. 
#' 
#' 
#' @seealso \code{\link{bordCRMsim}}
#' 
## #' @export


bordCRMSimDecisionRule <- function (MaxN,
                                    mnMTD,
                                    mnTOT,
                                    mnTAR,
                                    EWOC,
                                    Escalation = c("dose_max", "maximal", "target"),
                                    doseMaxEsc = 2
)
{


  Escalation <- match.arg(Escalation)
  if (any(doseMaxEsc <= 1))
    warning("Some dose escalation factors are smaller than or equal to 1.")
 
  
  function(currentRun, Ntot, currentDoseAdm, Doses){

    CDoseAdm <- currentDoseAdm

    datsummary <- currentRun$Data$data
    names(datsummary) <- toupper(names(datsummary))

    Pall <- currentRun$Pcat
    xrng <- max(which(Doses >= CDoseAdm))
    totNPAT <- sum(datsummary[, "NPAT"][datsummary[, "DOSEADM"] ==  CDoseAdm])
    totNTOX <- sum(datsummary[, "NTOX"][datsummary[, "DOSEADM"] ==  CDoseAdm])


    odprob <- Pall[, 3]
    tgprob <- Pall[, 2]
    dose.eligible <- odprob < EWOC #0.25
    lDeltaDose <- log(if(length(doseMaxEsc) == 1) rep(doseMaxEsc, length(Doses)) else doseMaxEsc)
    ldm <- log(Doses/CDoseAdm) <= lDeltaDose
    possible.doses <- (ldm  == 1)
    ndset <- (dose.eligible + possible.doses) == 2

    DM.H <- Doses >= CDoseAdm
    ndsetup <- (DM.H + dose.eligible + possible.doses) == 3


    if(sum(ndset)== 0) {
      too.toxic <- TRUE
      is.MTD <- FALSE
      DoseAdm <- NA
    }
    else{
      too.toxic <- FALSE
        if (Escalation == "dose_max"){
          dose_cur <- which(CDoseAdm == Doses)
          while (all(!ndset[dose_cur, ])) {
            dose_cur <- dose_cur - 1
          }

          dose_next <- dose_cur
          while (dose_next != length(Doses) && ndset[dose_next + 1]) {
            dose_next <- dose_next + 1
          }

          dose.pair <- dose_next
        }
      
      if (Escalation == "maximal"){
          odprob1 <- odprob
          odprob1[!ndset] <- NA
          dose.pair <- max(which(odprob1 == max(odprob1, na.rm = T), arr.ind = T))
      }
      
      if (Escalation == "target"){
        if (sum(ndsetup) > 0){
          tgprob1 <- tgprob
          tgprob1[!ndsetup] <- NA
          dose.pair <- max(which(tgprob1 == max(tgprob1, na.rm = T), arr.ind = T))
        }

      if(sum(ndsetup) == 0 & sum(ndset)!= 0){
        tgprob[ndset == FALSE] <- NA
        dose.pair <- max(which(tgprob == max(tgprob, na.rm = T), arr.ind = T))
      }
      
      }
      
      nextDose <- dose.pair
      next.equal <- (Doses[nextDose] == CDoseAdm)
      
      boundD <- CDoseAdm == max(Doses)
      mnMTD.check <- (totNPAT >= mnMTD)

        if (!boundD) {
          current.eligible <- odprob[which(Doses ==  CDoseAdm)] < EWOC
          current.highest <- min(odprob[(which(Doses == CDoseAdm) + 1):xrng]) >= EWOC
          #current.highest <- all(current.highestd)
        }

        if (boundD) {
          current.eligible <- odprob[which(Doses == CDoseAdm)] < EWOC
          current.highest <- odprob[which(Doses == CDoseAdm)] > EWOC
        }
        MTD.rule <- (Ntot >= mnTOT || tgprob[which(Doses == CDoseAdm)] >= mnTAR)


      if (mnMTD.check & current.eligible & MTD.rule & (current.highest | next.equal)) {
        is.MTD <- TRUE
        DoseAdm <- CDoseAdm
      }
      else {
        is.MTD <- FALSE
        DoseAdm <- Doses[nextDose]
      }
    }

    stop.trial <- too.toxic | is.MTD | Ntot >= MaxN


    return(list(too.toxic = too.toxic, stop.trial = stop.trial,
                is.MTD = is.MTD, DoseAdm = DoseAdm))
  }
}
