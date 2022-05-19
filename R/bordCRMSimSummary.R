
#------------------------------------------------------------------------------------------------------------------------------------
#summerize the objects and getting metrics for protocol
#------------------------------------------------------------------------------------------------------------------------------------


bordCRMSimSummary <- function(nTrials, Agent, Doses, Pcutoffs, TrueProbs, res.sim, Out.summary){

  Sims <- res.sim


under.dose.region <- TrueProbs < Pcutoffs[1]

target.dose.regions <- TrueProbs >= Pcutoffs[1] & TrueProbs < Pcutoffs[2]

over.dose.region <- TrueProbs >= Pcutoffs[2]


#Distribution of DLT in ntrials

MTD.dose <- Reduce("+", Sims$MTDsim)
aveMTD.dose <- round(MTD.dose/nTrials, 4)
row.names(aveMTD.dose) <- paste(Agent, "=", Doses, sep='')


#number of subjects in each dose in ntrials

totsubj.dose <- Reduce("+", Sims$NMATRX)
totsubj1.dose <- Reduce("+", Sims$NNMATRX)
avesubj.dose <- round(totsubj1.dose/nTrials,4)
row.names(avesubj.dose) <- paste(Agent, "=", Doses, sep='')

#number of DLT in each dose in ntrials

totDLT.dose <- Reduce("+", Sims$DLTMATRX)
aveDLT.dose <- round(totDLT.dose/nTrials, 4)
row.names(aveDLT.dose) <- paste(Agent, "=", Doses, sep='')


#Metric I: Average proportion of patients in target dose (>=Pcutoffs[1]- Pcutoffs[2])

totsubj.doseI <- totsubj.dose

totsubj.doseI[target.dose.regions== FALSE] <- 0

agvsubj.doseI <-sum(totsubj.doseI)/nTrials #sum(totsubj.dose)


#Metric II: Average proportion of patients in over dose (>= Pcutoffs[2])

totsubj.doseII <- totsubj.dose

totsubj.doseII[over.dose.region== FALSE] <- 0

agvsubj.doseII <-sum(totsubj.doseII)/nTrials #sum(totsubj.dose)


#Metric III: Average proportion of patients in under dose (< Pcutoffs[1])

totsubj.doseIII <- totsubj.dose

totsubj.doseIII[under.dose.region== FALSE] <- 0

agvsubj.doseIII <-sum(totsubj.doseIII)/nTrials #sum(totsubj.dose)

#-----------------------------------------------------------------------------------------------------------------
# How MTDs are distributed in each simulation
#----------------------------------------------------------------------------------------------------------------

MTDdist <- Reduce("+", Sims$MTDsim)

#Metric IV: Proportion of trials with MTD within target dose region (>=Pcutoffs[1]- Pcutoffs[2])
MTDdistIV <- MTDdist

MTDdistIV[target.dose.regions== FALSE] <- 0

avg.MTDIV <- sum(MTDdistIV)/nTrials


#Metric V : Proportion of trails with MTD within over dose region ( >= Pcutoffs[2])

MTDdistV <- MTDdist

MTDdistV[over.dose.region== FALSE] <- 0

avg.MTDV <- sum(MTDdistV)/nTrials


#Metric VI: Proportion of trails with MTD within under dose region (< Pcutoffs[1])
MTDdistVI <- MTDdist

MTDdistVI[under.dose.region== FALSE] <- 0

avg.MTDVI <- sum(MTDdistVI)/nTrials

#Stopped

stopd <- max(0, 1 - (avg.MTDIV + avg.MTDV + avg.MTDVI ))

#average sample size
avsamp1<- Reduce("+", Sims$NNMATRX)

avsamp <- sum(avsamp1)/nTrials

sink(Out.summary)
cat(paste("Metric I : Average proportion of patients in target dose (>=", Pcutoffs[1]*100,"%", "-", Pcutoffs[2]*100, "%)=",round(agvsubj.doseI,3),sep=""))
cat("\n")
cat(paste("Metric II : Average proportion of patients in over dose (>=",  Pcutoffs[2]*100, "%)=",round(agvsubj.doseII,3),sep=""))
cat("\n")
cat(paste("Metric III : Average proportion of patients in under dose (<", Pcutoffs[1]*100, "%)=",round(agvsubj.doseIII,3),sep=""))
cat("\n")
cat(paste("Metric IV : Proportion of trials with MTD within target dose region (>= ", Pcutoffs[1]*100,"%", "-", Pcutoffs[2]*100, "%)=",avg.MTDIV,sep=""))
cat("\n")
cat(paste("Metric V : Proportion of trials with MTD within over dose region (>= ",  Pcutoffs[2]*100, "%)=",avg.MTDV,sep=""))
cat("\n")
cat(paste("Metric VI : Proportion of trials with MTD within under dose region (< ", Pcutoffs[1]*100,"%)=",avg.MTDVI,sep=""))
cat("\n")
cat(paste("Stopped : ",stopd,sep=""))
cat("\n")
cat(paste("Average sample size : ",avsamp,sep=""))
cat("\n")
cat("Distribution of DLT")
cat("\n")
print(aveMTD.dose)
cat("\n")
cat("Average Number of patients in each dose  combination: ")
cat("\n")
print(avesubj.dose)
cat("\n")
cat("Average Number of DLT in each dose  combination:")
cat("\n")
print(aveDLT.dose)
sink()

sim.summ <- list(agvsubj.target=agvsubj.doseI,
                 agvsubj.over = agvsubj.doseII,
                 agvsubj.under  =  agvsubj.doseIII,
                 MTD.target = avg.MTDIV,
                 MTD.over   = avg.MTDV,
                 MTD.under  = avg.MTDVI,
                 trial.stopped = stopd,
                 average.sample = avsamp,
                 avesubj.dose = avesubj.dose,
                 aveDLT.dose = aveDLT.dose,
                 aveMTD.dose =aveMTD.dose
                 )

}
