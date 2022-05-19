BUGS.Table2Array.bord <- function (BUGSTable, Labels = NULL, Dim = NULL) 
{
  if (is.null(Labels) & !is.null(Dim)) {
    Labels <- list()
    for (j in 1:length(Dim)) {
      Labels[[j]] <- paste("X", j, ".", 1:Dim[j], sep = "")
    }
  }
  if (!is.null(Labels) & is.null(Dim)) {
    Dim <- sapply(Labels, length)
  }
  DimInd <- Dim
  Labels[[length(Labels) + 1]] <- dimnames(BUGSTable)[[2]]
  Dim <- c(Dim, ncol(BUGSTable))
  varName <- sub("(^[^\\[]*)(.*)", "\\1", rownames(BUGSTable)[1])
  ind <- paste(varName, "[", do.call(paste, c(do.call(expand.grid, 
                                                      lapply(DimInd, seq)), list(sep = ","))), "]", sep = "")
  # BUGSTable <- BUGSTable[ind, ]
  lenLabels <- length(Labels)
  BUGSArray <- array(BUGSTable, Dim)
  dimnames(BUGSArray) <- Labels
  return(BUGSArray = BUGSArray)
}