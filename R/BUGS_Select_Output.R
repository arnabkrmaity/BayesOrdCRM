BUGS.Select.Output <- function (parastr, BUGStable, cols = c("mean", "stdev", "2.5%", 
                                       "50%", "97.5%")) 
{
  vnames <- dimnames(BUGStable)[[1]]
  ns <- nchar(parastr)
  ix <- sapply(vnames, function(r) {
    e <- substring(paste(r, parastr, sep = ""), 1, ns) == 
      parastr
    if (nchar(r) > ns & substring(r, ns + 1, ns + 1) != "[") 
      e <- F
    return(e)
  })
  if (sum(ix) > 0) 
    return(BUGStable[ix, cols])
  else return(NULL)
}