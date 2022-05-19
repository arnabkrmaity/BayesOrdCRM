Plot.Pcat <- function (Pcat, xlim = NULL, ylim = NULL, RowValues = 1:nrow(Pcat),
          RowNames = rownames(Pcat), ColNames = colnames(Pcat), crit = rep(1, ncol(Pcat)),
          xlab = "Doses", ylab = "Probability", title = "Interval Probabilities by Dose",
          bw = F, blackwhite = c("black", "darkgrey"), color = c("red", "lightgreen"),
          wid = NULL, mar = c(5, 4, 4, 1), Outfile = NULL)
{
  if (is.vector(Pcat)) {
    cnames = names(Pcat)
    Pcat = matrix(Pcat, nrow = 1)
    colnames(Pcat) = cnames
  }
  mar.old = graphics::par()$mar
  graphics::par(mar = mar)
  nnn <- dim(Pcat)
  if (is.null(wid))
    wid = min(diff(RowValues)) * 0.8
  if (is.null(xlim))
    xlim = range(RowValues)
  if (bw)
    col <- blackwhite
  else col <- color
  # if (!is.null(Outfile))
  #   Plot2File(Outfile)
  graphics::plot(xlim, c(0, nnn[2]), type = "n", axes = F, xlab = xlab,
       ylab = ylab, ylim = ylim, main = title)
  graphics::axis(1, at = RowValues, labels = RowNames)
  graphics::axis(2, at = seq(0, nnn[2], by = 0.25), labels = c(rep(c("0", "", "0.5", ""), nnn[2]), ""), cex = 0.7)
  graphics::abline(h = 0:nnn[2] - 1)
  Nint = nnn[2]
  if (length(crit) == Nint)
    crit = c(crit, 1)
  for (i in 1:nnn[1]) {
    for (j in 1:nnn[2]) {
      if (j > (Nint - 2))
        color = ifelse(Pcat[i, j] > crit[j] | Pcat[i,
                                                   Nint - 1] + Pcat[i, Nint] > crit[Nint + 1],
                       col[1], col[2])
      if (j <= (Nint - 2))
        color = ifelse(Pcat[i, j] > crit[j], col[1],
                       col[2])
      graphics::rect(RowValues[i] - wid/2, j - 1, RowValues[i] +
             wid/2, j - 1 + Pcat[i, j], col = color)
    }
  }
  text.x = rev(seq(min(RowValues), max(RowValues), length = nnn[2] +
                     2))
  for (j in 1:nnn[2]) {
    graphics::text(text.x[j + 1], j - 1/3, ColNames[j], col = "darkgrey",
         cex = 1.5)
    if (crit[j] < 1)
      graphics::abline(h = j - 1 + crit[j], lty = 3)
  }
  graphics::par(mar = mar.old)
}
