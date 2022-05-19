Plot.P <- function (Est, Lower = "2.5%", Upper = "97.5%", Point = "mean",
          RowValues = 1:nrow(Est), RowNames = rownames(Est), ColNames = colnames(Est),
          xlim = NULL, ylim = NULL, xlab = "Doses", ylab = "DLT rate", title = "",
          bw = T, blackwhite = "black", color = "red", mar = c(5, 4, 4, 1), lty = 1,
          pch = 15, lwd = 1, las = 1, axes = F, Outfile = NULL)
{
  if (is.vector(Est)) {
    cnames = names(Est)
    Est = matrix(Est, nrow = 1)
    colnames(Est) = cnames
  }
  mar.old = graphics::par()$mar
  graphics::par(mar = mar)
  nnn <- dim(Est)
  if (!is.element(Lower, ColNames))
    stop("Lower is not a column in Est")
  if (!is.element(Upper, ColNames))
    stop("Upper is not a column in Est")
  if (!is.element(Point, ColNames))
    stop("Point is not a column in Est")
  if (is.null(xlim))
    xlim = range(RowValues)
  if (is.null(ylim))
    ylim = range(Est[, c(Lower, Upper)])
  if (bw)
    col <- blackwhite
  else col <- color
  #if (!is.null(Outfile))
    #Plot2File(Outfile)
  graphics::plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, main = title,
       las = las, axes = axes)
  for (j in 1:nnn[1]) {
    graphics::lines(rep(RowValues[j], 2), Est[j, c(Lower, Upper)],
          lty = lty, col = col, lwd = lwd)
    graphics::points(RowValues[j], Est[j, Point], pch = pch, col = col)
  }
  if (!axes) {
    graphics::axis(1, at = RowValues, labels = RowNames)
    graphics::axis(2)
  }
  graphics::par(mar = mar.old)
}
