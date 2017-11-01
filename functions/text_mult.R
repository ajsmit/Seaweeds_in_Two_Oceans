text.mult <- function (x, labels, display, choices = c(1, 2), rescale = TRUE,
          fill = 0.75, ...)
{
  X <- if (is.matrix(x)) {
    nc <- NCOL(x)
    if (nc != 2L) {
      stop("A 2-column matrix of coordinates is required & not supplied.")
    }
    x
  }
  else {
    if (inherits(x, "envfit")) {
      scores(x, display = "vectors", ...)[, 1:2]
    }
    else {
      scores(x, display = display, choices = choices,
             ...)
    }
    if (!rescale) {
      warning("Extracted scores usually need rescaling but you set 'rescale = FALSE'.\nConsider using 'rescale = TRUE', the default.")
    }
  }
  if (rescale) {
    mul <- ordiArrowMul(X, fill = fill)
    X <- X * mul
  }
  if (missing(labels)) {
    rnames <- rownames(X)
    labels <- if (is.null(rnames)) {
      paste("V", seq_len(NROW(X)))
    }
    else {
      rnames
    }
  }
  w <- strwidth(labels, units = "figure")
  h <- strheight(labels, units = "figure")
  b <- X[, 2]/X[, 1]
  off <- cbind(sign(X[, 1]) * (w/2 + h/4), 0.75 * h * sign(X[,
                                                             2]))
  for (i in seq_len(nrow(X))) {
    move <- off[i, 2]/b[i]
    if (is.finite(move) && abs(move) <= abs(off[i, 1]))
      off[i, 1] <- move
    else {
      move <- b[i] * off[i, 1]
      off[i, 2] <- move
    }
  }
  off + X
}
