##' Euler-Bernoulli bending with 3 parameters. EI is allowed to vary
##' linearly along the length of the beam.
##'
##' @title Euler-Bernoulli bending with 3 parameters
##' @param parms vector. 3 parameters for fitting (C1, C2, Offset). C1
##'   is the intercept, C2 is the slope, and the Offset allows the
##'   zero point to not be at \code{(0,0)}.
##' @param optim bernoulli. If \code{TRUE}, only the sum of squares is
##'   returned
##' @param n integer. Number of divisions for the theoretical beam.
##'   Defaults to 1000.
##' @param beam_length integer. Beam length in m.
##' @param data data.frame. Observed data.
##' @return If \code{optim = FALSE}, a list containing the theoretical
##'   bone, the observed data, and the sum of squares.
##' @author Kevin Middleton
##' @export
##' 
bending_fit_3 <- function(parms, optim = TRUE, n = 1000,
                          beam_length,
                          data){
  dat <- data
  
  # Set up constants
  if (length(parms) != 3){
    stop("Must supply 3 parameters.")
  } else {
    C1 <- parms[1]
    C2 <- parms[2]
    Offset <- parms[3]
  }
  
  if (!optim){
    message(paste0("C1: ", C1, "; C2: ", C2, "; Offset: ", Offset))
  }
  
  # Set up x's
  index <- 0:n
  x <- index * (beam_length / n)
  
  # Calculate EI
  C1_comp <- rep(C1, length(x))
  if (is.na(C2)){
    C2_comp <- rep(0, length(x))
  } else {
    C2_comp <- C2 * x
  }
  EI <- C1_comp + C2_comp
  
  # Calculate u''/P
  upp_P <- (x - beam_length) / EI
  
  # Calcuate u'/P
  up_P <- rep(0, length(x))
  for (i in 2:length(x)){
    up_P[i] <- up_P[i - 1] + 
      (x[i] - x[i - 1]) * (upp_P[i] + upp_P[i - 1]) / 2
  }
  
  # Calculate u/P
  u_P <- rep(NA, length(x))
  u_P[1] <- Offset
  for (i in 2:length(x)){
    u_P[i] <- u_P[i - 1] +
      (x[i] - x[i - 1]) * (up_P[i] + up_P[i - 1]) / 2
  }
  
  # Lookup predicted u/P for each observed x
  dat$pred <- NA
  for (i in 1:nrow(dat)){
    ind <- which(abs(x - dat$x[i]) == 
                   min(abs(x - dat$x[i])))
    dat$pred[i] <- u_P[ind]
  }
  
  # Squared difference
  dat$weight <- 1
  dat$error <- dat$weight * (dat$pred - dat$u_P)^2
  
  test_mat <- cbind(index, x, EI, upp_P, up_P, u_P)
  test_mat <- as.data.frame(test_mat)
  
  value <- sum(dat$error) * 100000000
  if (optim){
    return(value)
  } else {
    return(list(test_mat = test_mat,
                dat = dat,
                value = value))
  }
}
