##' Euler-Bernoulli bending with 2 parameters. EI is constant along
##' the length of the beam.
##'
##' @title Euler-Bernoulli bending with 2 parameters
##' @param parms vector. 2 parameters for fitting (C1, Offset). C1 is
##'   the intercept for EI, and the Offset allows the zero point to
##'   not be at \code{(0,0)}.
##' @param optim bernoulli. If \code{TRUE}, only the sum of squares is
##'   returned
##' @param n integer. Number of divisions for the theoretical beam.
##'   Defaults to 1000.
##' @param beam_length integer. Beam length in m.
##' @param data data.frame. Observed data.
##' @return If \code{optim = FALSE}, a list containing the theoretical
##'   bone (\code{test_mat}), the observed data (\code{dat}), and the
##'   sum of squares (\code{value}).
##' @author Kevin Middleton
##' @export
##' @examples
##' # Pixels per cm
##' scale_pix <- 107
##' 
##' # Length in m
##' beam_length <- 0.03
##' 
##' glosso$Trial <- as.factor(glosso$Trial)
##' 
##' glosso$x <- glosso$x_pix/scale_pix/100
##' glosso$u <- glosso$u_pix/scale_pix/100
##' 
##' fit_2_optim <- nmkb(c(0.001, 0),
##'                     bending_fit_2,
##'                     lower = rep(-1, 2),
##'                     upper = rep(1, 2),
##'                     beam_length = beam_length,
##'                     data = glosso)
##' fit_2_optim$value
##' 
##' fit_2 <- bending_fit_2(parms = c(fit_2_optim$par[1],
##'                                  fit_2_optim$par[2]),
##'                                  optim = FALSE,
##'                                  beam_length = beam_length,
##'                                  data = glosso)
##' fit_2$value
##' 
##' test_mat <- fit_2$test_mat
##' dat2 <- fit_2$dat
##' 
##' p1 <- ggplot(dat2, aes(x = x, y = u_P)) +
##'   geom_point(aes(color = Trial), size = 3) +
##'   geom_line(data = test_mat, aes(x = x, y = u_P),
##'             color = I("red")) +
##'   geom_point(data = data.frame(x = 0, y = 0), aes(x, y))
##' p2 <- ggplot(test_mat, aes(x = x, y = EI)) + geom_line()
##' grid.arrange(p1, p2, nrow = 2)
##' c(min(test_mat$EI), max(test_mat$EI))
##' 
bending_fit_2 <- function(parms, optim = TRUE, n = 1000,
                          beam_length,
                          data){
  dat <- check_data(data)
  
  # Set up constants
  if (length(parms) != 2){
    stop("Must supply 2 parameters.")
  } else {
    C1 <- parms[1]
    Offset <- parms[2]
  }
  
  if (!optim){
    message(paste0("C1: ", signif(C1, 5),
                   "; Offset: ", signif(Offset, 5)))
  }
  
  # Set up x's
  index <- 0:n
  x <- index * (beam_length / n)
  
  # Calculate EI
  C1_comp <- rep(C1, length(x))
  C2_comp <- rep(0, length(x))
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
