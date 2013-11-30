\name{bending_fit_2}
\alias{bending_fit_2}
\title{Euler-Bernoulli bending with 2 parameters}
\usage{
  bending_fit_2(parms, optim = TRUE, n = 1000, beam_length,
    data)
}
\arguments{
  \item{parms}{vector. 2 parameters for fitting (C1,
  Offset). C1 is the intercept for EI, and the Offset
  allows the zero point to not be at \code{(0,0)}.}

  \item{optim}{bernoulli. If \code{TRUE}, only the sum of
  squares is returned}

  \item{n}{integer. Number of divisions for the theoretical
  beam.  Defaults to 1000.}

  \item{beam_length}{integer. Beam length in m.}

  \item{data}{data.frame. Observed data.}
}
\value{
  If \code{optim = FALSE}, a list containing the
  theoretical bone (\code{test_mat}), the observed data
  (\code{dat}), and the sum of squares (\code{value}).
}
\description{
  Euler-Bernoulli bending with 2 parameters. EI is constant
  along the length of the beam.
}
\examples{
# Pixels per cm
scale_pix <- 107

# Length in m
beam_length <- 0.03

glosso$Trial <- as.factor(glosso$Trial)

glosso$x <- glosso$x_pix/scale_pix/100
glosso$u <- glosso$u_pix/scale_pix/100

fit_2_optim <- nmkb(c(0.001, 0),
                    bending_fit_2,
                    lower = rep(-1, 2),
                    upper = rep(1, 2),
                    beam_length = beam_length,
                    data = glosso)
fit_2_optim$value

fit_2 <- bending_fit_2(parms = c(fit_2_optim$par[1],
                                 fit_2_optim$par[2]),
                                 optim = FALSE,
                                 beam_length = beam_length,
                                 data = glosso)
fit_2$value

test_mat <- fit_2$test_mat
dat2 <- fit_2$dat

p1 <- ggplot(dat2, aes(x = x, y = u_P)) +
  geom_point(aes(color = Trial), size = 3) +
  geom_line(data = test_mat, aes(x = x, y = u_P),
            color = I("red")) +
  geom_point(data = data.frame(x = 0, y = 0), aes(x, y))
p2 <- ggplot(test_mat, aes(x = x, y = EI)) + geom_line()
grid.arrange(p1, p2, nrow = 2)
c(min(test_mat$EI), max(test_mat$EI))
}
\author{
  Kevin Middleton
}

