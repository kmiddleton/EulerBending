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
\author{
  Kevin Middleton
}

