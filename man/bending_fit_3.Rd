\name{bending_fit_3}
\alias{bending_fit_3}
\title{Euler-Bernoulli bending with 3 parameters}
\usage{
  bending_fit_3(parms, optim = TRUE, n = 1000, scale_pix,
    beam_length, data)
}
\arguments{
  \item{parms}{vector. 3 parameters for fitting (C1, C2,
  Offset)}

  \item{optim}{bernoulli. If \code{TRUE}, only the sum of
  squares is returned}

  \item{n}{integer. Number of divisions for the theoretical
  beam.  Defaults to 1000.}

  \item{scale_pix}{integer. Pixels per cm.}

  \item{beam_length}{integer. Beam length in m.}

  \item{data}{data.frame. Observed data.}
}
\value{
  If \code{optim = FALSE}, a list containing the
  theoretical bone, the observed data, and the sum of
  squares.
}
\description{
  Euler-Bernoulli bending with 3 parameters
}
\author{
  Kevin Middleton
}

