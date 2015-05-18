\name{create_genotype_prior_alt}
\alias{create_genotype_prior_alt}
\title{create_genotype_prior_alt}
\usage{
  create_genotype_prior_alt(scale = 0.5,
    R = exp(5 * scale), RA = exp(4.5 * scale),
    A = exp(4.5 * scale), RAB = exp(4 * scale),
    RABB = exp(3 * scale), RB = exp(3 * scale),
    AB = exp(3 * scale), RBB = exp(2 * scale),
    ABB = exp(2 * scale), B = exp(scale), BB = exp(scale))
}
\arguments{
  \item{scale}{}

  \item{R}{}

  \item{RA}{}

  \item{A}{}

  \item{RB}{}

  \item{AB}{}

  \item{B}{}

  \item{RAB}{}

  \item{RBB}{}

  \item{ABB}{}

  \item{RABB}{}

  \item{BB}{}
}
\value{
  prior
}
\description{
  create_genotype_prior_alt
}
\author{
  Martin Vincent
}

