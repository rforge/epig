\name{epiG.algorithm.config}
\alias{epiG.algorithm.config}
\title{epiG.algorithm.config}
\usage{
  epiG.algorithm.config(max_iterations = 1e+05,
    prior = list(create_genotype_prior(), create_genotype_prior_alt()),
    model = create_bisulfite_model(),
    unique_threshold = log(1.3), allele_match_prob = 0.7,
    max_alleles = 4, occupancy_weight = 0.01,
    chunk_size = 1000, chunk_method = "reads",
    reads_hard_limit = 1500, verbose = TRUE)
}
\arguments{
  \item{max_iterations}{}

  \item{prior}{}

  \item{model}{}

  \item{allele_match_prob}{}

  \item{max_alleles}{}

  \item{occupancy_weight}{}

  \item{reads_hard_limit}{}

  \item{verbose}{}
}
\value{
  configuration
}
\description{
  epiG.algorithm.config
}
\author{
  martin
}

