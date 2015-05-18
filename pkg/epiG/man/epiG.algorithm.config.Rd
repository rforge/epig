\name{epiG.algorithm.config}
\alias{epiG.algorithm.config}
\title{Create a epiG configuration}
\usage{
  epiG.algorithm.config(ref.file, alt.file,
    max_iterations = 1e+05,
    prior = list(create_genotype_prior_ref(), create_genotype_prior_alt()),
    model = create_bisulfite_model(),
    sequence_quality_adjust = 0.1, haplo_prior = 0.2,
    ref_prior = 0.9, min_overlap_length = 1,
    chunk_size = 1000, chunk_method = "reads",
    reads_hard_limit = 1500, verbose = TRUE)
}
\arguments{
  \item{ref.file}{}

  \item{alt.file}{}

  \item{max_iterations}{}

  \item{prior}{}

  \item{model}{}

  \item{sequence_quality_adjust}{}

  \item{haplo_prior}{modeled using a geometric
  distribution}

  \item{ref_prior}{}

  \item{min_overlap_length}{}

  \item{chunk_size}{}

  \item{chunk_method}{}

  \item{reads_hard_limit}{}

  \item{verbose}{}
}
\value{
  epiG configuration
}
\description{
  Create a epiG configuration
}
\author{
  martin
}

