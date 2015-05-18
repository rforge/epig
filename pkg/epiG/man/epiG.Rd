\name{epiG}
\alias{epiG}
\title{Fit an epiG model}
\usage{
  epiG(filename, refname, start, end, max_threads = 8L,
    config, refGenom_filename = config$ref.filename,
    altGenom_filename = config$alt.filename)
}
\arguments{
  \item{filename}{BAM file with aligned and ordered reads}

  \item{refGenom_filename}{}

  \item{altGenom_filename}{}

  \item{refname}{}

  \item{start}{}

  \item{end}{}

  \item{max_threads}{}

  \item{config}{}
}
\value{
  fitted model
}
\description{
  Fit an epiG model
}
\author{
  Martin Vincent
}

