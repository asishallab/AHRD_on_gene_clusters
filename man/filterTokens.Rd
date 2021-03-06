% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/protein_description_funks.R
\name{filterTokens}
\alias{filterTokens}
\title{Within each element of argument \code{tokens} all substrings matching
any of the regular expressions in argument \code{filter.regexs} is deleted.
Finally leading and trailing white spaces are also abolished.}
\usage{
filterTokens(tokens, filter.regexs = getOption("AHRD.token.filter",
  readLines(file.path(path.package("AhrdOnGeneClusters"), "filter_token.txt"))),
  ...)
}
\arguments{
\item{tokens}{A character vector of protein descriptions to be filtered.}

\item{filter.regexs}{A character vector of regular expressions to be used
within filtering. Default is \code{getOption('AHRD.prot.desc.filter',
readLines(file.path(path.package('AhrdOnGeneClusters'),
'filter_token.txt')))}}

\item{...}{optional additional parameters passed to function
\code{AhrdOnGeneClusters::filter}.}
}
\value{
A modified version of \code{tokens} in which all matches to any
regular expression in \code{filter.regexs}, along with leading or trailing
white spaces has been deleted.
}
\description{
Within each element of argument \code{tokens} all substrings matching
any of the regular expressions in argument \code{filter.regexs} is deleted.
Finally leading and trailing white spaces are also abolished.
}

