% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/protein_description_funks.R
\name{readProteinDescriptionDb}
\alias{readProteinDescriptionDb}
\title{Parses a FASTA file containg all or more proteins present in the gene
families that need annotations. The result is a database of protein
descriptions mapped to long and short accessions (IDs). A first filtering
step of the protein descriptions is done.}
\usage{
readProteinDescriptionDb(path.2.fasta,
  prot.id.regex = getOption("AHRD.prot.id.regex", "^(\\\\S+)"),
  prot.short.id.regex = getOption("AHRD.short.prot.id.regex",
  "^[^|]+\\\\|([^|]+)\\\\|.*$"),
  filter.prot.desc.funk = getOption("AHRD.filter.prot.desc.funk",
  filterProtDescs))
}
\arguments{
\item{path.2.fasta}{The valid path to the FASTA file containing the
references proteins and their descriptions.}

\item{prot.id.regex}{A regular expression expected to hold a single match
group. It is used to extract the protein IDs from the protein description
lines. Defaut is \code{getOption('AHRD.prot.id.regex', '^(\\S+)')}.}

\item{prot.short.id.regex}{In the case of some public protein databases like
Uniprot the short protein IDs are embedded within their longer form. This
optional regular expression is used to extract the short ID. Default is
\code{options('AHRD.short.prot.id.regex', '^[^|]+\\|([^|]+)\\|.*$')}. Set to
\code{NULL}, \code{NA}, or \code{c()} to skip extraction of short protein
IDs.}

\item{filter.prot.desc.funk}{A function used to filter the protein
descriptions. Default is \code{getOption('AHRD.filter.prot.desc.funk',
filterProtDescs)}.}
}
\value{
A data.frame with the following columns: \code{ID} holds the protein
IDs, \code{SHORT.ID} holds the protein short IDs or \code{NA},
\code{original.description} holds the infiltered original protein
descriptions, and \code{filtered.description} holds the filtered protein
descriptions.
}
\description{
Parses a FASTA file containg all or more proteins present in the gene
families that need annotations. The result is a database of protein
descriptions mapped to long and short accessions (IDs). A first filtering
step of the protein descriptions is done.
}

