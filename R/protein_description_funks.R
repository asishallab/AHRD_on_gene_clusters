#' Retains those argument \code{prot.descs} that do not match any regular
#' expression in the argument \code{blacklist.regexs}.
#'
#' @param prot.descs A character vector holding protein descriptions to be
#' submitted to blacklisting.
#' @param blacklist.regexs A character vector of regular expression to be used
#' in the process of blacklisting. Default is
#' \code{getOption('AHRD.prot.desc.blacklist',
#' readLines(file.path(path.package('AHRD.on.gene.clusters'),
#' 'blacklist_descline.txt')))}
#' @param ... optional additional parameters passed to function
#' \code{AHRD.on.gene.clusters::blacklist}
#'
#' @export
#' @return A subset of \code{prot.descs} of which no element matched any
#' regular expression in \code{blacklist.regexs}.
blacklistProtDescs <- function(prot.descs, blacklist.regexs = getOption("AHRD.prot.desc.blacklist", 
    readLines(file.path(path.package("AHRD.on.gene.clusters"), "blacklist_descline.txt"))), 
    ...) {
    blacklist(prot.descs, blacklist.regexs, ...)
}

#' Retains those argument \code{tokens} that do not match any regular
#' expression in the argument \code{blacklist.regexs}.
#'
#' @param tokens A character vector holding protein descriptions to be
#' submitted to blacklisting.
#' @param blacklist.regexs A character vector of regular expression to be used
#' in the process of blacklisting. Default is
#' \code{getOption('AHRD.on.gene.clusters',
#' readLines(file.path(path.package('AHRD.on.gene.clusters'),
#' 'blacklist_token.txt')))}
#' @param ... optional additional parameters passed to function
#' \code{AHRD.on.gene.clusters::blacklist}
#'
#' @export
#' @return A subset of \code{tokens} of which no element matched any regular
#' expression in \code{blacklist.regexs}.
blacklistTokens <- function(tokens, blacklist.regexs = getOption("AHRD.token.blacklist", 
    readLines(file.path(path.package("AHRD.on.gene.clusters"), "blacklist_token.txt"))), 
    ...) {
    blacklist(tokens, blacklist.regexs, ...)
}

#' Function returns all those elements of \code{strs} that do not match any
#' regular expression in \code{blacklist.regexs}.
#'
#' @param strs A character vector of strings to be submitted to blacklisting.
#' @param blacklist.regexs A character vector of regular expressions to be used
#' in blacklisting.
#' @param ... optional additional arguments passed to \code{base::grepl}.
#'
#' @export
#' @return A subset of \code{strs} of which no element matches any regular
#' expression in argument \code{blacklist.regexs}.
blacklist <- function(strs, blacklist.regexs, ...) {
    strs[!Reduce(`|`, lapply(blacklist.regexs, function(b.regex) grepl(b.regex, 
        strs, ...)))]
}

#' Within each element of argument \code{prot.descs} all substrings matching
#' any of the regular expressions in argument \code{filter.regexs} is deleted.
#' Finally leading and trailing white spaces are also abolished. 
#'
#' @param prot.descs A character vector of protein descriptions to be filtered.
#' @param filter.regexs A character vector of regular expressions to be used
#' within filtering. Default is \code{getOption('AHRD.prot.desc.filter',
#' readLines(file.path(path.package('AHRD.on.gene.clusters'),
#' 'filter_descline.txt')))}
#' @param ... optional additional parameters passed to function
#' \code{AHRD.on.gene.clusters::filter}.
#'
#' @export
#' @return A modified version of \code{prot.descs} in which all matches to any
#' regular expression in \code{filter.regexs}, along with leading or trailing
#' white spaces has been deleted.
filterProtDescs <- function(prot.descs, filter.regexs = getOption("AHRD.prot.desc.filter", 
    readLines(file.path(path.package("AHRD.on.gene.clusters"), "filter_descline.txt"))), 
    ...) {
    filter(prot.descs, filter.regexs, ...)
}

#' Within each element of argument \code{tokens} all substrings matching
#' any of the regular expressions in argument \code{filter.regexs} is deleted.
#' Finally leading and trailing white spaces are also abolished. 
#'
#' @param tokens A character vector of protein descriptions to be filtered.
#' @param filter.regexs A character vector of regular expressions to be used
#' within filtering. Default is \code{getOption('AHRD.prot.desc.filter',
#' readLines(file.path(path.package('AHRD.on.gene.clusters'),
#' 'filter_token.txt')))}
#' @param ... optional additional parameters passed to function
#' \code{AHRD.on.gene.clusters::filter}.
#'
#' @export
#' @return A modified version of \code{tokens} in which all matches to any
#' regular expression in \code{filter.regexs}, along with leading or trailing
#' white spaces has been deleted.
filterTokens <- function(tokens, filter.regexs = getOption("AHRD.token.filter", 
    readLines(file.path(path.package("AHRD.on.gene.clusters"), "filter_token.txt"))), 
    ...) {
    filter(tokens, filter.regexs, ...)
}

#' Within each element of argument \code{strs} all substrings matching
#' any of the regular expressions in argument \code{filter.regexs} is deleted.
#' Finally leading and trailing white spaces are also abolished. 
#'
#' @param strs A character vector of protein descriptions to be filtered.
#' @param filter.regexs A character vector of regular expressions to be used
#' within filtering. 
#' @param ... optional additional parameters passed to function
#' \code{base::gsub}.
#'
#' @export
#' @return A modified version of \code{strs} in which all matches to any
#' regular expression in \code{filter.regexs}, along with leading or trailing
#' white spaces has been deleted.
filter <- function(strs, filter.regexs, ...) {
    for (flt in filter.regexs) {
        strs <- gsub(flt, "", strs, ...)
    }
    sub("^\\s+", "", sub("\\s+$", "", strs))
}

#' Parses a FASTA file containg all or more proteins present in the gene
#' families that need annotations. The result is a database of protein
#' descriptions mapped to long and short accessions (IDs). A first filtering
#' step of the protein descriptions is done.
#'
#' @param path.2.fasta The valid path to the FASTA file containing the
#' references proteins and their descriptions.
#' @param prot.id.regex A regular expression expected to hold a single match
#' group. It is used to extract the protein IDs from the protein description
#' lines. Defaut is \code{getOption('AHRD.prot.id.regex', '^(\\S+)')}.
#' @param prot.short.id.regex In the case of some public protein databases like
#' Uniprot the short protein IDs are embedded within their longer form. This
#' optional regular expression is used to extract the short ID. Default is
#' \code{options('AHRD.short.prot.id.regex', '^[^|]+\\|([^|]+)\\|.*$')}. Set to
#' \code{NULL}, \code{NA}, or \code{c()} to skip extraction of short protein
#' IDs.
#' @param filter.prot.desc.funk A function used to filter the protein
#' descriptions. Default is \code{getOption('AHRD.filter.prot.desc.funk',
#' filterProtDescs)}.
#'
#' @export
#' @return A data.frame with the following columns: \code{ID} holds the protein
#' IDs, \code{SHORT.ID} holds the protein short IDs or \code{NA},
#' \code{original.description} holds the infiltered original protein
#' descriptions, and \code{filtered.description} holds the filtered protein
#' descriptions.
readProteinDescriptionDb <- function(path.2.fasta, prot.id.regex = getOption("AHRD.prot.id.regex", 
    "^(\\S+)"), prot.short.id.regex = getOption("AHRD.short.prot.id.regex", 
    "^[^|]+\\|([^|]+)\\|.*$"), filter.prot.desc.funk = getOption("AHRD.filter.prot.desc.funk", 
    filterProtDescs)) {
    prot.descs <- sub("^>", "", system(paste("grep '^>'", path.2.fasta), 
        intern = TRUE))
    prot.ids <- regmatches(prot.descs, regexpr(prot.id.regex, prot.descs))
    prot.short.ids <- if (!is.null(prot.short.id.regex) && !is.na(prot.short.id.regex) && 
        length(prot.short.id.regex) > 0) {
        regmatches(prot.descs, regexpr(prot.short.id.regex, prot.descs))
    } else NA
    data.frame(ID = prot.ids, SHORT.ID = prot.short.ids, original.description = prot.descs, 
        filtered.description = filter.prot.desc.funk(prot.descs), stringsAsFactors = FALSE)
}

#' Generates a human readable description (HRD) for the argument
#' \code{prot.fam} based on the protein descriptions found for the genes
#' contained within the family. Candidate descriptions are filtered using
#' \code{AHRD.on.gene.clusters::filterTokens} and blacklisted with
#' \code{AHRD.on.gene.clusters::blacklistProtDescs}. The most frequent of the
#' retained descriptions are assigned as the clusters HRD. Maximum frequency
#' and Shannon Entropy are also returned as a measure of confidence.
#'
#' @param prot.fam A character vector of gene accessions (IDs) forming the
#' argument protein family to annotate.
#' @param prot.desc.db A data.frame with protein descriptions for the family's
#' members. Typically the result of invoking
#' \code{AHRD.on.gene.clusters::readProteinDescriptionDb}.
#' @param prot.id.col A string or integer identifying the column of
#' \code{prot.desc.db} in which to find the protein IDs. Default is
#' \code{"SHORT.ID"}.
#' @param prot.desc.col A string or integer identifying the column of
#' \code{prot.desc.db} in which to find the protein descriptions. Default is
#' \code{"filtered.description"}.
#'
#' @export
#' @return A list with the following key-value pairs: \code{descriptions} a
#' character vector of the most frequent filtered and blacklisted protein
#' descriptions, \code{max.frequency} the maximum frequency used to identify
#' the former protein descriptions, \code{entropy} the Shannon Entropy of the
#' filtered and blacklisted protein descriptions found for this protein family
#' \code{prot.fam}. \code{NA} is returned if no protein descriptions are
#' retained after filtering and blacklisting.
annotateClusterWithProteinDescription <- function(prot.fam, prot.desc.db, 
    prot.id.col = "SHORT.ID", prot.desc.col = "filtered.description") {
    pf.descs <- blacklistProtDescs(filterTokens(prot.desc.db[which(prot.desc.col[, 
        prot.id.col] %in% prot.fam), prot.desc.col]))
    if (!is.null(pf.descs) && length(pf.descs) > 0) {
        pf.descs.tbl <- table(pf.descs)
        pf.descs.entropy <- entropy.empirical(pf.descs.tbl)
        pf.descs.tbl.max <- max(pf.descs.tbl)
        list(descriptions = names(pf.descs.tbl[which(pf.descs.tbl == pf.descs.tbl.max)]), 
            max.frequency = pf.descs.tbl.max, entropy = pf.descs.entropy)
    } else NA
}
