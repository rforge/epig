#
#     Description of this R script:
#     TODO
#
#     Intended for use with R.
#     Copyright (C) 2013 Martin Vincent
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' Compute chunk positions
#' 
#' @param filename 
#' @param refname 
#' @param start 
#' @param end 
#' @param chunk_size 
#' @return chunk positions
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_compute_chunk_positions
compute_chunk_positions <- function(filename, refname, start, end, chunk_size) {

    pos <- .Call(r_epiG_compute_chunk_positions, filename, refname, as.integer(start), as.integer(end), as.integer(chunk_size))

    return(pos);
}

#' fetch_reads_info
#' 
#' @param filename 
#' @param refname 
#' @param start 
#' @param end 
#' @return info 
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_fetch_reads_info
fetch_reads_info <- function(filename, refname, start, end) {

        res <- .Call(r_epiG_fetch_reads_info, filename, refname, as.integer(start), as.integer(end))

        info <- list()
        info$start <- as.integer(res$position)
        info$end <- as.integer(res$position + res$length - 1)
        info$length <- as.integer(res$length)
        info$nread <- length(res$position)

        return(info)
}

#' fetch_reads
#' 
#' @param object 
#' @return epiG model
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_fetch_reads
fetch.reads <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		reads <- .Call(r_epiG_fetch_reads, object$filename, object$refname, start(object), end(object))
		
		reads$positions <- reads$positions - object$offset
		object$reads <- reads
	
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		
		object <- lapply(object, fetch.reads)
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}
	
	return(object)
}

#' fetch_ref
#' 
#' @param object 
#' @param filename 
#' @return epiG model
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_read_fasta
fetch_ref <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		object$ref <- .Call(r_epiG_read_fasta, object$config$ref.filename, object$refname, start(object), length(object))
		
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		object <- lapply(object, function(x) fetch_ref(x))
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}
	
	return(object)
}

#' Read fasta
#' 
#' @param filename 
#' @param refname 
#' @param start 
#' @param len 
#' @return ??
#' 
#' @author martin
#' @export
read.fasta <- function(filename, refname, start, len) {
	return(.Call(r_epiG_read_fasta, filename, refname, as.integer(start), as.integer(len)))
}

#' fetch_alt
#' 
#' @param object 
#' @param filename 
#' @return epiG model
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_read_fasta
fetch_alt <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		object$alt <- .Call(r_epiG_read_fasta, object$config$alt.filename, object$refname, start(object), length(object))
	
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		object <- lapply(object, function(x) fetch_alt(x))
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}

	return(object)
}


