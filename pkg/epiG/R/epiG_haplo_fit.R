#
#     Description of this R script:
#     TODO
#
#     Intended for use with R.
#     Copyright (C) 2014 Martin Vincent
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http:#www.gnu.org/licenses/>
#

#' epiG
#' 
#' @param filename 
#' @param refGenom_filename 
#' @param altGenom_filename 
#' @param refname 
#' @param start 
#' @param end 
#' @param max_threads 
#' @param config 
#' @return fitted model
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_haplo_fit_filename
#' @useDynLib epiG r_epiG_haplo_fit_filename_chunks
epiG <- function(filename, refname, start, end, max_threads = 8L, config, refGenom_filename = config$ref.filename, altGenom_filename = config$alt.filename) {
	
	if(start >= end) {
		stop("end must be larger than start")
	}
	
	if(config$chunk.method == "none") {
		res <- .Call(r_epiG_haplo_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)
		
	} else if(config$chunk.method == "bases") {
		res <- .Call(r_epiG_haplo_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(config$chunk_size), config)
		
	} else if(config$chunk.method == "reads")  {
		s <- compute_chunk_positions(filename, refname, as.integer(start), as.integer(end), as.integer(config$chunk_size))
		
		if(length(s) >= 2) {
			chunks_start <- s[1:(length(s)-1)]
			chunks_end <- s[2:length(s)]-1L
			chunks_start[1] <- start
			refnames <- as.list(rep(refname, length(chunks_start)))
			
			res <- .Call(r_epiG_haplo_fit_filename_chunks, filename, refGenom_filename, altGenom_filename, refnames,  as.integer(chunks_start),  as.integer(chunks_end), as.integer(max_threads), config)
		
		} else {
			res <- .Call(r_epiG_haplo_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)
		}
		
	} else {
		stop("Unknown chunk method")
	}
	
	n_chunks <- res$number_of_chunks
	
	res.chunks <- list()

	for(i in 1:n_chunks) {
		
		res.chunks[[i]] <- list()
		
		res.chunks[[i]]$epiG_version <- packageVersion("epiG")
		res.chunks[[i]]$date <- date()
		
		res.chunks[[i]]$config <- config
		res.chunks[[i]]$filename <- filename
		res.chunks[[i]]$refname <- refname
		
		res.chunks[[i]]$offset <- res$chunks_start[i]
		
		res.chunks[[i]]$read_ids <- lapply(res$read_id[[i]], function(x) x + 1)
		res.chunks[[i]]$haplotype$chain <- as.integer(factor(res$haplotype[[i]]))
		res.chunks[[i]]$haplotype$start <- res$chain_start[[i]]
		res.chunks[[i]]$haplotype$end <- res$chain_end[[i]]
		
		if(length(res.chunks[[i]]$haplotype$chain) > 0) {
			res.chunks[[i]]$strands <- sapply(1:max(res.chunks[[i]]$haplotype$chain), function(chain) res$strands[[i]][res.chunks[[i]]$haplotype$chain == chain][1])
		} else {
			res.chunks[[i]]$strands <- as.character()
		}
		
		res.chunks[[i]]$strands <- factor(res.chunks[[i]]$strands)
		levels(res.chunks[[i]]$strands) <- c("fwd", "rev")
		
		res.chunks[[i]]$genotypes <-  lapply(res$genotypes[[i]], function(x) as.integer(x + 1))
				
		res.chunks[[i]]$length <- res$chunks_end[i] - res.chunks[[i]]$offset + 1
		
		class(res.chunks[[i]]) <- "epiG"
	}
	
	if(n_chunks == 1) {
		return(res.chunks[[1]])
	}	
	
	class(res.chunks) <- c("epiG", "chunks")
	return(res.chunks)
}

#' epiG.chunks
#' 
#' @param filename 
#' @param refGenom_filename 
#' @param altGenom_filename 
#' @param refname 
#' @param chunks_start 
#' @param chunks_end 
#' @param max_threads 
#' @param config 
#' @return fitted models
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_haplo_fit_filename_chunks
epiG.chunks <- function(filename, refnames, chunks_start, chunks_end, max_threads = 8L, config, refGenom_filename = config$ref.filename, altGenom_filename = config$alt.filename) {
	
	refnames <- as.list(refnames)
	
	if(length(refnames) != length(chunks_start)) {
		stop("length of refnames not equal to the number of chunks")
	}
	
	res <- .Call(r_epiG_haplo_fit_filename_chunks, filename, refGenom_filename, altGenom_filename, refnames,  as.integer(chunks_start),  as.integer(chunks_end), as.integer(max_threads), config)
	
	n_chunks <- res$number_of_chunks
	
	res.chunks <- list()
			
	for(i in 1:n_chunks) {

		res.chunks[[i]] <- list()
		
		res.chunks[[i]]$epiG_version <- packageVersion("epiG")
		res.chunks[[i]]$date <- date()
		
		res.chunks[[i]]$config <- config
		res.chunks[[i]]$filename <- filename
		res.chunks[[i]]$refname <- refnames[[i]]
		
		res.chunks[[i]]$offset <- res$chunks_start[i]
		
		res.chunks[[i]]$read_ids <- lapply(res$read_id[[i]], function(x) x + 1)
		res.chunks[[i]]$haplotype$chain <- as.integer(factor(res$haplotype[[i]]))
		res.chunks[[i]]$haplotype$start <- res$chain_start[[i]]
		res.chunks[[i]]$haplotype$end <- res$chain_end[[i]]
		
		if(length(res.chunks[[i]]$haplotype$chain) > 0) {
			res.chunks[[i]]$strands <- sapply(1:max(res.chunks[[i]]$haplotype$chain), function(chain) res$strands[[i]][res.chunks[[i]]$haplotype$chain == chain][1])
		} else {
			res.chunks[[i]]$strands <- as.character()
		}
		
		res.chunks[[i]]$strands <- factor(res.chunks[[i]]$strands)
		levels(res.chunks[[i]]$strands) <- c("fwd", "rev")
		
		res.chunks[[i]]$genotypes <-  lapply(res$genotypes[[i]], function(x) as.integer(x + 1))
		
		res.chunks[[i]]$length <- as.integer(res$chunks_end[i] - res.chunks[[i]]$offset + 1)
		
		class(res.chunks[[i]]) <- "epiG"
	}
	
	if(n_chunks == 1) {
		return(res.chunks[[1]])
	}	
	
	class(res.chunks) <- c("epiG", "chunks")
	return(res.chunks)
}


#' symbols
#' 
#' @param g 
#' @return symbols
#' 
#' @author martin
#' @export
symbols <- function(g) {
	return(unlist(sapply(g, function(x) c("N", "C", "G", "A", "T", "c", "g")[x+1])))
}

