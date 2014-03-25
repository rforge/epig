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
#' @param chunk.size 
#' @param chunk.method 
#' @param config 
#' @returnType epiG
#' @return fitted model
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_fit_filename
#' @useDynLib epiG r_epiG_fit_filename_chunks

epiG <- function(filename, refGenom_filename, altGenom_filename, refname, start, end, max_threads = 8L, chunk.size = 500, chunk.method = "reads", config = msbs.standard.config) {
	
        if(chunk.method == "none") {
            res <- .Call(r_epiG_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)

        } else if(chunk.method == "bases") {
            res <- .Call(r_epiG_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(chunk.size), config)

            } else if(chunk.method == "reads")  {
            s <- compute_chunk_positions(filename, refname, start, end, chunk.size)

            if(length(s) >= 2) {
                chunks_start <- s[1:(length(s)-1)]
                chunks_end <- s[2:length(s)]-1L
                chunks_start[1] <- start
                res <- .Call(r_epiG_fit_filename_chunks, filename, refGenom_filename, altGenom_filename, refname,  as.integer(chunks_start),  as.integer(chunks_end), as.integer(max_threads), config)
            } else {
                res <- .Call(r_epiG_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)
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
		res.chunks[[i]]$ref_name <- refname
		res.chunks[[i]]$number_of_chunks <- 1 #TODO remove
		
		res.chunks[[i]]$offset <- res$chunks_start[i]
		res.chunks[[i]]$g <- res$g[[i]]
		res.chunks[[i]]$occupancy <- res$occupancy[[i]]
		
		res.chunks[[i]]$length <- ncol(res$g[[i]])
		
		res.chunks[[i]]$read_id <- list(lapply(res$read_id[[i]], function(x) x + 1))
		res.chunks[[i]]$strands <- list(res$strands[[i]]) 	
		res.chunks[[i]]$alleles <- list(res$alleles[[i]]) 	
		res.chunks[[i]]$unique <- list(res$unique[[i]])
		
		res.chunks[[i]]$chunks_start <- 0
		res.chunks[[i]]$chunks_end <- res$chunks_end[i] - res.chunks[[i]]$offset
		
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
#' @returnType epiG
#' @return fitted models
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_fit_filename_chunks
epiG.chunks <- function(filename, refGenom_filename, altGenom_filename, refname, chunks_start, chunks_end, max_threads = 8L, config = msbs.standard.config) {
	
	res <- .Call(r_epiG_fit_filename_chunks, filename, refGenom_filename, altGenom_filename, refname,  as.integer(chunks_start),  as.integer(chunks_end), as.integer(max_threads), config)
	
	n_chunks <- res$number_of_chunks
	
	res.chunks <- list()
	
	for(i in 1:n_chunks) {
		
		res.chunks[[i]] <- list()
		
		res.chunks[[i]]$epiG_version <- packageVersion("epiG")
                res.chunks[[i]]$date <- date()

		res.chunks[[i]]$config <- config
		res.chunks[[i]]$filename <- filename
		res.chunks[[i]]$ref_name <- refname
		res.chunks[[i]]$number_of_chunks <- 1 #TODO remove
		
		res.chunks[[i]]$offset <- res$chunks_start[i]
		res.chunks[[i]]$g <- res$g[[i]]
		res.chunks[[i]]$occupancy <- res$occupancy[[i]]
		
		res.chunks[[i]]$length <- ncol(res$g[[i]])
		
		res.chunks[[i]]$read_id <- list(lapply(res$read_id[[i]], function(x) x + 1))
		res.chunks[[i]]$strands <- list(res$strands[[i]]) 	
		res.chunks[[i]]$alleles <- list(res$alleles[[i]]) 	
		res.chunks[[i]]$unique <- list(res$unique[[i]])
		
		res.chunks[[i]]$chunks_start <- 0
		res.chunks[[i]]$chunks_end <- res$chunks_end[i] - res.chunks[[i]]$offset
		
		class(res.chunks[[i]]) <- "epiG"
	}
	
	if(n_chunks == 1) {
		return(res.chunks[[1]])
	}	
	
	class(res.chunks) <- c("epiG", "chunks")
	return(res.chunks)
}

# Genotype prior coding:
#
#	1 R 
#	2 A1 
# 	3 A2 
#	4 A3
#	5 A1 A2 
#	6 A3 R 
#	7 A1 A3 
#	8 A2 R 
#	9 A1 R 
#	10 A2 A3 
#	11 A1 A3 R 
#	12 A1 A2 A3 
#	13 A2 R A3 
#	14 A1 R A2 
#	15 A2 A1 R A3 

#' create_genotype_prior
#' 
#' @param scale 
#' @param R 
#' @param RA 
#' @param A 
#' @param AA 
#' @param RAA 
#' @param AAA 
#' @param RAAA 
#' @returnType vector
#' @return prior
#' 
#' @author martin
#' @export
create_genotype_prior <- function(scale = 0.5, R = exp(6*scale), RA = exp(5*scale), A = exp(4*scale), AA = exp(3*scale) , RAA = exp(2*scale), AAA = exp(scale), RAAA = 1)  {

	prior <- c(R = R, 
			A1 = A, 
			A2 = A, 
			A3 = A, 
			A1A2 = AA, 
			RA3 = RA, 
			A1A3 = AA,  
			RA2 = RA,
			RA1 = RA, 
			A2A3 = AA, 
			RA1A3 = RAA,
			A1A2A3 = AAA,
			RA2A3 = RAA, 
			RA1A2 = RAA, 
			RA1A2A3 = RAAA)
	
	#Normalize
	prior <- prior/sum(prior)
	
	return(prior)
}

#' create_genotype_prior_alt
#' 
#' @param scale 
#' @param R 
#' @param RA 
#' @param A 
#' @param RB 
#' @param AB 
#' @param B 
#' @param RAB 
#' @param RBB 
#' @param ABB 
#' @param RABB 
#' @param BB 
#' @returnType vector
#' @return prior
#' 
#' @author martin
#' @export
create_genotype_prior_alt <- function(scale = 0.5, R = exp(11*scale), RA = exp(10*scale), A = exp(9*scale), 
		RB = exp(8*scale), AB = exp(7*scale), B = exp(6*scale),
		RAB = exp(5*scale), RBB = exp(4*scale), ABB = exp(3*scale), 
		RABB = exp(2*scale), BB = 1)  {
	
	prior_alt <- c(R = R, 
			A1 = A, 
			A2 = B, 
			A3 = B, 
			A1A2 = AB, 
			RA3 = RB, 
			A1A3 = AB,  
			RA2 = RB,
			RA1 = RA, 
			A2A3 = BB, 
			RA1A3 = RAB,
			A1A2A3 = ABB,
			RA2A3 = RBB, 
			RA1A2 = RAB, 
			RA1A2A3 = RABB)
	
	#Normalize
	prior_alt <- prior_alt/sum(prior_alt)
	
	return(prior_alt)
}

#' create_error_distributions
#' 
#' @param bisulfite_rate 
#' @param bisulfite_inap_rate 
#' @returnType list
#' @return bisulfite model
#' 
#' @author martin
#' @export
create_error_distributions <- function(bisulfite_rate = 0.94, bisulfite_inap_rate = 0.06) {
	
	#TODO split up into 2 functions one for fwd model and one for rev model
	
	bisulfite_model <- list()
	
	bisulfite_model$fwd <- matrix(nrow = 4, ncol = 6)
	bisulfite_model$rev <- matrix(nrow = 4, ncol = 6)
	
	rownames(bisulfite_model$fwd) <- c('C', 'G', 'A', 'T')	
	rownames(bisulfite_model$rev) <- c('C', 'G', 'A', 'T')	
	colnames(bisulfite_model$rev) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	colnames(bisulfite_model$fwd) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	
	bisulfite_model$fwd[,] <- 0
	bisulfite_model$rev[,] <- 0
	
	bisulfite_model$fwd[1, 1] <- 1 - bisulfite_rate #C C
	bisulfite_model$fwd[1, 5] <- 1 - bisulfite_inap_rate #C c
	bisulfite_model$fwd[4, 5] <- bisulfite_inap_rate #T c
	bisulfite_model$fwd[2, 2] <- 1 #G G
	bisulfite_model$fwd[2, 6] <- 1 #G g
	bisulfite_model$fwd[3, 3] <- 1 #A A
	bisulfite_model$fwd[4, 1] <- bisulfite_rate #T C
	bisulfite_model$fwd[4, 4] <- 1 #T T

	bisulfite_model$rev[1, 1] <- 1 #C C
	bisulfite_model$rev[1, 5] <- 1 #C c
	bisulfite_model$rev[2, 2] <- 1 - bisulfite_rate #G G
	bisulfite_model$rev[2, 6] <- 1 - bisulfite_inap_rate #G g
	bisulfite_model$rev[3, 6] <- bisulfite_inap_rate #A g
	bisulfite_model$rev[3, 3] <- 1 #A A
	bisulfite_model$rev[3, 2] <- bisulfite_rate #A G
	bisulfite_model$rev[4, 4] <- 1 #T T
	
	return(bisulfite_model)
}

#' exp_decay
#' 
#' @param lambda 
#' @param Lmax 
#' @param x 
#' @returnType numeric
#' @return function values
#' 
#' @author martin
#' @export
exp_decay <- function(lambda = 0.1, Lmax = 100, x = 0:(Lmax-1)) {
	
	c <- (1-exp(-lambda))/(1-exp(-lambda*(Lmax)))

	return(c*exp(-lambda*x))
}

#' create_bisulfite_model
#' 
#' @param bisulfite_rate 
#' @param bisulfite_inap_rate 
#' @param lambda 
#' @param Lmax 
#' @returnType list
#' @return bisulfite model
#' 
#' @author martin
#' @export
create_bisulfite_model <- function(bisulfite_rate = 0.94, bisulfite_inap_rate = 0.06, lambda = 0.1, Lmax = 100) {
	
	p <- exp_decay(lambda = lambda, Lmax = Lmax)
	
	bisulfite_rates <- 1 - (1-bisulfite_rate)*Lmax*p
	
	bisulfite_inap_rate <- rep(bisulfite_inap_rate, Lmax)
	
	model <- list()
	model$fwd <-lapply(1:length(p), function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rate[i])$fwd) 
	model$rev <-lapply(1:length(p), function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rate[i])$rev) 
	
	return(model)
}

#' epiG.algorithm.config
#' 
#' @param max_iterations 
#' @param prior 
#' @param model 
#' @param allele_match_prob 
#' @param max_alleles 
#' @param occupancy_weight 
#' @param reads_hard_limit 
#' @param verbose 
#' @returnType algorithm.config
#' @return configuration
#' 
#' @author martin
#' @export
epiG.algorithm.config <- function(max_iterations = 1e5, prior = list(create_genotype_prior(), create_genotype_prior_alt()), model = create_bisulfite_model(), allele_match_prob = 0.7, max_alleles = 4, occupancy_weight = 0.01, reads_hard_limit = 750, verbose = TRUE) {
	
	config <- list()
	
	config$max_iterations <- as.integer(max_iterations)
	
	config$log_prior <- lapply(prior, log)
	
	config$max_alleles <-  as.integer(max_alleles)
	
	config$allele_occupancy_weight <- occupancy_weight
	
	config$allele_match_prob <- allele_match_prob
	
	config$fwd_model <- model$fwd
	
	config$rev_model <- model$rev
		
        config$reads_hard_limit <- as.integer(reads_hard_limit)

	config$verbose <- verbose
	
	class(config) <- "epiG.config"
	
	return(config)
}

epiG.standard.config <- epiG.algorithm.config()

#' symbols
#' 
#' @param g 
#' @returnType 
#' @return symbols
#' 
#' @author martin
#' @export
symbols <- function(g) {
	return(unlist(sapply(g, function(x) c("N", "C", "G", "A", "T", "c", "g")[x+1])))
}
