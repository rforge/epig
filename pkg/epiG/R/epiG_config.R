# TODO: Add comment
# 
# Author: martin
###############################################################################


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
#' @return prior
#' 
#' @author Martin Vincent
#' @export
create_genotype_prior_ref <- function(scale = 0.5, R = exp(6*scale), RA = exp(5*scale), A = exp(4*scale), AA = exp(3*scale) , RAA = exp(2*scale), AAA = exp(scale), RAAA = 1)  {
	
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
#' @return prior
#' 
#' @author Martin Vincent
#' @export
create_genotype_prior_alt <- function(scale = 0.5, 
		R = exp(5*scale), 
		RA = exp(4.5*scale), 
		A = exp(4.5*scale), 
		RAB = exp(4*scale), 
		RABB = exp(3*scale), RB = exp(3*scale), AB = exp(3*scale),
		RBB = exp(2*scale), ABB = exp(2*scale), 
		B = exp(scale), BB = exp(scale))  {
	
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

#create_genotype_prior_alt <- function(scale = 0.5, R = exp(11*scale), RA = exp(10*scale), A = exp(9*scale), 
#		RB = exp(8*scale), AB = exp(7*scale), B = exp(6*scale),
#		RAB = exp(5*scale), RBB = exp(4*scale), ABB = exp(3*scale), 
#		RABB = exp(2*scale), BB = 1)  {
#	
#	prior_alt <- c(R = R, 
#			A1 = A, 
#			A2 = B, 
#			A3 = B, 
#			A1A2 = AB, 
#			RA3 = RB, 
#			A1A3 = AB,  
#			RA2 = RB,
#			RA1 = RA, 
#			A2A3 = BB, 
#			RA1A3 = RAB,
#			A1A2A3 = ABB,
#			RA2A3 = RBB, 
#			RA1A2 = RAB, 
#			RA1A2A3 = RABB)
#	
#	#Normalize
#	prior_alt <- prior_alt/sum(prior_alt)
#	
#	return(prior_alt)
#}

#' create_error_distributions
#' 
#' @param bisulfite_rate 
#' @param bisulfite_inap_rate 
#' @return bisulfite model
#' 
#' @author Martin Vincent
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
	bisulfite_model$fwd[2, 6] <- 0 #G g
	bisulfite_model$fwd[3, 3] <- 1 #A A
	bisulfite_model$fwd[4, 1] <- bisulfite_rate #T C
	bisulfite_model$fwd[4, 4] <- 1 #T T
	
	bisulfite_model$rev[1, 1] <- 1 #C C
	bisulfite_model$rev[1, 5] <- 0 #C c
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
#' @param bisulfite_rates 
#' @param bisulfite_inap_rate 
#' @param lambda 
#' @param Lmax 
#' @return ...
#' 
#' @author Martin Vincent
#' @export
create_bisulfite_model <- function(bisulfite_rates = 0.94, bisulfite_inap_rate = 0.06, lambda = 0.1, Lmax = 100) {
	
	#TODO postion decay
	#p <- exp_decay(lambda = lambda, Lmax = Lmax)
	#bisulfite_rates <- 1 - (1-bisulfite_rate)*Lmax*p
	
	bisulfite_rates <- rep(bisulfite_rates, Lmax)
	bisulfite_inap_rate <- rep(bisulfite_inap_rate, Lmax)
	
	model <- list()
	model$fwd <-lapply(1:Lmax, function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rate[i])$fwd) 
	model$rev <-lapply(1:Lmax, function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rate[i])$rev) 
	
	return(model)
}

#' Create a epiG configuration
#' 
#' @param ref.file 
#' @param alt.file 
#' @param max_iterations 
#' @param prior 
#' @param model 
#' @param sequence_quality_adjust 
#' @param haplo_prior modeled using a geometric distribution
#' @param ref_prior 
#' @param min_overlap_length 
#' @param chunk_size 
#' @param chunk_method 
#' @param reads_hard_limit 
#' @param verbose 
#' @return epiG configuration
#' 
#' @author martin
#' @export
epiG.algorithm.config <- function(ref.file, alt.file, max_iterations = 1e5, prior = list(create_genotype_prior_ref(), create_genotype_prior_alt()), model = create_bisulfite_model(), sequence_quality_adjust = 0.1, haplo_prior = 0.2, ref_prior = 0.9, min_overlap_length = 1, chunk_size = 1000, chunk_method = "reads", reads_hard_limit = 1500, verbose = TRUE) {
	
	#TODO check config valid
	# 1) chunk_size < reads_hard_limit
	
	config <- list()
	
	config$ref.filename <- ref.file
	
	config$alt.filename <- alt.file
	
	config$max_iterations <- as.integer(max_iterations)
	
	config$fwd_model <- model$fwd
	
	config$rev_model <- model$rev
	
	config$sequence_quality_adjust <- sequence_quality_adjust
	
	config$log_prior <- lapply(prior, log)
	
	config$chunk_size <- chunk_size
	
	config$chunk.method <- chunk_method
	
	config$reads_hard_limit <- as.integer(reads_hard_limit)
	
	config$haplo_prior <- haplo_prior
	
	config$ref_prior <- ref_prior
	
	config$min_overlap_length <- as.integer(min_overlap_length)
	
	config$verbose <- verbose
			
	class(config) <- "epiG.config"
	
	return(config)
}
