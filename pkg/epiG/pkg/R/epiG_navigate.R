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

#' subregion
#' 
#' @param object 
#' @param start_pos 
#' @param end_pos 
#' @param relative 
#' @returnType 
#' @return epiG model
#' 
#' @author martin
#' @export
subregion <- function(object, start_pos, end_pos, relative = TRUE) {
	
	if(relative == FALSE) {
		start_pos <- start_pos - object$offset 
		end_pos <- end_pos - object$offset		
	}
	
	new_model <- NULL
	new_model$filename <- object$filename
	new_model$config <- object$config
	new_model$offset <- as.integer(object$offset + start_pos)
	new_model$g <- object$g[,(start_pos+1):(end_pos+1)]
	new_model$occupancy <- object$occupancy[,(start_pos+1):(end_pos+1)]
	new_model$loglike <- object$loglike[(start_pos+1):(end_pos+1)]
	new_model$ref_name <- object$ref_name
	new_model$length <- ncol(new_model$g)

	if("ref" %in% names(object)) {
		new_model$ref <-  object$ref[(start_pos+1):(end_pos+1)]
	}

	if("alt" %in% names(object)) {
		new_model$alt <-  object$alt[(start_pos+1):(end_pos+1)]
	}
	
	## Compute new chunk regions
	new_model$chunks_start <- NULL
	new_model$chunks_end <- NULL
	
	new_model$read_id <- list()
	
	if("reads" %in% names(object)) {
		new_model$reads <- list()
	}
	
	for(i in 1:object$number_of_chunks) {
	
		if(object$chunks_start[i] > end_pos || object$chunks_end[i] < start_pos) {
			next
		}
			
		new_chunk_start <- max(start_pos, object$chunks_start[i])
		new_chunk_end <- min(end_pos, object$chunks_end[i])
				
		new_model$chunks_start <- c(new_model$chunks_start, new_chunk_start + object$offset - new_model$offset)
		new_model$chunks_end <- c(new_model$chunks_end, new_chunk_end + object$offset - new_model$offset)
		
		new_read_id <- object$read_id[[i]][(new_chunk_start - object$chunks_start[i]+1):(new_chunk_end + 1)]
		read_ids <- unique(unlist(new_read_id))
		new_model$alleles[[length(new_model$alleles)+1]] <- object$alleles[[i]][read_ids]
		new_model$strands[[length(new_model$strands)+1]] <- object$strands[[i]][read_ids]
		new_model$unique[[length(new_model$unique)+1]] <- object$unique[[i]][read_ids]
		
		if("reads" %in% names(object)) {
			
			tmp <- list()
			
			tmp$positions <- object$reads[[i]]$positions[read_ids]-start_pos
			tmp$lengths <- object$reads[[i]]$lengths[read_ids]
			tmp$reads <- object$reads[[i]]$reads[read_ids]
			tmp$quality <- object$reads[[i]]$quality[read_ids]
			
			new_model$reads[[length(new_model$reads)+1]] <- tmp
		}		
		
		#Recalibrate read_ids
		new_id <- 1:length(read_ids)
		names(new_id) <-sort(read_ids)
		new_read_id <- lapply(new_read_id, function(x) as.integer(new_id[as.character(x)]))
		new_model$read_id[[length(new_model$read_id)+1]] <- new_read_id 
		
	}	
		
	new_model$number_of_chunks <- length(new_model$chunks_start)
	
	#FIXME class
	class(new_model) <- "epiG"
	
	return(new_model)
}

#' start
#' 
#' @param object 
#' @returnType 
#' @return numeric
#' 
#' @author martin
#' @export
start.epiG <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(object$offset)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, start))
	}
	
	stop("Unknown class")
}

#' end
#' 
#' @param object 
#' @returnType numeric
#' @return numeric
#' 
#' @author martin
#' @export
end.epiG <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(object$offset + object$length - 1)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, end))
	}
	
	stop("Unknown class")
}

#' toRel
#' 
#' @param position 
#' @param object 
#' @returnType 
#' @return vector
#' 
#' @author martin
#' @export
toRel <- function(position, object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(position - object$offset)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(position - min(start(object)))
	}
	
	stop("Unknown class")
	
}

#' genotype
#' 
#' @param x 
#' @returnType 
#' @return genotype
#' 
#' @author martin
#' @export
genotype <- function(x) {
	
	if(paste(class(x), collapse = ".") == "integer") {
		return(c(0,1,2,3,4,1,2)[x+1])
	}
	
	if(paste(class(x), collapse = ".") == "matrix") {
		return(apply(x, c(1,2), genotype))
	}
	
	if(paste(class(x), collapse = ".") == "epiG") {
		return(genotype(x$g))
	}
	
        if(paste(class(x), collapse = ".") == "epiG.chunks") {
                s <- min(start(x))
                e <- max(end(x))
                n <- x[[1]]$config$max_alleles

                r <- matrix(nrow = n, ncol = e-s+1)

                for(model in x) {
                        r[,start(model):end(model)-s+1] <- genotype(model)
                }

                return(r)
        }

	stop("Unknown class")
}

#' methylation
#' 
#' @param x 
#' @returnType 
#' @return methylation
#' 
#' @author martin
#' @export
methylation <- function(x) {
	
	if(paste(class(x), collapse = ".") == "integer") {
		return(c(NA,FALSE,FALSE,NA,NA,TRUE,TRUE)[x+1])
		#return(c(0,1,2,0,0,5,6)[x+1])
	}
	
	if(paste(class(x), collapse = ".") == "matrix") {
		return(apply(x, c(1,2), methylation))
	}
	
	if(paste(class(x), collapse = ".") == "epiG") {
		return(methylation(x$g))
	}
	
	if(paste(class(x), collapse = ".") == "epiG.chunks") {
		s <- min(start(x))
		e <- max(end(x))
		n <- x[[1]]$config$max_alleles
		
		r <- matrix(nrow = n, ncol = e-s+1)
		
		for(model in x) {
			r[,start(model):end(model)-s+1] <- methylation(model)
		}
		
		return(r)
	}
	
	stop("Unknown class")
}

#' occupancy
#' 
#' @param x 
#' @returnType 
#' @return occupancy
#' 
#' @author martin
#' @export
occupancy <- function(x) {
	
	if(paste(class(x), collapse = ".") == "epiG") {
		return(x$occupancy)
	}
	
	if(paste(class(x), collapse = ".") == "epiG.chunks") {
		s <- min(start(x))
		e <- max(end(x))
		n <- x[[1]]$config$max_alleles+1
		
		r <- matrix(nrow = n, ncol = e-s+1)
		
		for(model in x) {
			r[,start(model):end(model)-s+1] <- occupancy(model)
		}
		
		return(r)
	}
	
	stop("Unknown class")
}


#' locate_variation
#' 
#' @param x 
#' @returnType 
#' @return locate_variation
#' 
#' @author martin
#' @export
locate_variation <- function(x) {
	
	if(paste(class(x), collapse = ".") == "matrix") {
		return(which(apply(x, 2, function(x) sum(c(1,2,3,4,5,6) %in% x) > 1 )))
	}
	
	if(paste(class(x), collapse = ".") == "epiG") {
		return(locate_variation(x$g))
	}
	
	#FIXME class list
		
	stop("Unknown class")
} 

#' coverage
#' 
#' @param object 
#' @returnType 
#' @return coverage
#' 
#' @author martin
#' @export
coverage <- function(object) {
	lapply(object$read_id, function(x) sapply(x, function(y) length(y)))
	
	#FIXME class list
	
}

#' HCG
#' 
#' @param seq 
#' @returnType 
#' @return TRUE FALSE
#' 
#' @author martin
#' @export
HCG <- function(seq) {
	return(seq[1] %in% c(1,5,3,4) & seq[2] %in% c(1,5) & seq[3] %in% c(2,6))	
}

#' GCH
#' 
#' @param seq 
#' @returnType 
#' @return TRUE FALSE
#' 
#' @author martin
#' @export
GCH <- function(seq) {
	return(seq[3] %in% c(1,5,3,4) & seq[2] %in% c(1,5) & seq[1] %in% c(2,6))	
}

#FIXME more flexiable locater and fast c++ impl using openmp
#' locate
#' 
#' @param object 
#' @param .what 
#' @returnType 
#' @return positions
#' 
#' @author martin
#' @export
locate <- function(object, .what) {
	
	if(paste(class(object), collapse = ".") == "epiG") {

                if(is.null(object[["ref"]])) {
                    stop("No ref genom found")
                }

                return(which(sapply(1:length(object$ref), function(i) .what(object$ref[i:(i+4)])))+object$offset)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
                return(unlist(lapply(object, function(x) locate(x, .what))))
	}
	
	stop("Unknown class")
}
