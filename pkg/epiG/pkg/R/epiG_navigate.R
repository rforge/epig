# TODO: Add comment
# 
# Author: martin
###############################################################################

#' start
#' 
#' @param object 
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
		return(start(object) + length(object)-1L)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, end))
	}
	
	stop("Unknown class")
}

length.epiG <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(object$length)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sum(sapply(object, length)))
	}
	
	stop("Unknown class")
}

nread <- function(object, ... ) UseMethod("nread")

nread.epiG <- function(object, ...)  {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(length(unique(unlist(object$read_ids))))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, nread))
	}
	
	stop("Unknown class")
	
}

genotype <- function(object, pos, remove.meth, ... ) UseMethod("genotype")

# codeing C = 1, G = 2, A = 3, T = 4
genotype.epiG <- function(object, pos, remove.meth = FALSE, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
	
		if(start(object) > pos || end(object) < pos) {
			stop("Position not in range")
		}
	
		collected <- sapply(1:length(object$genotype), function(i) if((pos - object$haplotype$start)[i] >= 0 && (object$haplotype$end[i] - pos) >= 0) object$genotype[[i]][(pos - object$haplotype$start)[i]+1] else NA)
		
		if(remove.meth) {
			collected[collected %in% c(5,6)] <-  collected[collected %in% c(5,6)] %% 4
		}
		
		return(collected[!is.na(collected)])
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		stop("Not yet implemented for chunks")
	}
	
	stop("Unknown class")
	
}

methylation <- function(object, pos, ... ) UseMethod("methylation")

methylation.epiG <- function(object, pos, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		g <- genotype(object, pos)
		s <- strand(object, pos)
		
		# Remove chains where methylation is not possible
		g[!((s == "fwd" & g %in% c(1,5)) | (s == "rev" & g %in% c(2,6)))] <- NA
		
		return( g == 5 | g == 6 )
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		stop("Not yet implemented for chunks")
	}
	
	stop("Unknown class")
	
}

strand <- function(object, pos, ... ) UseMethod("strand")

strand.epiG <- function(object, pos, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(start(object) > pos || end(object) < pos) {
			stop("Position not in range")
		}
		
		collected <- sapply(1:nchain(object), function(i) if((pos - object$haplotype$start)[i] >= 0 && (object$haplotype$end[i] - pos) >= 0) as.character(object$strands[i]) else NA)
		
		return(factor(collected[!is.na(collected)], levels = c("fwd", "rev")))
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		stop("Not yet implemented for chunks")
	}
	
	stop("Unknown class")
	
}

coverage <- function(object, pos = NULL, ... ) UseMethod("coverage")

coverage.epiG <- function(object, pos = NULL, ...) {
	
	if(is.null(pos)) {
		return(sapply(start(object):end(object), function(pos) coverage(object, pos)))
	}
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(start(object) > pos || end(object) < pos) {
			stop("Position not in range")
		}
		
		return(length(object$read_ids[[pos - start(object)+1]]))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		
		stop("Not yet implemented")
		#TODO FIX
		tmp <- sapply(object, function(x) if(start(x) > pos || end(x) < pos) coverage(x) else NA)
		return(tmp[!is.na(tmp)])
	}
	
	stop("Unknown class")
	
}

position.info <- function(object, pos, ... ) UseMethod("position.info")

position.info.epiG <- function(object, pos, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(length(pos) > 1) {
			
			info.df <- NULL
			
			for(p in pos) {
				info.df <- rbind(info.df, position.info(object, p))
			}
			
			return(info.df)
		}
		
		if(coverage(object, pos) == 0) {
			#Return data.frame
			return(data.frame(position = pos, chain.id = NA, ref = NA, genotype = NA, methylated = NA, strand = NA, coverage = 0))
		}	
		
		chains <- sort(object$haplotype$chain[object$read_ids[[pos - start(object)+1]]])
		
		if(!is.null(object[["ref"]])) {
			info.df <- data.frame(position = pos, chain.id = unique(chains), ref = symbols(object$ref[pos - object$offset + 1]), genotype = symbols(genotype(object, pos, remove.met = TRUE)), methylated = methylation(object, pos), strand = strand(object, pos), coverage = sapply(unique(chains), function(x) sum(chains == x)))
		} else {
			info.df <- data.frame(position = pos, chain.id = unique(chains), ref = NA, genotype = symbols(genotype(object, pos, remove.met = TRUE)), methylated = methylation(object, pos), strand = strand(object, pos), coverage = sapply(unique(chains), function(x) sum(chains == x)))
		}
		
		return(info.df)
		
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		
		tmp <- lapply(object, function(x) position.info(x, pos[pos %in% start(x):end(x)]), ...)		
		
		# Adjust chain.id
		chain.id.offset = 0
		for(i in 1:length(tmp)) {
			tmp[[i]]$chain.id <- tmp[[i]]$chain.id + chain.id.offset
			chain.id.offset <- max(tmp[[i]]$chain.id, na.rm = TRUE)
		}
		
		return(do.call("rbind", tmp))
	}
	
	stop("Unknown class")
	
}

chain.info <- function(object, ... ) UseMethod("chain.info")

chain.info.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(data.frame(start = object$haplotype$start, end = object$haplotype$end, length = object$haplotype$end - object$haplotype$start + 1, nreads = as.vector(table(object$haplotype$chain)), strand = object$strands))
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		tmp <- lapply(object, function(x) chain.info(x, ...))		
		return(do.call("rbind", tmp))
	}
	
	stop("Unknown class")
	
}

nchunks <- function(object, pos, ... ) UseMethod("nchunks")
nchunks.epiG <- function(object, ...) {
	if(paste(class(object), collapse = ".") == "epiG") {
		return(1)
	}
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sum(sapply(object, nchunks)))
	}
		
}

nchain <- function(object, pos, ... ) UseMethod("nchain")
nchain.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(length(object$haplotype$start))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, nchains))
	}
	
}

subregion <- function(object, start, end, chop.reads = TRUE, ... ) UseMethod("subregion")
subregion.epiG <- function(object, start, end, chop.reads = TRUE, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {

		if(end < start(object) || start > end(object)) {
			stop("Out of range")
		}
				
		rel_start_pos <- max(start, start(object)) - object$offset
		rel_end_pos <- min(end, end(object)) - object$offset
		
		new_object <- NULL
		
		new_object$epiG_version <- object$epiG_version
		new_object$date <- object$date
		new_object$config <- object$config
		new_object$filename <- object$filename
		new_object$refname <- object$refname
		
		new_object$offset <- as.integer(start)
		new_object$length <- as.integer(rel_end_pos - rel_start_pos + 1)		
		
		new_object$read_ids <- object$read_ids[rel_start_pos:rel_end_pos+1]
		remaining_read_ids <- sort(unique(unlist(new_object$read_ids)))
	
		#Recalibrate read_ids
		new_id <- 1:length(remaining_read_ids)
		names(new_id) <- remaining_read_ids
		new_object$read_ids <- lapply(new_object$read_ids, function(x) as.integer(new_id[as.character(x)]))
		
		# Ref and alt
		
		if("ref" %in% names(object)) {
			new_object$ref <-  object$ref[rel_start_pos:rel_end_pos+1]
		}
		
		if("alt" %in% names(object)) {
			new_object$alt <-  object$alt[rel_start_pos:rel_end_pos+1]
		}
		
		
		# Haplotype and genotype
		new_object$haplotype$chain <- object$haplotype$chain[remaining_read_ids]
		remaining_chains <- sort(unique(new_object$haplotype$chain))
		new_object$haplotype$chain <- as.integer(factor(new_object$haplotype$chain))
					
		new_object$strands <- object$strands[remaining_chains]
		
		new_object$haplotype$start <- sapply(object$haplotype$start[remaining_chains], function(x) max(start, x))
		new_object$haplotype$end <- sapply(object$haplotype$end[remaining_chains], function(x) min(end, x))
		
		new_object$genotypes <- list()
		for(i in remaining_chains) {
				chain_rel_start_pos <- max(start, object$haplotype$start[i]) - object$haplotype$start[i]
				chain_rel_end_pos <- min(end, object$haplotype$end[i]) - object$haplotype$start[i]
				
				new_object$genotypes[[length(new_object$genotypes)+1]] <- object$genotypes[[i]][chain_rel_start_pos:chain_rel_end_pos+1]				
		}
	
		# Reads
		if("reads" %in% names(object)) {
			new_object$reads$reads <- object$reads$reads[remaining_read_ids]
			new_object$reads$quality <- object$reads$quality[remaining_read_ids]
			new_object$reads$positions <- object$reads$positions[remaining_read_ids] - rel_start_pos
			new_object$reads$lengths <- object$reads$lengths[remaining_read_ids]			

			if(chop.reads) {
			
				for(i in 1:length(new_object$reads$lengths)) {
					
					read_rel_start_pos <- max(0,-new_object$reads$positions[i]) 
					read_length <- new_object$reads$lengths[i] - read_rel_start_pos - max(0, new_object$reads$positions[i] + new_object$reads$lengths[i] - new_object$length)
												
					new_object$reads$positions[i] <- new_object$reads$positions[i] + read_rel_start_pos
					new_object$reads$lengths[i] <- read_length
					new_object$reads$reads[[i]] <- new_object$reads$reads[[i]][(read_rel_start_pos+1):(read_rel_start_pos+read_length)]
					new_object$reads$quality[[i]] <- new_object$reads$quality[[i]][(read_rel_start_pos+1):(read_rel_start_pos+read_length)]
					
				}
			}
		}
		
		class(new_object) <- "epiG"
		return(new_object)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		stop("Not yet implemented")
	}
	
}
	