# TODO: Add comment
# 
# Author: martin
###############################################################################

#' start
#' 
#' @param object 
#' @return numeric
#' 
#' @author Martin Vincent
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
#' @return numeric
#' 
#' @author Martin Vincent
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

#' Length of model in base pairs
#' @param object 
#' @return Length of model in base pairs
#' 
#' @author Martin Vincent
#' @export
length.epiG <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(object$length)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sum(sapply(object, length)))
	}
	
	stop("Unknown class")
}

#' Number of reads in model
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
nread <- function(object, ... ) UseMethod("nread")

#' Number of reads in model
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
nread.epiG <- function(object, ...)  {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(length(unique(unlist(object$read_ids))))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, nread))
	}
	
	stop("Unknown class")
	
}

#' genotype
#' @param object 
#' @param pos 
#' @param remove.meth 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
genotype <- function(object, pos, remove.meth, ... ) UseMethod("genotype")

#' genotype
#' codeing C = 1, G = 2, A = 3, T = 4
#' @param object 
#' @param pos 
#' @param remove.meth 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
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

#' methylation
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
methylation <- function(object, pos, ... ) UseMethod("methylation")

#' methylation
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
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

#' strand
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
strand <- function(object, pos, ... ) UseMethod("strand")

#' strand
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
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

#' coverage
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
coverage <- function(object, pos = NULL, ... ) UseMethod("coverage")

#' coverage
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
coverage.epiG <- function(object, pos = NULL, ...) {
	
	if(is.null(pos)) {
		return(sapply(start(object):end(object), function(pos) coverage(object, pos)))
	}
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(start(object) > pos || end(object) < pos) {
			stop("Position not in range")
		}
		
		if(length(object$read_ids) == 0) {
			return(0)
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

#' position.info
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
position.info <- function(object, pos, ... ) UseMethod("position.info")

#' position.info
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
position.info.epiG <- function(object, pos, ...) {
	
	if(length(pos) == 0) {
		return(NULL)
	}
	
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
			return(data.frame(position = pos, chain.id = NA, ref = NA, alt = NA, genotype = NA, methylated = NA, strand = NA, coverage = 0))
		}	
		
		chains <- sort(object$haplotype$chain[object$read_ids[[pos - start(object)+1]]])

		info.df <- data.frame(position = pos, chain.id = unique(chains), ref = NA, alt = NA, genotype = symbols(genotype(object, pos, remove.met = TRUE)), methylated = methylation(object, pos), strand = strand(object, pos), coverage = sapply(unique(chains), function(x) sum(chains == x)))
				
		if(!is.null(object[["ref"]])) {
			info.df$ref <- symbols(object$ref[pos - object$offset + 1])
		} 
		
		if(!is.null(object[["alt"]])) {
			info.df$alt <- symbols(object$alt[pos - object$offset + 1])
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

#' chain.info
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
chain.info <- function(object, ... ) UseMethod("chain.info")

#' chain.info
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
chain.info.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(data.frame(chain.id = sort(unique(fit$haplotype$chain)), start = object$haplotype$start, end = object$haplotype$end, length = object$haplotype$end - object$haplotype$start + 1, nreads = as.vector(table(object$haplotype$chain)), strand = object$strands))
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		tmp <- lapply(object, function(x) chain.info(x, ...))		
		return(do.call("rbind", tmp))
	}
	
	stop("Unknown class")
	
}

#' read.info
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
read.info <- function(object, ... ) UseMethod("read.info")

#' read.info
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
read.info.epiG <- function(object, ...) {
	
	if(!("reads" %in% names(object))) {
		object <- fetch.reads(object)
	}
		
	if(paste(class(object), collapse = ".") == "epiG") {

		info <- NULL
		
		for(idx in 1:nread(object)) {
			tmp <- data.frame(position = object$reads$position[idx]:(object$reads$position[idx]+object$reads$length[idx]-1) + object$offset, 
					symbol = symbols(object$reads$reads[[idx]]), read.id = idx, quality = object$reads$quality[[idx]], chain.id=object$haplotype$chain[idx], strand=object$strands[object$haplotype$chain[idx]])
			
			info <- rbind(info, tmp)	
		}
		
		return(info)
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		tmp <- lapply(object, function(x) read.info(x, ...))
		
		read.id.offset <- 0
		for(i in 1:length(tmp)) {
			tmp[[i]]$read.id <- tmp[[i]]$read.id + read.id.offset
			read.id.offset <- max(tmp[[i]]$read.id)
		}
		
		return(do.call("rbind", tmp))
	}
	
	stop("Unknown class")
	
}

#' Number of chunks 
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
nchunks <- function(object, pos, ... ) UseMethod("nchunks")

#' Number of chunks 
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
nchunks.epiG <- function(object, ...) {
	if(paste(class(object), collapse = ".") == "epiG") {
		return(1)
	}
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sum(sapply(object, nchunks)))
	}
		
}

#' Number of chains
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
nchain <- function(object, pos, ... ) UseMethod("nchain")

#' Number of chains
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
nchain.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(length(object$haplotype$start))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, nchains))
	}
	
}

#' subregion
#' @param object 
#' @param start 
#' @param end 
#' @param chop.reads 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
subregion <- function(object, start, end, chop.reads = TRUE, ... ) UseMethod("subregion")

#' subregion
#' @param object 
#' @param start 
#' @param end 
#' @param chop.reads 
#' @param ... 
#' @return ??
#' 
#' @author martin
#' @export
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
	