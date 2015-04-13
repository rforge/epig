/*
 * msbs.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: martin
 */

//FIXME auto detect openmp
//FIXME openmp not available warning

//Progress monitor
#include <progress.hpp>

#include <RcppCommon.h>
#include <Rconfig.h>
#include <RcppArmadilloConfig.h>

//#define DO_TIMING
#define EPIG_DEBUG
#define EPIG_USE_OPENMP

// Debugging
#ifdef EPIG_DEBUG
// Do debugging
#ifdef ARMA_NO_DEBUG
#undef ARMA_NO_DEBUG
#endif
#ifdef NDEBUG
#undef NDEBUG
#endif
#else
// Do no debugging
#define ARMA_NO_DEBUG
#define NDEBUG
#endif

#include <string>
#include <armadillo>
#include <Rcpp.h>
#include "arma_additions.h"
#include "rtools.h"
#include "simple_timer.h"

#include "math_tools.h"
#include "types.h"
#include "seq_tools.h"
#include "bam_reader.h"
#include "epiG_algorithm_config.h"
#include "alignment_data.h"
#include "prior.h"
#include "chain_opt.h"
#include "genotype_opt.h"
#include "chunked_chain_opt.h"

//// R interface
#include "rinterface/epiG.h"
#include "rinterface/fetch.h"

// vector search
#include "vector_search.h"
