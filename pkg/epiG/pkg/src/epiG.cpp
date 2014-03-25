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
#define PRINT_BACKTRACE

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

#ifdef PRINT_BACKTRACE
#include "Backtrace.h"
#endif

#include <string>
#include <armadillo>
#include <Rcpp.h>
#include "arma_additions.h"
#include "rtools.h"
#include "simple_timer.h"

#include "math_tools.hpp"
#include "types.hpp"
#include "bamReader.hpp"
#include "epiG_algorithm_config.h"
#include "reference_genome_prior.hpp"
#include "alignment_data.hpp"
#include "genotype_iterator.hpp"
#include "optimizer.hpp"
#include "chunk_optimizer.hpp"

//R interface
#include "epiG_R_interface.hpp"
