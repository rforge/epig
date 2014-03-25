/*
 * types.hpp
 *
 *  Created on: Nov 3, 2013
 *      Author: martin
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <armadillo>
using namespace arma;

typedef uword t_index;
typedef Col<t_index> t_indices;

typedef u32 t_count;
typedef Col<t_count> t_counts;

typedef t_count t_length;
typedef t_counts t_lengths;

typedef int t_position;
typedef Col<int> t_positions;

typedef char t_base; //coding 100 -N, 0 - C, 1 - G, 2 - A, 3 - T
typedef Col<t_base> t_seq_bases;

typedef Col<double> t_epsilon_quality;

typedef char t_epi_base; //coding 0 - Ø (No allele), 1 - C, 2 - G, 3 - A, 4 - T, 5 - C^me, 6 - G_me
typedef Col<t_epi_base> t_genotype; // vector of length (max number of alleles)
typedef Mat<t_epi_base> t_methylome; //(max number of alleles) x (length of sequence) //TODO try with sparse matrix

/*
  	Subset coding

 	0   = No coverage
	1 C = cytosine
	2 G = guanine
 	3 A = adenine
	4 T = thymine
	5 R = G A (purine)
	6 Y = T C (pyrimidine)
	7 K = G T (keto)
	8 M = A C (amino)
	9 S = G C (strong bonds)
	10 W = A T (weak bonds)
	11 B = G T C (all but A)
	12 D = G A T (all but C)
	13 H = A C T (all but G)
	14 V = G C A (all but T)
	15 N = A G C T (any)
*/
typedef char t_base_subset;
typedef Col<t_base_subset> t_base_subsets;

typedef char t_strand; // Coding : 0 - forward 1 - reverse
typedef Col<t_strand> t_strands; // vector of length (number of reads)

typedef unsigned int t_allele; //Extended allele indicator, element in 0, ..., max_alleles, max_alleles -> read unmatched
typedef Col<t_allele> t_alleles;

typedef Col<double> t_loglike_vector;

typedef Col<int> t_local_allele_occupancy;
typedef Mat<int> t_allele_occupancy; //matrix of size (max_allele+1) x (length of sequence)

typedef Mat<double> t_model; //{C, G, A, T} x {C, G, A, T, C^me, G_me}
typedef field<t_model> t_models;

typedef Col<double>::fixed<15> t_prior_vector; //TODO

#endif /* TYPES_HPP_ */
