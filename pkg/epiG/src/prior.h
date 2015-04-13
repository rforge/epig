/*
 * reference_genome_prior.hpp
 *
 *  Created on: Jan 14, 2014
 *      Author: martin
 */

#ifndef REFERENCE_GENOME_PRIOR_HPP_
#define REFERENCE_GENOME_PRIOR_HPP_

//TODO simplify
//TODO priors when no ref available - prior(pos, genotype)


class reference_genome_prior {

t_seq_bases ref;
t_seq_bases alt;
const t_prior_vector priors;
const t_prior_vector priors_alt;

arma::Mat<double> log_binom_table;

char refmap(char ref, char alt, char base) const;

public:

reference_genome_prior(t_seq_bases const& ref, t_seq_bases const& alt, AlgorithmConfiguration const& config);
double getPrior(t_position pos, t_genotype const& genotype) const;
double getPrior(t_position pos, t_genotype const& genotype, t_counts const& read_counts, t_strands const& strands) const;

};

reference_genome_prior::reference_genome_prior(t_seq_bases const& ref, t_seq_bases const& alt, AlgorithmConfiguration const& config) :
		ref(ref), alt(alt), priors(config.genotype_priors(0)), priors_alt(config.genotype_priors(1)),
		log_binom_table(make_log_binom_half_table(config.reads_hard_limit)) {

	//TODO remove
//	cout << log_binom_table(0,0) << endl;
//	cout << log_binom_table(1,0) << endl;
//	cout << log_binom_table(10,1) << endl;
//	cout << log_binom_table(10,5) << endl;
//	cout << log_binom_table(11,5) << endl;

}

// Return new coding: Ref = 1, Alt = 2, others 3 and 4.
// If alt = N then: Ref = 1, others 2, 3, 4
// base coding 0 - C, 1 - G, 2 - A, 3 - T, 4 - C^me, 5 - G_me
inline char reference_genome_prior::refmap(char ref, char alt, char base) const {

	char g = base % 4 + 1; //remove met,

	if(ref == 0) {
		return  1; //FIXME thing this through + warning msg
		//throw std::runtime_error("No ref found"); //TODO update error msg
	}

	if(alt == 0) {

		if(g == ref) {
			return 1;
		}

		if(g == 1) {
			return ref;
		}

		return g;
	}

	// alt != 0

	if(g == ref) {
		return 1;
	}

	if(g == alt) {
		return 2;
	}

	if(g == 1) {
		return ref;
	}

	if(g == 2) {
		return alt;
	}

	return g;
}

inline double reference_genome_prior::getPrior(t_position pos, t_genotype const& genotype, t_counts const& read_counts, t_strands const& strands) const {

	unsigned int n_rev = sum(strands);
	if(n_rev == 0 || n_rev == strands.n_elem) {
		return getPrior(pos, genotype);
	}

	t_counts fwd_counts(4, arma::fill::zeros);
	t_counts rev_counts(4, arma::fill::zeros);

	// Collect statistics
	t_genotype::const_iterator g = genotype.begin();
	t_strands::const_iterator s = strands.begin();
	t_counts::const_iterator rc = read_counts.begin();

	for(;g != genotype.end(); ++g, ++s, ++rc) {

		switch (*s) {
			case strand_fwd:
				fwd_counts(*g % 4) += *rc;
				break;
			case strand_rev:
				rev_counts(*g % 4) += *rc;
				break;
			default:
				throw std::runtime_error("getPrior : Error");
				break;

		}
	}

	// Compute fwd rev balance prior
	t_counts::const_iterator fwd_c = fwd_counts.begin();
	t_counts::const_iterator rev_c = rev_counts.begin();

	double balance_prior = 0;
	for(;fwd_c != fwd_counts.end(); ++fwd_c, ++rev_c) {
		balance_prior += log_binom_table(*fwd_c + *rev_c, min(*fwd_c, *rev_c));
	}

	//Return total prior balance prior + genotype prior

	return balance_prior + getPrior(pos, genotype);
}


//TODO pre-compute prior based on ref and genotype
inline double reference_genome_prior::getPrior(t_position pos, t_genotype const& genotype) const {

	t_base_subset subset = 0;
	t_genotype::const_iterator g = genotype.begin();
	for(;g != genotype.end(); ++g) {
		subset = base_union(subset, refmap(ref(pos), alt(pos), *g));
	}

	//TODO remove
//	if(ref(pos) == 1) {
//	cout << static_cast<int>(ref(pos))<< " : " << static_cast<int>(alt(pos)) << endl;
//	cout << trans(conv_to<uvec>::from(genotype)) << endl;
//	cout << static_cast<int>(subset) << " : " << priors_alt(subset-1) << " : " << priors(subset-1) << endl;
//	}

	if(alt(pos) != 0) {
		return priors_alt(subset-1);
	}

	return priors(subset-1);
}


#endif /* REFERENCE_GENOME_PRIOR_HPP_ */
