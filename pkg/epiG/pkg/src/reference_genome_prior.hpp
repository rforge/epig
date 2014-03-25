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

#include "seq_tools.hpp"

class reference_genome_prior {

t_seq_bases ref;
t_seq_bases alt;
const t_prior_vector priors;
const t_prior_vector priors_alt;

char refmap(char ref, char alt, char base) const;


char switch_code(char c1, char c2, char g) const {

	if(g == c1) {
		return c2;
	}

	if(g == c2) {
		return c1;
	}

	return g;
}

public:

reference_genome_prior(std::string const& string_ref, std::string string_alt, t_prior_vector const& priors, t_prior_vector const& priors_alt);
double getPrior(t_position pos, t_genotype const& g) const;
double getPrior(t_position pos, t_genotype const& g, t_local_allele_occupancy const& occupancy) const;
double getMinPrior() const;


};

reference_genome_prior::reference_genome_prior(std::string const& string_ref, std::string string_alt, t_prior_vector const& priors, t_prior_vector const& priors_alt) :
		ref(create_bases_vector(string_ref)), alt(create_bases_vector(string_alt)), priors(priors), priors_alt(priors_alt) {}

inline double reference_genome_prior::getMinPrior() const {
	return priors.max();
}

inline char reference_genome_prior::refmap(char ref, char alt, char base) const {

	char g = (base-1) % 4 + 1; //remove met

	if(ref == 0) {
		return g;
	}

	//--ref; //Adjust coding
	g = switch_code(1, ref, g);

	if(alt != 0) {
		//--alt; //Adjust coding
		g = switch_code(2, switch_code(1, ref, alt), g);
	}

	return g;
}

//TODO pre-compute prior based on ref and genotype
inline double reference_genome_prior::getPrior(t_position pos, t_genotype const& genotype) const {

	t_base_subset subset = 0;
	t_genotype::const_iterator g = genotype.begin();
	for(;g != genotype.end(); ++g) {

		if(*g != 0) {
			subset = base_union(subset, refmap(ref(pos), alt(pos), *g));
		}
	}

	//TODO remove
	//cout << static_cast<int>(ref(pos))<< " : " << static_cast<int>(alt(pos)) << endl;
	//cout << trans(conv_to<uvec>::from(genotype)) << " : " << priors_alt(subset-1) << " : " << priors(subset-1) << endl;

	if(alt(pos) != 0) {
		return priors_alt(subset-1);
	}

	return priors(subset-1);
}

inline double reference_genome_prior::getPrior(t_position pos,
		t_genotype const& genotype,
		t_local_allele_occupancy const& occupancy) const {

	//TODO remove
	//cout << pos << " : " << ref.n_elem << " : " << alt.n_elem << endl;

	Col<int>::fixed < 4 > supp;
	supp.zeros();

	t_base_subset subset = 0;
	t_allele a = 0;
	t_genotype::const_iterator g = genotype.begin();
	for (; g != genotype.end(); ++g, ++a) {

		if (*g != 0) {
			supp((*g - 1) % 4) += occupancy(a);
			subset = base_union(subset, refmap(ref(pos), alt(pos), *g));
		}
	}

	//TODO remove
	//cout << static_cast<int>(ref(pos))<< " : " << static_cast<int>(alt(pos)) << endl;
	//cout << trans(conv_to<uvec>::from(genotype)) << " : " << priors_alt(subset-1) << " : " << priors(subset-1) << endl;

	double supp_sum = 0;
	double supp_terms = 0;
	for (int i = 0; i < 4; ++i) {
		if (supp(i) > 0) {
			supp_sum += supp(i);
			++supp_terms;
		}
	}

	double d = 0;
	for (int i = 0; i < 4; ++i) {
		if (supp(i) > 0) {
			d += square(1/supp_terms - supp(i)/supp_sum);
		}
	}

	d = sqrt(d);

	//TODO remove
//	if(pos == 1098) {
//	//if(ref(pos) == 1 && alt(pos) == 4 && supp_sum == 15 && supp(0) == 3 && supp(3) == 12 && genotype(0) == 1) {
//	cout << " -------------------------- " << endl;
//	cout << trans(conv_to<uvec>::from(genotype)) << endl;
//	cout << trans(supp) << endl;
//	cout << pos << " : " << supp_terms << " : " << d << endl;
//	cout << priors_alt(subset - 1) + 0.5*log(1-d) << " : " << priors(subset - 1) + 0.5*log(1-d) << endl;
//	}

	if (alt(pos) != 0) {
		return priors_alt(subset - 1) + log(1-d);
	}

	return priors(subset - 1) + log(1-d);
}

#endif /* REFERENCE_GENOME_PRIOR_HPP_ */
