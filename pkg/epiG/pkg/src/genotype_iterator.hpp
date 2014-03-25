/*
 * genotype_iterator.hpp
 *
 *  Created on: Jan 26, 2014
 *      Author: martin
 */

#ifndef GENOTYPE_ITERATOR_HPP_
#define GENOTYPE_ITERATOR_HPP_

using namespace arma;

class genotype_iterator {

	t_count const max_alleles;

	field<Col<int> > iterators; //genotype iterators, one for each subset
	t_base_subsets const base_subsets; //Bases found at position

	t_index itr_index;
	t_base_subset position_subset;

	field<t_counts> n_occ;

	bool is_genotype_achievable(t_genotype const& g, t_base_subset subset) const;
	bool is_base_achievable(t_epi_base const base, t_base_subset const subset) const;

public:

	genotype_iterator(AlgorithmConfiguration const& config, alignment_data const& data);

	void set_position(t_position pos);

	genotype_iterator & operator++ ();

	bool has_next() const;

	t_genotype operator* () const;

	t_index get_index() const {
		return itr_index;
	}

	t_index get_iterator_size(t_position const position) const {
		return iterators(base_subsets(position)).n_elem;
	}

	t_count get_occ() const {
		return n_occ(position_subset)(itr_index);
	}
};

inline bool genotype_iterator::has_next() const {
	return itr_index < iterators(position_subset).n_elem;
}

inline void genotype_iterator::set_position(t_position const pos) {
	position_subset = base_subsets(pos);
	itr_index = 0;
}

inline genotype_iterator& genotype_iterator::operator ++() {
	++itr_index;
	return *this;
}

inline t_genotype genotype_iterator::operator *() const {
	t_genotype g(max_alleles);

	toBase<char, 7>(g, iterators(position_subset)(itr_index));

	return g;
}

inline genotype_iterator::genotype_iterator(AlgorithmConfiguration const& config, alignment_data const& data) :
			max_alleles(config.max_alleles), iterators(16), base_subsets(data.base_subsets), itr_index(0), position_subset(base_subsets(0)), n_occ(16) {

	//Compute position specific iterators
	t_genotype g(max_alleles);

	for (t_base_subset subset = 0; subset <= 15; ++subset) {

		t_index j = 0;

		//Compute size
		for (t_count i = 0; i < pow<7>(max_alleles); ++i) {

			toBase<char, 7>(g, i);

			if (is_genotype_achievable(g, subset)) {
				++j;
			}
		}

		iterators(subset).set_size(j);
		n_occ(subset).set_size(j);

	//Compute indices
		j = 0;
		for (t_count i = 0; i < pow<7>(max_alleles); ++i) {

			toBase<char, 7>(g, i);

			if (is_genotype_achievable(g, subset)) {
				iterators(subset)(j) = i;
				n_occ(subset)(j) = sum(g != 0);
				++j;
			}
		}

	}

}

bool genotype_iterator::is_genotype_achievable(t_genotype const& g, t_base_subset const subset) const {

	for(u32 i = 0; i < g.n_elem; ++i) {
		if(!is_base_achievable(g(i), subset)) {
			return false;
		}
	}

	if(subset != 0 && all(g == 0)) {
		return false;
	}

	return true;
}

bool genotype_iterator::is_base_achievable(t_epi_base const base, t_base_subset const p) const {

	if (p == 0 && base != 0) {
		return false;
	}

	switch (base) {
	case 0:
		return true;

	case 1:
		if (p == 2 || p == 3 || p == 5) {
			return false;
		}

		return true;

	case 2:
		if (p == 1 || p == 4 || p == 6) {
			return false;
		}

		return true;
	case 3:
		if (p == 1 || p == 2 || p == 4 || p == 6 || p == 7 || p == 9
				|| p == 11) {
			return false;
		}

		return true;
	case 4:
		if (p == 1 || p == 2 || p == 3 || p == 5 || p == 8 || p == 9
				|| p == 14) {
			return false;
		}

		return true;
	case 5:
		if (p == 2 || p == 3 || p == 4 || p == 5 || p == 7 || p == 10
				|| p == 12) {
			return false;
		}

		return true;
	case 6:
		if (p == 1 || p == 3 || p == 4 || p == 6 || p == 8 || p == 10
				|| p == 13) {
			return false;
		}

		return true;
	default:
		throw std::runtime_error("genotype_iterator -- internal error");
		break;
	}

}

#endif /* GENOTYPE_ITERATOR_HPP_ */
