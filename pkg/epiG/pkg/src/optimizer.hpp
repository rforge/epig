/*
 * optimizer.hpp
 *
 *  Created on: Jan 26, 2014
 *      Author: martin
 */

#ifndef OPTIMIZER_HPP_
#define OPTIMIZER_HPP_

class optimizer {

	//Algorithm configurations
	t_count const max_iterations;

	//Allele configuartion
	t_count const max_alleles;
	double const delta; //match threshold, delta in [0,1]

	//Model parameters
	t_methylome methylome;
	t_strands strands;
	t_alleles alleles;
	t_strands unique_assigned; //TODO type


	//Helper classes
	mutable genotype_iterator genotype_itr;
	reference_genome_prior genotype_prior;

	//Allele occ weights
	arma::Col<double> allele_occupancy_penalty_weights;

	//Data
	arma::field<arma::Col<double> > loglike_sums; // field of size seq length
	t_allele_occupancy allele_occupancy; //mat of size max_alleles x seq length

	alignment_data const& data;

	template<typename abort_checker>
	t_count optimize_profile(abort_checker const& ac);

	template<typename abort_checker>
	void optimize_genotype(abort_checker const& ac);

	t_loglike_vector compute_profile_loglike(t_index const read_number) const;

	double local_sum(t_position const pos, t_index const genotype_index) const;
	t_local_allele_occupancy local_allele_occupancy(t_position const pos) const;
	void update(t_index const read_number, t_strand const old_strand, t_allele const old_allele);

	double compute_loglike_term(t_index const read_number, t_position const pos, t_genotype const& genotype, t_allele const allele, t_strand const strand) const;
	double compute_genotype_prior(t_position const pos, t_genotype const& genotype, t_local_allele_occupancy const& occupancy) const;

	double compute_allele_occupancy_penalty(t_local_allele_occupancy const& occupancy) const;
	t_local_allele_occupancy compute_updated_allele_occupancy(t_position const pos, t_index const read_number, t_allele const old_allele, t_allele const new_allele) const;

	bool is_feasible(t_genotype const& genotype, t_local_allele_occupancy const& local_allele_occupancy) const;

	t_position read_start_postion(t_index const read_number) const;
	t_position read_end_postion(t_index const read_number) const;

public:

	optimizer(AlgorithmConfiguration const& config, alignment_data const& data, std::string const& string_ref, std::string const& string_alt);

	template<typename abort_checker>
	void run(abort_checker const& ac);

	double compute_objective() const;

	t_methylome get_methylome() const {
		return methylome;
	}

	t_alleles get_alleles() const {
		return alleles;
	}

	t_strands get_strands() const {
		return strands;
	}

	t_strands get_unique_assigned() const {
		return unique_assigned;
	}

	t_allele_occupancy get_allele_occupancy() const {
		return allele_occupancy;
	}
};

inline optimizer::optimizer(AlgorithmConfiguration const& config, alignment_data const& data, std::string const& string_ref, std::string const& string_alt) :
		max_iterations(config.max_iterations), max_alleles(config.max_alleles), delta(config.allele_match_prob_threshold), methylome(max_alleles, data.sequence_length),
		strands(data.n_reads), alleles(data.n_reads), unique_assigned(data.n_reads), genotype_itr(config, data), genotype_prior(string_ref, string_alt, config.log_prior(0), config.log_prior(1)),
		allele_occupancy_penalty_weights(max_alleles+1), loglike_sums(data.sequence_length), allele_occupancy(max_alleles+1, data.sequence_length), data(data) {


	//Initialize:

	//init strands
	strands.zeros();

	//TODO
	unique_assigned.zeros();

	//init alleles
	alleles.fill(max_alleles); //set all reads to unmatched

	//init allele_occupancy
	allele_occupancy.zeros();
	for(t_index i = 0; i < data.n_reads; ++i) {
		allele_occupancy.submat(alleles(i), read_start_postion(i), alleles(i), read_end_postion(i)) = allele_occupancy.submat(alleles(i), read_start_postion(i), alleles(i), read_end_postion(i)) + 1;
	}

	//init allele_occupancy_penalty_weights
	for(u32 i = 0; i < max_alleles+1; ++i) {
		allele_occupancy_penalty_weights(i) = i * config.allele_occupancy_weight;
	}

	//init loglike_sums
	for(t_position j = 0; j < data.sequence_length; ++j) {
		loglike_sums(j).set_size(genotype_itr.get_iterator_size(j));
		loglike_sums(j).zeros();
	}

	for(t_index i = 0; i < data.n_reads; ++i) {
		for(t_position j = read_start_postion(i); j <= read_end_postion(i); ++j) {
			for (genotype_itr.set_position(j); genotype_itr.has_next();
						++genotype_itr) {

					t_genotype genotype = *genotype_itr;

					loglike_sums(j)(genotype_itr.get_index()) += compute_loglike_term(
							i, j, genotype, alleles(i), strands(i));
			}
		}
	}
}

template<typename abort_checker>
inline void optimizer::run(const abort_checker& ac) {

    TIMER_START

	t_count itr = 1;

	while (itr > 0) {
		itr = optimize_profile(ac);
	}

	optimize_genotype(ac);

}

template<typename abort_checker>
inline t_count optimizer::optimize_profile(const abort_checker& ac) {

	t_count i = 0;
	for (; i < max_iterations; i++) {

        if (ac.is_aborted()) {
			return 0;
		}

		t_count changes = 0;

		t_alleles::iterator allele_itr = alleles.begin();
		t_strands::iterator strand_itr = strands.begin();

		for (t_index read_number = 0; read_number < data.n_reads;
				++read_number, ++allele_itr, ++strand_itr) {

			//Check for ctrl C
			if (ac.check_abort()) {
				break;
			}

			// index = 2*allele + strand, allele = 0, .., max_allele - 1, stand = 0, 1
			t_loglike_vector loglike = compute_profile_loglike(read_number);

			t_index i = argmax(loglike);

			// Uniquely assigned
			arma::uvec indices = arma::abs((loglike - loglike(i))/loglike(i)) < 1e-5;

			bool unique_allele = 0;
			bool unique_strand = 0;

			for(t_index j = 0; j < indices.n_elem; ++j) {
				//unique allele
				unique_allele = i/2 == indices(j)/2;

				//unique strand
				unique_strand = i % 2 == indices(j) % 2;
			}

			unique_assigned(read_number) = 2*unique_allele + unique_strand;

			//TODO this is a DEBUG CHECK
			if( ! is_finite(loglike(i))) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index old_i = *allele_itr != max_alleles ? (*allele_itr) * 2 + *strand_itr : 2*max_alleles;

			if (is_finite(loglike(old_i)) && abs(loglike(i) - loglike(old_i)) < 1e-5) {
				continue;
			}

			t_allele old_allele = *allele_itr;
			t_strand old_strand = *strand_itr;

			if(i < 2*max_alleles) {
				*strand_itr = i % 2;
				*allele_itr = i/2;
			}

			else {
				*allele_itr = max_alleles;
			}

			//TODO remove
			//cout << old_allele << " : " << *allele_itr << endl;
			//cout << loglike(old_i) << " : " << loglike(i) << endl;

			update(read_number, old_strand, old_allele);

			//add change
			changes++;

			//cout << read_number << " : " << compute_objective() << endl;

		}

		if (changes == 0) {
			break;
		}
	}

	if (i == max_iterations) {
		report_error("Max iteration limit reached");
	}

	return i;

}

inline t_loglike_vector optimizer::compute_profile_loglike(
		t_index const read_number) const {

    TIMER_START;

	t_position const read_start = read_start_postion(read_number);
	t_position const read_end = read_end_postion(read_number);

	t_allele const old_allele = alleles(read_number);
	t_strand const old_strand = strands(read_number);

    t_loglike_vector loglike(max_alleles * 2 + 1, arma::fill::zeros);
	t_loglike_vector tmp_max(max_alleles * 2 + 1);
	field<t_local_allele_occupancy> occupancy(max_alleles+1);
	Col<double> occupancy_penalty(max_alleles+1);

	for (t_position pos = read_start; pos <= read_end; pos++) {

		if(data.base_subsets(pos) == 0) {
			//No covarge -- all observed bases overlapping position are N
			continue;
		}

		tmp_max.fill(-std::numeric_limits<double>::infinity());

		//Compute new allele occupancies
		unsigned int n_occupancy = accu(allele_occupancy.submat(0, pos, max_alleles-1, pos) != 0);

		for(t_allele a = 0; a <= max_alleles; ++a) {
			occupancy(a) = compute_updated_allele_occupancy(pos, read_number, old_allele, a);
			occupancy_penalty(a) = compute_allele_occupancy_penalty(occupancy(a));
		}

		/////////////////////////
		//
		// Compute genotype minimizer
		//

		for (genotype_itr.set_position(pos); genotype_itr.has_next(); ++genotype_itr) {

			//Chek if genotype is possible
			if(genotype_itr.get_occ() > n_occupancy + 1 || genotype_itr.get_occ() + 1 < n_occupancy) {
				continue;
			}

			t_genotype genotype = *genotype_itr;

			double local_loglike = local_sum(pos, genotype_itr.get_index());

			//Remove term coming from current read
			local_loglike -= compute_loglike_term(read_number, pos, genotype, old_allele, old_strand);

			//TODO choose type of genotype prior

			for (t_allele a = 0; a < max_alleles; ++a) {

				//Check feasibility
				if(!is_feasible(genotype, occupancy(a))) {
					continue;
				}

				double prior = compute_genotype_prior(pos, genotype, occupancy(a)) + occupancy_penalty(a);

				//Allele = a, Strand = 0
				double s0 = local_loglike
						+ compute_loglike_term(read_number, pos, genotype, a, 0)  + prior;

				//Allele = a, Strand = 1
				double s1 = local_loglike
						+ compute_loglike_term(read_number, pos, genotype, a, 1) + prior;

				if (s0 > tmp_max(2 * a)) {
					tmp_max(2 * a) = s0;
				}

				if (s1 > tmp_max(2 * a + 1)) {
					tmp_max(2 * a + 1) = s1;
				}

			}

			//Unmatched reads
			t_allele a = max_alleles;

			//Check feasibility
			if(!is_feasible(genotype, occupancy(a))) {
				continue;
			}

			double prior = compute_genotype_prior(pos, genotype, occupancy(a)) + occupancy_penalty(a);

			double s = local_loglike + compute_loglike_term(read_number, pos, genotype, a, old_strand) + prior;

			if (s > tmp_max(2 * max_alleles)) {
				tmp_max(2 * max_alleles) = s;
			}

		}

		//summing over positions
		loglike += tmp_max;

	}

	return loglike;

}

inline double optimizer::compute_loglike_term(t_index const read_number, t_position const pos,
		t_genotype const& genotype, t_allele const allele, t_strand const strand) const {

	if(allele == max_alleles) {
		//unmatched read
		return data.log_pmax(read_number) + log(delta);
	}

	if(genotype(allele) == 0) {
		//not feasible
		return 0;
	}

	//Matched read
	return data.loglike_terms(read_number, strand)(pos - read_start_postion(read_number), genotype(allele) - 1);
}

inline double optimizer::compute_genotype_prior(t_position const pos,
		t_genotype const& genotype, t_local_allele_occupancy const& occupancy) const {
	return genotype_prior.getPrior(pos, genotype, occupancy);
}

inline bool optimizer::is_feasible(t_genotype const& genotype,
		t_local_allele_occupancy const& local_allele_occupancy) const {

	for (t_allele a = 0; a < max_alleles; ++a) {
		if((genotype(a) == 0 && local_allele_occupancy(a) != 0) || (genotype(a) != 0 && local_allele_occupancy(a) == 0) ) {
			return false;
		}
	}

	return true;
}

inline double optimizer::local_sum(t_position const pos, t_index const genotype_index) const {
	return loglike_sums(pos)(genotype_index);
}

inline t_local_allele_occupancy optimizer::local_allele_occupancy(
		t_position const pos) const {
	return allele_occupancy.col(pos);
}

inline void optimizer::update(t_index const read_number, t_strand const old_strand,
		t_allele const old_allele) {

	//cout << read_number << " : " << old_strand << " : " << old_allele << endl;

	//Check if sum needs to be updated
	if (old_strand == strands(read_number)
			&& old_allele == alleles(read_number)) {
		return;
	}

	t_position start = read_start_postion(read_number);
	t_position end = read_end_postion(read_number);

	for (t_position pos = start; pos <= end; ++pos) {

		//// allele occupancy
		allele_occupancy(old_allele, pos) = allele_occupancy(old_allele, pos) - 1;
		allele_occupancy(alleles(read_number), pos) = allele_occupancy(alleles(read_number), pos) + 1;

		//TODO debug check
		if(accu(allele_occupancy < 0) != 0) {
			throw std::runtime_error("update - internal error");
		}

		for (genotype_itr.set_position(pos); genotype_itr.has_next();
				++genotype_itr) {

			t_genotype genotype = *genotype_itr;

			//// loglike sums
			//Remove old term for read_updated from the sum
			loglike_sums(pos)(genotype_itr.get_index()) -= compute_loglike_term(
					read_number, pos, genotype, old_allele, old_strand);
			//Add current
			loglike_sums(pos)(genotype_itr.get_index()) += compute_loglike_term(
					read_number, pos, genotype, alleles(read_number),
					strands(read_number));
		}
	}
}

template<typename abort_checker>
inline void optimizer::optimize_genotype(const abort_checker& ac) {

	methylome.zeros();

	for (t_position pos = 0; pos < data.sequence_length; ++pos) {

			//Check for ctrl C
			if (ac.check_abort()) {
				return;
			}

			if(data.base_subsets(pos) == 0) {
				//No covarge -- all observed bases overlapping position are N
				continue;
			}

			if (sum(local_allele_occupancy(pos).subvec(0, max_alleles-1) != 0) == 0) {
				continue;
			}

			t_genotype new_g(max_alleles, arma::fill::zeros);

            //TODO remove
            //double old_max = -std::numeric_limits<double>::infinity();
			double tmp_max = -std::numeric_limits<double>::infinity();

			for (genotype_itr.set_position(pos); genotype_itr.has_next(); ++genotype_itr) {

				t_genotype genotype = *genotype_itr;

				if( ! is_feasible(genotype, local_allele_occupancy(pos))) {
					continue;
				}

				double local_loglike = local_sum(pos, genotype_itr.get_index());

				//TODO Genotype prior
				//local_loglike += compute_genotype_prior(pos, genotype);
				local_loglike += compute_genotype_prior(pos, genotype, local_allele_occupancy(pos));

				if (local_loglike > tmp_max) {
					tmp_max = local_loglike;
					new_g = genotype;
				}

			}

			//TODO remove
			//cout << trans(conv_to<uvec>::from(new_g)) << endl;

			//TODO debug check
			if(is_zero(new_g)) {
				throw std::runtime_error("optimize_genotype - internal error");
			}

			methylome.col(pos) = new_g;
		}
}

inline t_position optimizer::read_start_postion(t_index const read_number) const {
	return data.reads_start_postions(read_number);
}

inline t_position optimizer::read_end_postion(t_index const read_number) const {
	return data.reads_end_postions(read_number);
}

inline t_local_allele_occupancy optimizer::compute_updated_allele_occupancy(
		t_position const pos, t_index const read_number, t_allele const old_allele, t_allele const new_allele) const {

	t_local_allele_occupancy occupancy(allele_occupancy.col(pos));

	occupancy(old_allele) = occupancy(old_allele) - 1;
	occupancy(new_allele) = occupancy(new_allele) + 1;

	return occupancy;
}

inline double optimizer::compute_objective() const {

	double objective = 0;

	for (t_position pos = 0; pos < data.sequence_length; ++pos) {

        //TODO remove
        //double old_max = -std::numeric_limits<double>::infinity();
		double tmp_max = -std::numeric_limits<double>::infinity();

		for (genotype_itr.set_position(pos); genotype_itr.has_next();
				++genotype_itr) {

			t_genotype genotype = *genotype_itr;
			t_local_allele_occupancy occupancy = local_allele_occupancy(pos);

			if (!is_feasible(genotype, occupancy)) {
				continue;
			}

			double local_loglike = local_sum(pos, genotype_itr.get_index());

			//Genotype prior
			//local_loglike += compute_genotype_prior(pos, genotype);
			local_loglike += compute_genotype_prior(pos, genotype, local_allele_occupancy(pos));

			local_loglike += compute_allele_occupancy_penalty(allele_occupancy.col(pos));


			if (local_loglike > tmp_max) {
				tmp_max = local_loglike;
			}

		}

		objective += tmp_max;

	}

	return objective;
}

inline double optimizer::compute_allele_occupancy_penalty(
	t_local_allele_occupancy const& occupancy) const {

	double c = sum(occupancy); //FIXME pre-compute

	return - 1/c*dot(allele_occupancy_penalty_weights.subvec(0,max_alleles-1), sort(occupancy.subvec(0,max_alleles-1), 1))
					- 1/c*allele_occupancy_penalty_weights(sum(occupancy.subvec(0,max_alleles-1) != 0))*occupancy(max_alleles);


}

#endif /* OPTIMIZER_HPP_ */
