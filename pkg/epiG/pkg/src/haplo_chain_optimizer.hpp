#ifndef HAPLO_CHAIN_OPTIMIZER_HPP_
#define HAPLO_CHAIN_OPTIMIZER_HPP_

class haplo_chain_optimizer {

	//Algorithm configurations
	t_count const max_iterations;

	t_position const min_overlap_length;

	//Model parameters
	t_haplotype haplo; //vector of size n_reads
	t_strands strands;

	//Data
	alignment_data const& data;

	//Ref
	t_seq_bases const& ref;
	t_seq_bases const& alt;


	//Priors
	double ref_no_match_prior_log;
	double ref_match_prior_log;
	arma::Col<double> haplochain_log_prior;

	template<typename abort_checker>
	t_count optimize_profile(abort_checker const& ac);

	template<typename abort_checker>
	t_count optimize_haplochains(abort_checker const& ac);

	//TODO rename : read_profile_posterior
	t_loglike_vector profile_posterior(t_index const read_number,
			t_haplotype const& feasible_haplotypes) const;

	t_loglike_vector profile_posterior(t_indices const& reads, t_haplotype const& feasible_haplotypes) const;

	//TODO rename to: find_overlapping_haplochains
	t_haplotype compute_feasible_haplotypes(t_position const start, t_position const end) const;

	double compute_haplochain_prior(t_count const unique_chains) const;

	t_loglike_vector compute_logsum(t_indices const& reads, t_position const pos) const;

	double compute_posterior(t_haplotype const& h) const;

	double compute_delta_posterior(t_index const read_number, t_haplochain const old_chain) const;
	double compute_delta_posterior(t_indices const& reads, t_haplochain const old_chain) const;

	double compute_chain_loglike(t_haplotype const& h, t_haplochain const chain) const;
	double compute_prior(t_haplotype const& h, t_position const start, t_position const end) const;

	//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
	void update_haplotype_chain(t_haplotype & haplotype, t_haplochain const old_chain) const;
	//bool is_haplotype_chain_feasiable(t_haplochain const chain) const;

	t_strand chain_strand(t_haplochain const chain) const;

	t_indices get_overlapping_reads_in_chain(t_haplotype const& h, t_position const pos,
			t_haplochain const chain) const;

	t_position read_start_postion(t_index const read_number) const;
	t_position read_end_postion(t_index const read_number) const;

	//TODO note duplicate in genotype_optimizer
	t_position haplo_chain_start(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		return data.reads_start_postions(min(reads_in_chain));
	}

	//TODO note duplicate in genotype_optimizer
	t_count haplo_chain_end(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		t_positions end_pos = data.reads_end_postions(reads_in_chain);
		return max(end_pos);
	}

public:

	haplo_chain_optimizer(AlgorithmConfiguration const& config,
			alignment_data const& data, t_seq_bases const& ref, t_seq_bases const& alt);

	template<typename abort_checker>
	void run(abort_checker const& ac);

	t_haplotype get_haplotype_chains() const {
		return haplo;
	}

	t_strands get_chain_strands() const;

};

//TODO make it possible to switch between max_logsum + prior and logsum
inline t_loglike_vector haplo_chain_optimizer::compute_logsum(
		t_indices const& reads, t_position const pos) const {

	t_loglike_vector max_loglike(2);
	max_loglike(0) = -1 * std::numeric_limits<float>::infinity();
	max_loglike(1) = -1 * std::numeric_limits<float>::infinity();

	for (t_epi_base g = 0; g < 6; ++g) {

		//Ref prior
		double prior = g % 4 + 1 == ref(pos) || g % 4 + 1 == alt(pos) ? ref_match_prior_log : ref_no_match_prior_log;

		double loglike_fwd = prior;
		double loglike_rev = prior;

		for (t_index i = 0; i < reads.n_elem; ++i) {

			loglike_fwd += data.loglike_terms(reads(i), 0)(
					pos - read_start_postion(reads(i)), g);
			loglike_rev += data.loglike_terms(reads(i), 1)(
					pos - read_start_postion(reads(i)), g);
		}

		if (loglike_fwd > max_loglike(0)) {
			max_loglike(0) = loglike_fwd;
		}

		if (loglike_rev > max_loglike(1)) {
			max_loglike(1) = loglike_rev;
		}
	}

	if (!max_loglike.is_finite()) {
		throw std::runtime_error("compute_max_logsum : result not finite");
	}

	return max_loglike;

}

// DO NOT REMOVE
//inline t_loglike_vector haplo_chain_optimizer::compute_logsum(
//		t_indices const& reads, t_position const pos) const {
//
//	t_loglike_vector logsum(2, arma::fill::zeros);
//
//	for (t_epi_base g = 0; g < 6; ++g) {
//
//		double prod_fwd = 1;
//		double prod_rev = 1;
//
//		for (t_index i = 0; i < reads.n_elem; ++i) {
//
//			prod_fwd *= data.like_terms(reads(i), 0)(
//					pos - read_start_postion(reads(i)), g);
//			prod_rev *= data.like_terms(reads(i), 1)(
//					pos - read_start_postion(reads(i)), g);
//		}
//
//		logsum(0) += prod_fwd;
//		logsum(1) += prod_rev;
//
//	}
//
//	logsum = log(logsum);
//
//	if (!logsum.is_finite()) {
//		throw std::runtime_error("compute_logsum : result not finite");
//	}
//
//	return logsum;
//
//}


inline double haplo_chain_optimizer::compute_haplochain_prior(
		t_count const unique_chains) const {

	if (unique_chains >= haplochain_log_prior.n_elem || unique_chains <= 0) {
		throw std::runtime_error("compute_haplochain_prior : invalid input");
	}

	return haplochain_log_prior(unique_chains - 1);
}

//TODO if we use logsum instead of max_logsum then we do not need ref and alt
inline haplo_chain_optimizer::haplo_chain_optimizer(
		AlgorithmConfiguration const& config, alignment_data const& data, t_seq_bases const& ref, t_seq_bases const& alt) :
		max_iterations(config.max_iterations), min_overlap_length(config.min_overlap_length),
		haplo(data.n_reads), data(data), ref(ref), alt(alt), haplochain_log_prior(10000) {

	///////////////// Initialize:

	///// Init haplo chains

	//haplo.zeros();
	for(t_index i = 0; i < haplo.n_elem; ++i) {
		haplo(i) = i;
	}

	///// Init haplochain priors

	//geometric distribution
	haplochain_log_prior(0) = config.haplo_prior;
	for (t_count i = 1; i < haplochain_log_prior.n_elem; ++i) {
		haplochain_log_prior(i) = (1 - config.haplo_prior) * haplochain_log_prior(i - 1);
	}

	haplochain_log_prior = log(haplochain_log_prior);

	//Init ref prior
	ref_no_match_prior_log = log(1-config.ref_prior);
	ref_match_prior_log = log(config.ref_prior);

}

template<typename abort_checker>
inline void haplo_chain_optimizer::run(const abort_checker& ac) {

	TIMER_START

	t_count change_count = 1;

	while(change_count != 0) {
		change_count = optimize_profile(ac);
		change_count += optimize_haplochains(ac);
	}
}

//Find overlapping haplo chains
inline t_haplotype haplo_chain_optimizer::compute_feasible_haplotypes(t_position const start, t_position const end) const {

	t_index const read_number_start = min(data.read_numbres(start));
	t_index const read_number_end = max(data.read_numbres(end));

	t_haplotype feasible_haplotypes;

	for (t_index read = read_number_start; read <= read_number_end; read++) {

		//Check that read overlap haplo chain
		if (start <= read_start_postion(read)) {

			t_position overlap_length = 1 + min(end, read_end_postion(read)) - max(start, read_start_postion(read));

			//TODO this is a debug check -- should not be active in release versions
			if(overlap_length < 1) {
				//TODO remove
				//cout << overlap_length << endl;
				//cout << read_end << " : " << read_end_postion(read) << " : " << read_start << " : " << read_start_postion(read) << endl;
				throw std::runtime_error("Internal error");
			}

			if(overlap_length >= min_overlap_length) {
				feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
				feasible_haplotypes(feasible_haplotypes.n_elem - 1) = haplo(read);
			}
		}

	}

	return unique(feasible_haplotypes);
}

template<typename abort_checker>
inline t_count haplo_chain_optimizer::optimize_profile(
		const abort_checker& ac) {

	t_count i = 0;
	for (; i < max_iterations; i++) {

		if (ac.is_aborted()) {
			return 0;
		}

		t_count changes = 0;

		for (t_index read_number = 0; read_number < data.n_reads;
				++read_number) {

			//Check for ctrl C
			if (ac.check_abort()) {
				break;
			}

			t_haplotype feasible_haplotypes = compute_feasible_haplotypes(read_start_postion(read_number), read_end_postion(read_number));

			//Add free haplo chain , i.e. a haplo chain with no reads
			//TODO take first free haplotype
			feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
			feasible_haplotypes(feasible_haplotypes.n_elem - 1) = max(haplo) + 1;

			t_loglike_vector loglike = profile_posterior(read_number, feasible_haplotypes);

			if (!is_finite(loglike)) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index i = argmax(loglike);

			//check if any improvement
			if (1e-10 >= loglike(i)) {
				continue;
			}

			//TODO remove
//			cout << " ------------------------------- " << endl;
//			cout << read_number << " : " << compute_posterior(haplo) << endl;
//			cout << haplo(read_number) << " -> " << feasible_haplotypes(i) << endl;
//			cout << trans(feasible_haplotypes) << endl;
//			cout << trans(loglike) << endl;

			//Get current haplotype chain of read
			t_haplochain old_chain = haplo(read_number);

			//Update haplotype chain of read
			haplo(read_number) = feasible_haplotypes(i);

			//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
			update_haplotype_chain(haplo, old_chain);

			//add change
			changes++;

		}

//TODO remove
//		t_haplotype u = unique(haplo);
//		cout << "ch : " << changes << " : " << max(haplo) << " : " << u.n_elem
//				<< " : " << compute_posterior(haplo) << endl;

		if (changes == 0) {
			break;
		}
	}

	if (i == max_iterations) {
		report_error("Max iteration limit reached");
	}

	return i;

}

template<typename abort_checker>
inline t_count haplo_chain_optimizer::optimize_haplochains(
		const abort_checker& ac) {

	t_count i = 0;
	for (; i < max_iterations; i++) {

		if (ac.is_aborted()) {
			return 0;
		}

		t_count changes = 0;

		t_haplotype unique_chains = unique(haplo);

		for (t_index j = 0; j < unique_chains.n_elem; ++j) {

			//Check for ctrl C
			if (ac.check_abort()) {
				break;
			}

			t_haplochain chain = unique_chains(j);

			//Find chain start and end position
			t_indices reads_in_chain = find(haplo == chain);

			t_position chain_start = read_start_postion(reads_in_chain(0));
			t_position chain_end = max(data.reads_end_postions(reads_in_chain));

			t_haplotype feasible_haplotypes = compute_feasible_haplotypes(chain_start, chain_end);

			if(feasible_haplotypes.is_empty()) {
				continue;
			}

			t_loglike_vector posterior = profile_posterior(reads_in_chain, feasible_haplotypes);

			//TODO remove
			//cout << chain_start << " : " << trans(feasible_haplotypes) << endl;
			//cout << chain_start << " : " << trans(posterior) << endl;

			if (!is_finite(posterior)) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index i = argmax(posterior);

			//check if any improvement
			if (1e-10 >= posterior(i)) {
				continue;
			}

			//Update haplotype chain of read
			haplo(reads_in_chain).fill(feasible_haplotypes(i));

			//add change
			changes++;
		}

		//TODO remove
		//t_haplotype u = unique(haplo);
		//cout << "ch : " << changes << " : " << max(haplo) << " : " << u.n_elem
		//				<< " : " << compute_posterior(haplo) << endl;

		if (changes == 0) {
			break;
		}
	}

	if (i == max_iterations) {
		report_error("Max iteration limit reached");
	}

	return i;

}


void haplo_chain_optimizer::update_haplotype_chain(t_haplotype & haplotype, t_haplochain const old_chain) const {

	if (sum(haplotype == old_chain) == 0) {
		//empty chain
		return;
	}

	t_indices reads_in_chain = find(haplotype == old_chain);

	t_haplochain chain = old_chain;
	t_position chain_end = read_end_postion(reads_in_chain(0));

	for (t_index i = 1; i < reads_in_chain.n_elem; ++i) {

		t_index read = reads_in_chain(i);

		if (read_start_postion(read) > chain_end - min_overlap_length + 1) {
			//reads not overlapping => new chain
			chain = max(haplotype) + 1;
		}

		haplotype(read) = chain;
		chain_end = read_end_postion(read);
	}
}

double haplo_chain_optimizer::compute_prior(t_haplotype const& h, t_position const start, t_position const end) const {

	double prior = 0;

	for (t_position pos = start; pos < end; ++pos) {

			if (data.read_numbres(pos).n_elem == 0) {
				continue;
			}

			t_haplotype unique_chains = unique(h(data.read_numbres(pos)));
			prior += compute_haplochain_prior(unique_chains.n_elem);
		}

	return prior;
}

double haplo_chain_optimizer::compute_chain_loglike(t_haplotype const & h,
		t_haplochain const chain) const {

	t_indices reads_in_chain = find(h == chain);

	if (reads_in_chain.n_elem == 0) {
		return 0;
	}

	//find start, end pos of chain
	t_position chain_start = read_start_postion(min(reads_in_chain));
	t_position chain_end = read_end_postion(reads_in_chain(0));

	for (t_index j = 1; j < reads_in_chain.n_elem; ++j) {
		if (read_end_postion(reads_in_chain(j)) > chain_end) {
			chain_end = read_end_postion(reads_in_chain(j));
		}
	}

	t_loglike_vector logsum(2, arma::fill::zeros);

	for (t_position pos = chain_start; pos <= chain_end; ++pos) {

		t_indices reads = get_overlapping_reads_in_chain(h, pos, chain);

		//Vector of size 2 logsum of (fwd, rev)
		logsum += compute_logsum(reads, pos);

	}

	return max(logsum);
}

int haplo_chain_optimizer::chain_strand(
		t_haplochain const chain) const {

	t_indices reads_in_chain = find(haplo == chain);

	if (reads_in_chain.n_elem == 0) {
		return 0;
	}

	//find start, end pos of chain
	t_position chain_start = read_start_postion(min(reads_in_chain));
	t_position chain_end = read_end_postion(reads_in_chain(0));

	for (t_index j = 1; j < reads_in_chain.n_elem; ++j) {
		if (read_end_postion(reads_in_chain(j)) > chain_end) {
			chain_end = read_end_postion(reads_in_chain(j));
		}
	}

	t_loglike_vector logsum(2, arma::fill::zeros);

	for (t_position pos = chain_start; pos <= chain_end; ++pos) {

		t_indices reads = get_overlapping_reads_in_chain(haplo, pos, chain);

		//Vector of size 2 logsum of (fwd, rev)
		logsum += compute_logsum(reads, pos);

	}

	if(logsum(0) > logsum(1)) {
		return strand_fwd;
	}

	return strand_rev;
}

t_strands haplo_chain_optimizer::get_chain_strands() const {

	t_strands strands(haplo.n_elem);

	for(t_index i = 0; i < haplo.n_elem; ++i) {
		strands(i) = chain_strand(haplo(i));
	}

	return strands;
}


double haplo_chain_optimizer::compute_posterior(t_haplotype const& h) const {

	//Compute prior

	double posterior = compute_prior(h, 0, data.sequence_length);

	// Compute likelihood
	t_haplotype unique_chains = unique(h);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {

		posterior += compute_chain_loglike(h, unique_chains(i));

	}

	return posterior;
}

double haplo_chain_optimizer::compute_delta_posterior(t_index const read_number, t_haplochain const new_chain) const {

	t_haplotype new_haplo(haplo); //TODO this is alot of copying
	t_haplochain current_chain = haplo(read_number);

	t_position start = haplo_chain_start(current_chain);
	t_position end = haplo_chain_end(current_chain);

	if(current_chain == new_chain) {
		return 0;
	}

	double delta_posterior = 0;

	//Compute posterior
	new_haplo(read_number) = new_chain;
	update_haplotype_chain(new_haplo, current_chain);

	delta_posterior += compute_prior(new_haplo, start, end);
	delta_posterior -= compute_prior(haplo, start, end);

	delta_posterior += compute_chain_loglike(new_haplo, new_chain);
	delta_posterior -= compute_chain_loglike(haplo, new_chain);

	delta_posterior += compute_chain_loglike(new_haplo, current_chain);
	delta_posterior -= compute_chain_loglike(haplo, current_chain);

	for(t_haplochain chain = max(haplo) + 1; chain <= max(new_haplo); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		delta_posterior += compute_chain_loglike(new_haplo, chain);
	}

	//TODO debug guards
//	if(abs(delta_posterior - compute_posterior(new_haplo) + compute_posterior(haplo)) > 1e-10) {
//		throw std::runtime_error("compute_delta_posterior : Error");
//	}

	return delta_posterior;
}

//Merge old_chain = (reads_in_chain) and new_chain
//TODO rename - chain merge
double haplo_chain_optimizer::compute_delta_posterior(t_indices const& reads_in_chain, t_haplochain const new_chain) const {

	t_position start = data.reads_start_postions(reads_in_chain(0));
	t_position end =  max(data.reads_end_postions(reads_in_chain));

	t_haplotype new_haplo(haplo);
	t_haplochain current_chain = haplo(reads_in_chain(0));

	if(current_chain == new_chain) {
		return 0;
	}

	double posterior = 0;

	//current posterior
	posterior -= compute_prior(haplo, start, end);
	posterior -= compute_chain_loglike(haplo, current_chain);
	posterior -= compute_chain_loglike(haplo, new_chain);

	//new posterior
	new_haplo(reads_in_chain).fill(new_chain);

	posterior += compute_prior(new_haplo, start, end);
	posterior += compute_chain_loglike(new_haplo, current_chain);
	posterior += compute_chain_loglike(new_haplo, new_chain);

	return posterior;
}

t_loglike_vector haplo_chain_optimizer::profile_posterior(
		t_index const read_number,
		t_haplotype const& feasible_haplotypes) const {

	t_loglike_vector loglike(feasible_haplotypes.n_elem, arma::fill::zeros);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(i) = compute_delta_posterior(read_number, feasible_haplotypes(i));
	}

	return loglike;
}

t_loglike_vector haplo_chain_optimizer::profile_posterior(t_indices const& reads_in_chain, t_haplotype const& feasible_haplotypes) const {

	t_loglike_vector loglike(feasible_haplotypes.n_elem, arma::fill::zeros);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(i) = compute_delta_posterior(reads_in_chain, feasible_haplotypes(i));
	}

	return loglike;
}


inline t_indices haplo_chain_optimizer::get_overlapping_reads_in_chain(t_haplotype const& h,
		t_position const pos, t_haplochain const chain) const {

	t_indices reads_at_pos = data.read_numbres(pos);
	arma::uvec chain_indicator = (h(reads_at_pos) == chain);

	if (sum(chain_indicator) == 0) {
		//Return empty vector
		return t_indices();
	}

	//Return read numbers of reads in chain
	return reads_at_pos(find(chain_indicator));
}

inline t_position haplo_chain_optimizer::read_start_postion(
		t_index const read_number) const {
	return data.reads_start_postions(read_number);
}

inline t_position haplo_chain_optimizer::read_end_postion(
		t_index const read_number) const {
	return data.reads_end_postions(read_number);
}

#endif /* HAPLO_CHAIN_OPTIMIZER_HPP_ */
