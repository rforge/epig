/*
 * alignment_data.hpp
 *
 *  Created on: Jan 26, 2014
 *      Author: martin
 */

#ifndef ALIGNMENT_DATA_HPP_
#define ALIGNMENT_DATA_HPP_

class alignment_data { //TODO better name for this class

public:

	t_count const n_reads;

	t_models const fwd_model;
	t_models const rev_model;

	t_positions const reads_start_postions;
	t_positions const reads_end_postions;

	t_position const offset;
	t_position const sequence_length; //length of sequence

	t_base_subsets const base_subsets; //vector of size (length of sequence)

	arma::field<arma::mat> const loglike_terms; //field of size (number of reads) x 2
	arma::field<arma::mat> const like_terms; //field of size (number of reads) x 2

	t_loglike_vector const log_pmax; //vector of size (number of reads)

	field<t_indices> const read_numbres; //Reads overlapping position

    alignment_data(AlgorithmConfiguration const& config, std::vector<aligned_read> const& reads, t_count hard_limit);

};

inline alignment_data::alignment_data(AlgorithmConfiguration const& config, std::vector<aligned_read> const& reads, t_count hard_limit) :
        n_reads(min(static_cast<t_count>(reads.size()), hard_limit)), fwd_model(config.fwd_model), rev_model(config.rev_model), offset(0), sequence_length(0) {

	if(reads.empty()) {
		throw std::runtime_error("No reads");
	}

	t_positions starts(n_reads);
	t_positions ends(n_reads);

	for (t_count i = 0; i < n_reads; ++i) {

		aligned_read read = reads[i];

		starts(i) = read.position;
		ends(i) = read.position + read.length - 1;
	}

	//Set variables
	const_cast<t_position&>(this->offset) = starts.min();
	const_cast<t_position&>(this->sequence_length) = ends.max() - starts.min() + 1;

	const_cast<t_positions&>(this->reads_start_postions) = starts - offset;
	const_cast<t_positions&>(this->reads_end_postions) = ends - offset;

	t_counts coverage(sequence_length, arma::fill::zeros);
	for (t_count i = 0; i < n_reads; ++i) {
		coverage.subvec(reads_start_postions(i), reads_end_postions(i)) = coverage.subvec(reads_start_postions(i), reads_end_postions(i)) + 1;
	}

	t_base_subsets subsets(sequence_length, arma::fill::zeros);
	arma::field<arma::mat> tmp_loglike_terms(n_reads, 2);
	arma::field<arma::mat> tmp_like_terms(n_reads, 2);

	arma::field<t_indices> tmp_read_numbers(sequence_length);

	t_loglike_vector tmp_pmax(n_reads);

	t_counts pos_read_count(sequence_length, arma::fill::zeros);

	for (t_count i = 0; i < n_reads; ++i) {

		aligned_read read = reads[i];
		t_length L = read.length;

		//Check L <= max length
		if(L > fwd_model.n_elem) {
			throw std::runtime_error("Read exceeds max length");
		}

		tmp_loglike_terms(i, 0).set_size(L, 6);
		tmp_loglike_terms(i, 1).set_size(L, 6);

		tmp_like_terms(i, 0).set_size(L, 6);
		tmp_like_terms(i, 1).set_size(L, 6);

		for (t_position j = reads_start_postions(i);
				j <= reads_end_postions(i); ++j) {

			//read numbers
			if (tmp_read_numbers(j).is_empty()) {
				tmp_read_numbers(j).set_size(coverage(j));
			}

			tmp_read_numbers(j)(pos_read_count(j)) = i;
			++pos_read_count(j);

			//log sums
			t_position k = j - reads_start_postions(i);

			t_base base = read.bases(k);
			double epsilon = read.epsilon(k) * config.sequence_quality_adjust; //FIXME configable

			if(base == 0) {
				//base = N
				tmp_loglike_terms(i, 0).row(k).fill(log(0.25));
				tmp_loglike_terms(i, 1).row(k).fill(log(0.25));

				tmp_like_terms(i, 0).row(k).fill(0.25);
				tmp_like_terms(i, 1).row(k).fill(0.25);

			} else {
				//fwd strand
				tmp_loglike_terms(i, 0).row(k) = log((1-4/static_cast<double>(3)*epsilon)*fwd_model(k).row(base-1) + epsilon/static_cast<double>(3));

				//rev strand
				tmp_loglike_terms(i, 1).row(k) = log((1-4/static_cast<double>(3)*epsilon)*rev_model(L-1-k).row(base-1) + epsilon/static_cast<double>(3));

				tmp_like_terms(i, 0).row(k) = (1-4/static_cast<double>(3)*epsilon)*fwd_model(k).row(base-1) + epsilon/static_cast<double>(3);
				tmp_like_terms(i, 1).row(k) = (1-4/static_cast<double>(3)*epsilon)*rev_model(L-1-k).row(base-1) + epsilon/static_cast<double>(3);
			}

			//Position marks
			subsets(j) = base_union(subsets(j), base);
		}

		tmp_pmax(i) = max(as_scalar(mean(max(tmp_loglike_terms(i,0), 1))), as_scalar(mean(max(tmp_loglike_terms(i,1), 1))));
	}

	//Set variables
	const_cast<t_base_subsets&>(this->base_subsets) = subsets;
	const_cast<arma::field<arma::mat>&>(this->loglike_terms) = tmp_loglike_terms;
	const_cast<arma::field<arma::mat>&>(this->like_terms) = tmp_like_terms;
	const_cast<field<t_indices>&>(this->read_numbres) = tmp_read_numbers;
	const_cast<t_loglike_vector&>(this->log_pmax) = tmp_pmax;

}




#endif /* ALIGNMENT_DATA_HPP_ */
