#ifndef GENOTYPE_OPTIMIZER_HPP_
#define GENOTYPE_OPTIMIZER_HPP_

class genotype_optimizer {

	//*** Fixed model parameters
	t_haplotype const haplo; //vector of size n_reads
	t_strands const& strands; // vector of size n_reads

	//*** Data
	alignment_data const& data;
	field<t_genotype> achievable_epi_bases; //field of size 15 (number of possible subsets of bases)

	//*** Refs
	t_seq_bases const& ref;
	t_seq_bases const& alt;

	// Prior
	reference_genome_prior prior;

	//TODO note duplicate in haplo_chain_optimizer
	t_position haplo_chain_start(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		return data.reads_start_postions(min(reads_in_chain));
	}

	//TODO note duplicate in haplo_chain_optimizer
	t_count haplo_chain_end(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		t_positions end_pos = data.reads_end_postions(reads_in_chain);
		return max(end_pos);
	}

	t_strand haplo_chain_strand(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain, 1));
		return strands(reads_in_chain(0));
	}

	t_strands haplo_chain_strand(t_haplotype const chains) const {

		t_strands s(chains.n_elem);

		for(t_index i = 0; i < chains.n_elem; ++i) {
			s(i) = haplo_chain_strand(chains(i));
		}

		return s;
	}

	double compute_genotype_loglike(t_genotype const& genotype,
			t_haplotype const& haplo_at_pos, t_indices const& reads_at_pos,
			t_position const pos) const {

		double loglike = 0;
		for (t_index i = 0; i < reads_at_pos.n_elem; ++i) {
			loglike += compute_loglike_of_read(reads_at_pos(i), haplo_at_pos, genotype, pos);
		}

		return loglike;
	}

	double compute_loglike_of_read(t_index const read_number, t_haplotype const& haplo_at_pos, t_genotype const& genotype, t_position const pos) const {

		t_indices haplo_index(find(haplo_at_pos == haplo(read_number),1));
		t_epi_base g = genotype(haplo_index(0));

		return data.loglike_terms(read_number, strands(read_number))(pos - data.reads_start_postions(read_number), g);
	}

	double compute_genotype_prior(t_genotype const& genotype, t_counts const& read_counts, t_strands const& strands, t_position const pos) const {

		return prior.getPrior(pos, genotype, read_counts, strands);
	}

	//return true if changes where made
	bool fit_genotype_at_pos(t_position const pos, t_haplotype const& haplo_at_pos, t_indices const& haplo_to_opt) {

		t_indices reads_at_pos = data.read_numbres(pos);

		t_strands haplo_strands = haplo_chain_strand(haplo_at_pos);

		if(haplo_to_opt.is_empty()) {
			return false;
		}

		//Get current genotype and read count
		t_counts read_counts(haplo_at_pos.n_elem);
		t_genotype current_genotype(haplo_at_pos.n_elem);

		for(t_index i = 0; i < haplo_at_pos.n_elem; ++i) {
			current_genotype(i) = genotypes(haplo_at_pos(i))(pos - haplo_chain_start(haplo_at_pos(i)));
			read_counts(i) = sum(haplo(reads_at_pos) == haplo_at_pos(i));
		}

		t_genotype aeb(achievable_epi_bases(data.base_subsets(pos)));
		base_counter g(aeb.n_elem, haplo_to_opt.n_elem);

		//Compute old posterior
		double old_posterior = compute_genotype_prior(current_genotype, read_counts, haplo_strands, pos);
		old_posterior += compute_genotype_loglike(current_genotype, haplo_at_pos, reads_at_pos, pos);

		bool change_made = false;

		g.set_zero();
		do {

			t_genotype genotype(current_genotype);
			genotype(haplo_to_opt) = aeb(g.get_number());

			double posterior = compute_genotype_prior(genotype, read_counts, haplo_strands, pos);
			posterior += compute_genotype_loglike(genotype, haplo_at_pos, reads_at_pos, pos);

			if(posterior > old_posterior) {

				old_posterior = posterior;
				change_made = true;

				//TODO we only need to update the genotypes of haplo_to_opt
				for(t_index i = 0; i < haplo_at_pos.n_elem; ++i) {
					genotypes(haplo_at_pos(i))(pos - haplo_chain_start(haplo_at_pos(i))) = genotype(i);
				}

			}

			++g;

		} while(!g.is_zero());


		return change_made;
	}

	void fit_genotype_at_pos(t_position const pos) {

		t_indices reads_at_pos = data.read_numbres(pos);
		t_haplotype haplo_at_pos = unique(haplo(reads_at_pos));

		if(haplo_at_pos.is_empty()) {
			return;
		}

		if(haplo_at_pos.n_elem <= 2) {

			t_indices haplo_to_opt;
			haplo_to_opt << 0 << 1 << endr;

			fit_genotype_at_pos(pos, haplo_at_pos, haplo_to_opt.subvec(0,haplo_at_pos.n_elem-1));

			return;
		}

		if(haplo_at_pos.n_elem <= 3) {

			t_indices haplo_to_opt;
			haplo_to_opt << 0 << 1 << 2 << endr;

			fit_genotype_at_pos(pos, haplo_at_pos, haplo_to_opt.subvec(0,haplo_at_pos.n_elem-1));

			return;
		}

		int changes;
		do {
			t_indices haplo_to_opt;
			haplo_to_opt << 0 << 1 << 2 << 3 << endr; //TODO number of haplos to opt configable

			changes = 0;
			for(t_index i = 0; i < haplo_at_pos.n_elem - 3; ++i) {
				changes += fit_genotype_at_pos(pos, haplo_at_pos, haplo_to_opt);
				haplo_to_opt = haplo_to_opt + 1;
			}

			//TODO configable
			if(changes != 0) {
				haplo_at_pos = shuffle(haplo_at_pos);
			}

		} while(changes != 0);

	}


	bool is_base_achievable(t_epi_base const base, t_base_subset const p) const {

		return true;

		if (p == 0) {
			return false;
		}

		switch (base) {
		case 0:
			if (p == 2 || p == 3 || p == 5) {
				return false;
			}

			return true;

		case 1:
			if (p == 1 || p == 4 || p == 6) {
				return false;
			}

			return true;

		case 2:
			if (p == 1 || p == 2 || p == 4 || p == 6 || p == 7 || p == 9
					|| p == 11) {
				return false;
			}

			return true;

		case 3:
			if (p == 1 || p == 2 || p == 3 || p == 5 || p == 8 || p == 9
					|| p == 14) {
				return false;
			}

			return true;

		case 4:
			if (p == 2 || p == 3 || p == 4 || p == 5 || p == 7 || p == 10
					|| p == 12) {
				return false;
			}

			return true;

		case 5:
			if (p == 1 || p == 3 || p == 4 || p == 6 || p == 8 || p == 10
					|| p == 13) {
				return false;
			}

			return true;
		default:
			throw std::runtime_error("internal error");
			break;
		}

	}

public:

	// genotype
	arma::field<arma::Col<t_epi_base> > genotypes; //field of size = size of unique(haplo)

	genotype_optimizer(alignment_data const& data, t_haplotype const& haplo_unfactored,
			t_strands const& strands, t_seq_bases const& ref,
			t_seq_bases const& alt, AlgorithmConfiguration const& config) :
			haplo(factor(haplo_unfactored)), strands(strands), data(data),
			achievable_epi_bases(16), ref(ref), alt(alt),
			prior(ref, alt, config),
			genotypes(max(haplo)+1)
	{
		//note haplo chains ids are 0, 1, ..., number of chains - 1

		for(t_index i = 0; i < genotypes.n_elem; ++i) {
			genotypes(i).set_size(haplo_chain_end(i) - haplo_chain_start(i)+1);
			genotypes(i).zeros();
		}


		for (t_base_subset subset = 0; subset <= 15; ++subset) {
			for(t_epi_base g = 0; g < 6; ++g) {
				if(is_base_achievable(g, subset)) {
					achievable_epi_bases(subset).set_size(achievable_epi_bases(subset).n_elem+1);
					achievable_epi_bases(subset)(achievable_epi_bases(subset).n_elem-1) = g;
				}
			}
		}

	}

	void fit() {

		for(t_position pos = 0; pos < data.sequence_length; ++pos)  {
			fit_genotype_at_pos(pos);
		}

	}

	t_positions haplo_chain_start() const {
		t_positions pos(max(haplo)+1);
		for(t_index i = 0; i < pos.n_elem; ++i) {
			pos(i) = haplo_chain_start(i);
		}

		return pos;
	}

	t_positions haplo_chain_end() const {
		t_positions pos(max(haplo)+1);
		for(t_index i = 0; i < pos.n_elem; ++i) {
			pos(i) = haplo_chain_end(i);
		}

		return pos;
	}

};

#endif /* GENOTYPE_OPTIMIZER_HPP_ */
