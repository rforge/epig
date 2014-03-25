/*
 * openmp_optimizer.hpp
 *
 *  Created on: Nov 15, 2013
 *      Author: martin
 */

#ifndef CHUNK_OPTIMIZER_HPP_
#define CHUNK_OPTIMIZER_HPP_

using namespace arma;

class chunk_optimizer {

public:
	t_count const number_of_chunks;


	std::string refGenom_filename;
	std::string altGenom_filename;
	std::string const refName;

	t_position const offset;

	t_count const max_threads;

	AlgorithmConfiguration const config;

	t_prior_vector priors;
	t_prior_vector priors_alt;

private:
	std::string bam_file;

	field<t_methylome> methylome; //coding 0 -N, 1 - C, 2 - G, 3 - A, 4 - T, 5 - c, 6 - g

	field<t_strands> strands; // one for each chunk
	field<t_alleles> alleles; // one for each chunk
	field<t_strands> unique_assigned; //TODO type
	field<field<t_indices> > read_ids; // one for each chunk
	field<t_allele_occupancy> allele_occupancy; // one for each chunk

	t_positions chunk_start_pos;
	t_positions chunk_end_pos;

    long unsigned int n_reads_hard_limit;

public:

    chunk_optimizer(std::string const& bam_file, std::string const& refGenom_filename, std::string const& altGenom_filename,
                    std::string const& refName, t_positions const& chunk_start_positions, t_positions const& chunk_end_positions,
			 t_count max_threads, AlgorithmConfiguration const& config);

	void run();

	field<t_methylome> get_methylome() const {
		return methylome;
	}

	field<t_alleles> get_alleles() const {
		return alleles;
	}

	field<t_strands> get_strands() const {
		return strands;
	}

	field<t_strands> get_unique_assigned() const {
		return unique_assigned;
	}

	field<t_allele_occupancy> get_allele_occupancy() const {
		return allele_occupancy;
    }

	field<field<t_indices> > get_readIDs() const {
		return read_ids;
	}

	t_positions get_chunks_start() const {
		return chunk_start_pos;
	}

	t_positions get_chunks_end() const {
		return chunk_end_pos;
	}
};

chunk_optimizer::chunk_optimizer(std::string const& bam_file, std::string const& refGenom_filename, std::string const& altGenom_filename,
        std::string const& refName, t_positions const& chunk_start_positions, t_positions const& chunk_end_positions,
		t_count max_threads, AlgorithmConfiguration const& config) :
		number_of_chunks(chunk_start_positions.n_elem),
		refGenom_filename(refGenom_filename), altGenom_filename(altGenom_filename),
		refName(refName), offset(chunk_start_positions.min()), max_threads(max_threads), config(config),
				priors(config.log_prior(0)), priors_alt(config.log_prior(1)),
				bam_file(bam_file), methylome(number_of_chunks), strands(
				number_of_chunks), alleles(number_of_chunks), unique_assigned(number_of_chunks), read_ids(
				number_of_chunks), allele_occupancy(number_of_chunks),
                chunk_start_pos(chunk_start_positions), chunk_end_pos(chunk_end_positions), n_reads_hard_limit(config.reads_hard_limit) {

    //TODO domain check consistency of chunk positions

	//Check refs avaialbe and build fai index files if needed
	read_fasta(refGenom_filename, refName, 1, 2);
	read_fasta(altGenom_filename, refName, 1, 2);
}


void chunk_optimizer::run() {

	// create progress monitor
	Progress p(number_of_chunks, config.verbose);

    //Warnings
    omp_rwarn warnings;

#ifdef EPIG_USE_OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel for schedule(dynamic)
#endif
	for (t_count i = 0; i < number_of_chunks; ++i) {

		if (!p.is_aborted()) {

			bamReader reader(bam_file, refName, chunk_start_pos(i),
                    chunk_end_pos(i));

#ifdef EPIG_USE_OPENMP
#pragma omp critical
#endif
           {
            //TODO does this need to be in a critical region ??
            //TODO use har limit when tetching reads

			reader.fetch(); //Load reads

			//Check that chunk is non empty
			if (reader.get_reads().empty()) {
                    methylome(i).set_size(config.max_alleles, chunk_end_pos(i) - chunk_start_pos(i) + 1);
					methylome(i).zeros();
					strands(i).set_size(0);
					alleles(i).set_size(0);
					unique_assigned(i).set_size(0);
					read_ids(i).set_size(0);
					allele_occupancy(i).set_size(config.max_alleles + 1, chunk_end_pos(i) - chunk_start_pos(i) + 1);
					allele_occupancy(i).zeros();
				}

			}

            if(reader.get_reads().empty()) {
                p.increment();
                continue;
            }

            if(reader.get_reads().size() > n_reads_hard_limit) {
                std::ostringstream msg;
                msg << reader.get_reads().size() << " reads are covering chunk (position: " << chunk_start_pos(i) << " : " << chunk_end_pos(i) << ") this exceeds the hard limit. Only the first " << n_reads_hard_limit << " reads are used.";
                warnings.add(msg.str());
            }

            std::vector<aligned_read> const& reads = reader.get_reads();
            alignment_data data(config, reads, n_reads_hard_limit);

			//Load refs
            std::string refGenom = read_fasta(refGenom_filename, refName, data.offset,
					data.sequence_length);
            std::string altGenom = read_fasta(altGenom_filename, refName, data.offset,
					data.sequence_length);

			if (refGenom.size() != static_cast<unsigned int>(data.sequence_length)) {
				throw std::runtime_error("Problem with refGenom"); //TODO error msg
			}

			if (altGenom.size() != static_cast<unsigned int>(data.sequence_length)) {
				throw std::runtime_error("Problem with altGenom"); //TODO error msg
			}

			optimizer opt(config, data, refGenom, altGenom);

			opt.run(p);

#ifdef EPIG_USE_OPENMP
#pragma omp critical
#endif
			{

				chunk_start_pos(i) = max(chunk_start_pos(i),
						data.reads_start_postions.min() + data.offset);
				chunk_end_pos(i) = min(chunk_end_pos(i),
						data.reads_end_postions.max() + data.offset);

                methylome(i) = opt.get_methylome().cols(
                        chunk_start_pos(i) - data.offset,
                        chunk_end_pos(i) - data.offset);

                strands(i) = opt.get_strands();
                alleles(i) = opt.get_alleles();
                unique_assigned(i) = opt.get_unique_assigned();

                read_ids(i) = data.read_numbres.subfield(
                        chunk_start_pos(i) - data.offset, 0,
                        chunk_end_pos(i) - data.offset, 0);

                allele_occupancy(i) = opt.get_allele_occupancy().cols(
                        chunk_start_pos(i) - data.offset,
                        chunk_end_pos(i) - data.offset);
            }

            p.increment();
		}
	}

}

chunk_optimizer create_base_chunk_optimizer(std::string bam_file, std::string refGenom_filename, std::string altGenom_filename,
            std::string refName, t_position start_position, t_position end_position, t_count max_threads,
            t_position chunk_size, AlgorithmConfiguration const& config) {

            t_count number_of_chunks = ceil((end_position - start_position + 1)/ static_cast<double>(chunk_size));

            t_positions chunk_start_pos(number_of_chunks);
            t_positions chunk_end_pos(number_of_chunks);

        //Compute chunk positions
        chunk_start_pos(0) = start_position;
        chunk_end_pos(0) = chunk_start_pos(0) + chunk_size - 1;

        for (t_count i = 1; i < number_of_chunks; ++i) {
            chunk_start_pos(i) = chunk_end_pos(i - 1) + 1;
            chunk_end_pos(i) = chunk_start_pos(i) + chunk_size - 1;
        }

        //Correct last end position
        chunk_end_pos(number_of_chunks - 1) = end_position;

        return chunk_optimizer(bam_file, refGenom_filename, altGenom_filename, refName, chunk_start_pos, chunk_end_pos, max_threads, config);
    }


#endif /* CHUNK_OPTIMIZER_HPP_ */
