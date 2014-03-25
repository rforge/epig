/*
 * epiG_R_interface.hpp
 *
 *  Created on: Nov 11, 2013
 *      Author: martin
 */

#ifndef EPIG_R_INTERFACE_HPP_
#define EPIG_R_INTERFACE_HPP_

extern "C" {

SEXP r_epiG_fit_filename(SEXP r_filename, SEXP r_ref_genom_filename, SEXP r_alt_genom_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_max_threads, SEXP r_max_chunk_size, SEXP r_config);
SEXP r_epiG_fit_filename_chunks(SEXP r_filename, SEXP r_ref_genom_filename, SEXP r_alt_genom_filename, SEXP r_refName, SEXP r_chunks_start, SEXP r_chunks_end, SEXP r_max_threads, SEXP r_config);

SEXP r_epiG_fetch_reads(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end);
SEXP r_epiG_fetch_reads_info(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end);
SEXP r_epiG_compute_chunk_positions(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_chunk_size);

SEXP r_epiG_read_fasta(SEXP r_filename, SEXP r_ref, SEXP r_position, SEXP r_length);

}

SEXP epiG_read_fasta(SEXP r_filename, SEXP r_ref, SEXP r_position, SEXP r_length) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string ref = get_value<std::string>(r_ref);
	const int pos = get_value<int>(r_position);
	const int length = get_value<int>(r_length);

	//Note the first base is at postion 0
	return rObject(create_bases_vector(read_fasta(filename, ref, pos, length)));
}


SEXP r_epiG_read_fasta(SEXP r_filename, SEXP r_ref, SEXP r_position, SEXP r_length) {

	try {

		return epiG_read_fasta(r_filename, r_ref, r_position, r_length);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}

SEXP epiG_fetch_reads(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string refName = get_value<std::string>(r_refName);
	const t_position start = get_value<t_position>(r_start);
	const t_position end = get_value<t_position>(r_end);

	bamReader reader(filename, refName, start, end);
	reader.fetch();

    std::vector<aligned_read> const& reads = reader.get_reads();

	field<t_seq_bases> bases(reads.size());
	field<t_epsilon_quality> epsilon(reads.size());
	t_positions pos(reads.size());
	t_lengths len(reads.size());

    for(unsigned int i = 0; i < reads.size(); ++i) {

		aligned_read read = reads[i];

		bases(i) = read.bases;
		epsilon(i) = read.epsilon;
		pos(i) = read.position;
		len(i) = read.length;
	}

    rList res;

	res.attach(rObject(bases), "reads");
	res.attach(rObject(epsilon), "quality");
	res.attach(rObject(pos), "positions");
	res.attach(rObject(len), "lengths");

    return rObject(res);
}

SEXP r_epiG_fetch_reads(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

	try {

		return epiG_fetch_reads(r_filename, r_refName, r_start, r_end);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}

SEXP epiG_fetch_reads_info(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

    const std::string filename = get_value<std::string>(r_filename);
    const std::string refName = get_value<std::string>(r_refName);
    const t_position start = get_value<t_position>(r_start);
    const t_position end = get_value<t_position>(r_end);

    bamReader reader(filename, refName, start, end);
    reader.fetch_reads_info();
    std::vector<read_info> const& infos = reader.get_infos();

    t_positions pos(infos.size());
    t_lengths len(infos.size());

    for(unsigned int i = 0; i < infos.size(); ++i) {

        read_info info = infos[i];

        pos(i) = info.position;
        len(i) = info.length;
    }

    rList res;

    res.attach(rObject(pos), "positions");
    res.attach(rObject(len), "lengths");

    return rObject(res);
}

SEXP r_epiG_fetch_reads_info(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end) {

    try {

        return epiG_fetch_reads_info(r_filename, r_refName, r_start, r_end);

        //Catch unhandled exceptions

    } catch (std::exception & e) {

        if (e.what() != NULL) {
            report_error(e.what());
        }

        else {
            report_error("Unknown error");
        }

    } catch (...) {
        report_error("Unknown error");
    }

    return R_NilValue; //Avoid compiler warnings
}

SEXP epiG_compute_chunk_positions(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_chunk_size) {

    const std::string filename = get_value<std::string>(r_filename);
    const std::string refName = get_value<std::string>(r_refName);
    const t_position start = get_value<t_position>(r_start);
    const t_position end = get_value<t_position>(r_end);
    const t_count chunk_size = get_value<t_count>(r_chunk_size);

    bamReader reader(filename, refName, start, end);
    reader.fetch_reads_info();
    std::vector<read_info> const& infos = reader.get_infos();

    t_positions start_pos(infos.size());
    t_positions end_pos(infos.size());

    unsigned int n_reads = infos.size();

    for(unsigned int i = 0; i < n_reads; ++i) {

        read_info info = infos[i];

        start_pos(i) = info.position;
        end_pos(i) = info.position + info.length - 1;
    }

    std::vector<t_position> pos;
    pos.push_back(0);

    int end_sum = 0;

    for(unsigned int i = 0; i < n_reads; ++i) {

        if( i+1 - end_sum >= chunk_size) {
                pos.push_back(end_pos(i)+1);
                end_sum = sum(end_pos.subvec(0,i-1) < pos.back());
        }
    }

    //Correct ends
    pos[0] = start;
    pos.back() = end;

    return rObject(pos);
}

SEXP r_epiG_compute_chunk_positions(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_chunk_size) {

    try {

        return epiG_compute_chunk_positions(r_filename, r_refName, r_start, r_end, r_chunk_size);

        //Catch unhandled exceptions

    } catch (std::exception & e) {

        if (e.what() != NULL) {
            report_error(e.what());
        }

        else {
            report_error("Unknown error");
        }

    } catch (...) {
        report_error("Unknown error");
    }

    return R_NilValue; //Avoid compiler warnings
}


SEXP epiG_fit_filename(SEXP r_filename, SEXP r_ref_genom_filename, SEXP r_alt_genom_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_max_threads, SEXP r_max_chunk_size, SEXP r_config) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string ref_genom_filename = get_value<std::string>(r_ref_genom_filename);
	const std::string alt_genom_filename = get_value<std::string>(r_alt_genom_filename);
	const std::string refName = get_value<std::string>(r_refName);
	const t_position start = get_value<t_position>(r_start);
	const t_position end = get_value<t_position>(r_end);
	const t_count max_threads = get_value<t_count>(r_max_threads);
	const t_count max_chunk_size = get_value<t_count>(r_max_chunk_size);

	AlgorithmConfiguration const config(r_config);


    chunk_optimizer opt = create_base_chunk_optimizer(filename, ref_genom_filename, alt_genom_filename, refName, start, end, max_threads, max_chunk_size, config);
    //TODO remove
    //chunk_optimizer opt = create_read_chunk_optimizer(filename, ref_genom_filename, alt_genom_filename, refName, start, end, max_threads, max_chunk_size, config);


	opt.run();

    rList res;

	res.attach(rObject(opt.offset), "offset");
	res.attach(rObject(opt.number_of_chunks), "number_of_chunks");
	res.attach(rObject(opt.get_chunks_start()), "chunks_start");
	res.attach(rObject(opt.get_chunks_end()), "chunks_end");
	res.attach(rObject(opt.get_methylome()), "g");
	//res.attach(rObject(opt.compute_partial_likelihood()), "loglike");
	res.attach(rObject(opt.get_alleles()), "alleles");
	res.attach(rObject(opt.get_strands()), "strands");
	res.attach(rObject(opt.get_unique_assigned()), "unique");
	res.attach(rObject(opt.get_readIDs()), "read_id");
	res.attach(rObject(opt.get_allele_occupancy()), "occupancy");


    return rObject(res);
}

SEXP r_epiG_fit_filename(SEXP r_filename, SEXP r_ref_genom_filename, SEXP r_alt_genom_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_max_threads, SEXP r_max_chunk_size, SEXP r_config) {

	try {

		return epiG_fit_filename(r_filename, r_ref_genom_filename, r_alt_genom_filename, r_refName, r_start, r_end, r_max_threads, r_max_chunk_size, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}

SEXP epiG_fit_filename_chunks(SEXP r_filename, SEXP r_ref_genom_filename, SEXP r_alt_genom_filename, SEXP r_refName, SEXP r_chunks_start, SEXP r_chunks_end, SEXP r_max_threads, SEXP r_config) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string ref_genom_filename = get_value<std::string>(r_ref_genom_filename);
	const std::string alt_genom_filename = get_value<std::string>(r_alt_genom_filename);
	const std::string refName = get_value<std::string>(r_refName);
	const t_positions chunks_start = get_value<t_positions>(r_chunks_start);
	const t_positions chunks_end = get_value<t_positions>(r_chunks_end);
	const t_count max_threads = get_value<t_count>(r_max_threads);

	AlgorithmConfiguration config(r_config);

	chunk_optimizer opt(filename, ref_genom_filename, alt_genom_filename, refName, chunks_start, chunks_end, max_threads, config);

	opt.run();

    rList res;

	res.attach(rObject(opt.offset), "offset");
	res.attach(rObject(opt.number_of_chunks), "number_of_chunks");
	res.attach(rObject(opt.get_chunks_start()), "chunks_start");
	res.attach(rObject(opt.get_chunks_end()), "chunks_end");
	res.attach(rObject(opt.get_methylome()), "g");
	//res.attach(rObject(opt.compute_partial_likelihood()), "loglike");
	res.attach(rObject(opt.get_alleles()), "alleles");
	res.attach(rObject(opt.get_strands()), "strands");
	res.attach(rObject(opt.get_unique_assigned()), "unique");
	res.attach(rObject(opt.get_readIDs()), "read_id");
	res.attach(rObject(opt.get_allele_occupancy()), "occupancy");

    return rObject(res);
}

SEXP r_epiG_fit_filename_chunks(SEXP r_filename, SEXP r_ref_genom_filename, SEXP r_alt_genom_filename, SEXP r_refName, SEXP r_chunks_start, SEXP r_chunks_end, SEXP r_max_threads, SEXP r_config) {

	try {

		return epiG_fit_filename_chunks(r_filename, r_ref_genom_filename, r_alt_genom_filename, r_refName, r_chunks_start, r_chunks_end, r_max_threads, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}


#endif /* EPIG_R_INTERFACE_HPP_ */
