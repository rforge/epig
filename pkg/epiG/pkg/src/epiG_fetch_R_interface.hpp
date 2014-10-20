#ifndef EPIG_FETCH_R_INTERFACE_HPP_
#define EPIG_FETCH_R_INTERFACE_HPP_

extern "C" {

SEXP r_epiG_fetch_reads(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end);
SEXP r_epiG_fetch_reads_info(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end);

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

#endif /* EPIG_FETCH_R_INTERFACE_HPP_ */
