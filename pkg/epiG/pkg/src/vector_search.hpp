#ifndef VECTOR_SEARCH_HPP_
#define VECTOR_SEARCH_HPP_

extern "C" {
SEXP r_epiG_locate(SEXP pattern, SEXP text);
}


arma::uvec locate(arma::uvec const& pattern, arma::uvec const& text) {

	arma::uvec pos;

	arma::u32 p = 0;
	arma::uvec::const_iterator pi = pattern.begin();
	for(arma::uvec::const_iterator ti = text.begin(); ti != text.end(); ++ti, ++p) {

		if (*pi == *ti) {

			++pi;

			if(pi == pattern.end()) {
				pos.resize(pos.n_elem+1);
				pos(pos.n_elem - 1) = p - pattern.n_elem + 1;
				pi = pattern.begin();
			}
		}

		else {
			pi = pattern.begin();
		}

	}

	return pos;

}

SEXP epiG_locate(SEXP r_pattern, SEXP r_text) {

	const arma::uvec pattern = get_value<arma::uvec>(r_pattern);
	const arma::uvec text = get_value<arma::uvec>(r_text);

	arma::uvec pos = locate(pattern, text);

    return rObject(pos);
}

SEXP r_epiG_locate(SEXP r_pattern, SEXP r_text) {

    try {

        return epiG_locate(r_pattern, r_text);

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


#endif /* VECTOR_SEARCH_HPP_ */
