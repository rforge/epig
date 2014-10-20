/* Routines for TODO
 Intended for use with R.
 Copyright (C) 2012 Martin Vincent

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef EPIG_ALGORITHM_CONFIG_H_
#define EPIG_ALGORITHM_CONFIG_H_

template<typename type>
static type getConfigAttribute(rList const& config, std::string const& name) {

	int index;
	if (index = config.getIndex(name), index >= 0) {

		return get_value<type>(config.get(index));

	} else {

		std::string msg = "Missing configuration parameter : ";
		throw std::domain_error(msg.append(name).c_str());
		return type(); //avoid compiler warnings
	}
}

//TODO do we need this
template<typename type>
static field<type> getConfigList(rList const& config, std::string const& name) {

	int index;
	if (index = config.getIndex(name), index >= 0) {

		return get_field<type>(config.get(index));

	} else {

		std::string msg = "Missing configuration parameter : ";
		throw std::domain_error(msg.append(name).c_str());
		return field<type>(); //avoid compiler warnings
	}
}

class AlgorithmConfiguration {

public:

	const t_count max_iterations; //max number of iterations per optim loop

	const t_models fwd_model;
	const t_models rev_model;

	const double sequence_quality_adjust;

    const t_count reads_hard_limit;

    const double haplo_prior;
    const double ref_prior;

    const field<arma::vec> genotype_priors;

    const t_position min_overlap_length;

	bool const verbose;

	AlgorithmConfiguration(rList const& config) :

			max_iterations(getConfigAttribute<t_count>(config, "max_iterations")),

			fwd_model(getConfigList<t_model>(config, "fwd_model")),

			rev_model(getConfigList<t_model>(config, "rev_model")),

			sequence_quality_adjust(getConfigAttribute<double>(config, "sequence_quality_adjust")),

            reads_hard_limit(getConfigAttribute<t_count>(config, "reads_hard_limit")),

            haplo_prior(getConfigAttribute<double>(config, "haplo_prior")),

            ref_prior(getConfigAttribute<double>(config, "ref_prior")),

            genotype_priors(getConfigList<arma::vec>(config, "log_prior")),

            min_overlap_length(getConfigAttribute<t_position>(config, "min_overlap_length")),

			verbose(getConfigAttribute<bool>(config, "verbose")) {
	}
};

#endif /* EPIG_ALGORITHM_CONFIG_H_ */
