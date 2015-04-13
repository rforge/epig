/*
 * math_tools.hpp
 *
 *  Created on: Jan 26, 2014
 *      Author: martin
 */

#ifndef MATH_TOOLS_HPP_
#define MATH_TOOLS_HPP_

#include <armadillo>
using namespace arma;

template<unsigned int a>
u32 pow(u32 n) {

	u32 r = 1;

	for(u32 i = 0; i < n; ++i) {
		r *=a;
	}

	return r;
}

template<typename T>
T max(T a, T b) {
	return a > b ? a : b;
}

template<typename T>
T min(T a, T b) {
	return a < b ? a : b;
}

template<typename T>
T square(T x) {
	return x*x;
}

template<typename T, unsigned int base>
void toBase(Col<T> & base_int, u32 decimal_int) {

	int c = decimal_int;

	for(int i = 0; i < base_int.n_elem; ++i) {
		base_int(i) = static_cast<T>(c % base);
		c = c / base;
	}
}

template<typename T>
arma::u32 argmin(arma::Col<T> const& x) {

	T xmin = x(0);
	arma::u32 j = 0;
	for(arma::u32 i = 1; i < x.n_elem; ++i) {
		if(x(i) < xmin) {
			xmin = x(i);
			j = i;
		}
	}

	return j;
}

template<typename T>
arma::u32 argmax(arma::Col<T> const& x) {

	T xmax = x(0);
	arma::u32 j = 0;
	for(arma::u32 i = 1; i < x.n_elem; ++i) {
		if(x(i) > xmax) {
			xmax = x(i);
			j = i;
		}
	}

	return j;
}

// Compute binomial distribution with p = 0.5 table
// Row - total trials
// Col - success
arma::Mat<double> make_log_binom_half_table(unsigned int const N_max) {

	arma::Col<double> tmp(N_max);
	for(unsigned int n = 0; n < N_max; ++n) {
		tmp(n) = log(n+1);
	}

	arma::Mat<double> table(N_max+1, N_max+1);

	for(unsigned int n = 0; n <= N_max; ++n) {

		table.row(n).fill(n * log(0.5));

		for(unsigned int k = 1; k < n; ++k) {
			table(n,k) += sum(tmp.subvec(k,n-1))-sum(tmp.subvec(0,n-k-1));
		}
	}

	return table;
}


arma::uvec factor(arma::uvec const& x) {

	arma::uvec levels = unique(x);

	arma::uvec r(x.n_elem);

	for(arma::uword i = 0; i < x.n_elem; ++i) {
		arma::uvec level = find(levels == x(i), 1); //this will be a vector with 1 element
		r(i) = level(0);
	}

	return r;
}

class base_counter {

	arma::uvec number;

public:

	unsigned int const base;
	unsigned int const n_dec; //number of decimals

	base_counter(unsigned int const base, unsigned int n_dec) : number(n_dec), base(base), n_dec(n_dec) {}

	void set_zero() {
		number.zeros();
	}

	base_counter const& operator++() {

		for(arma::uword i = 0; i < n_dec; ++i) {
			number(i) = (number(i) + 1) % base;

			if(number(i) != 0) {
				break;
			}
		}

		return *this;
	}

	unsigned int operator()(unsigned int i) const {
		return number(i);
	}

	bool is_zero() const {
		return accu(number != 0) == 0;
	}

	arma::uvec const get_number() const {
		return number;
	}

};

#endif /* MATH_TOOLS_HPP_ */
