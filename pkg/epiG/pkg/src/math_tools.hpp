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

#endif /* MATH_TOOLS_HPP_ */
