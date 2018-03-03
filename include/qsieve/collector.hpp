/*Collect the information needed for sieveing*/

#pragma once

//local includes
#include <utils.hpp>
#include <matrix.hpp>

class Collector {
	private:
		Matrix exp_matrix;

		int_vector is_bsmooth(mpz_class n, int_vector smooth_primes);
		mpz_vector generate_poly_values(mpz_class n, long poly_start, long block_size);
		long get_bound(mpz_class n);
		int_vector get_factor_base(mpz_class n, long bound);

	public:
		mpz_vector collect(mpz_class n);
		Matrix get_matrix();
};
