/*Collect the information needed for sieveing*/

#pragma once

//local includes
#include <utils.hpp>
#include <matrix.hpp>

//std includes
#include <fstream>
#include <string>
#include <unistd.h>

class Collector {
	private:
		Matrix exp_matrix;

		int_vector is_bsmooth(mpz_class n, int_vector smooth_primes);
		mpz_vector generate_poly_values(mpz_class n, mpz_class poly_start, long block_size);
		long get_bound(mpz_class n);
		int_vector get_factor_base(mpz_class n, long bound);
		int get_cache_size(int cache_index);
		int_vector msqrt(mpz_class a, int p);
		int_vector msqrt_pow(mpz_class n, int p, int e);
		mpz_vector sieve(mpz_class n, int_vector factor_base, int block_size);
		std::vector<std::vector<int>> generate_solutions(mpz_class n, int_vector factor_base, int power);
	public:
		mpz_vector collect(mpz_class n);
		Matrix get_matrix();
};

/*#pragma once

//local includes
#include <utils.hpp>
#include <matrix.hpp>

//std includes
#include <fstream>
#include <string>
#include <bitset>

class Collector {
	private:
		Matrix exp_matrix;

		int_vector is_bsmooth(mpz_class n, int_vector smooth_primes);
		mpz_vector generate_poly_values(mpz_class n, mpz_class poly_start, long block_size);
		long get_bound(mpz_class n);
		int_vector get_factor_base(mpz_class n, long bound);
		int get_cache_size(int cache_index);
		int_vector msqrt(mpz_class a, int p);
		int_vector msqrt_pow(mpz_class n, int p, int e);
		mpz_vector sieve(mpz_class n, int_vector factor_base, int block_size);
		int_vector generate_solutions(mpz_class n, int_vector factor_base);
		void filter_smooth_values(mpz_class n, mpz_vector poly_chunk, std::vector<double> sieve, int_vector factor_base, mpz_vector *bsmooth_values, std::vector<std::vector<int>> *exp_mat, double cutoff);
	public:
		mpz_vector collect(mpz_class n);
		Matrix get_matrix();
};
*/
