/*Factor a number with the quadratic sieve*/

#pragma once

//local includes
#include <utils.hpp>
#include <collector.hpp>
#include <processor.hpp>

class QuadraticSieve {
	public:
		mpz_vector get_factors(mpz_class value);
};
