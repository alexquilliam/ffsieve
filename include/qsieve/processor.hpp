/*Process the data gathered in the collection phase*/

#pragma once

//local includes
#include <utils.hpp>
#include <matrix.hpp>

class Processor {
	private:
		mpz_vector get_multiples(mpz_class num, mpz_vector bsmooth_values, int_vector null_space);
	public:
		mpz_vector process(mpz_class num, mpz_vector bsmooth_values, Matrix mat);
};
