/*Finds the prime factors of a number.*/

#pragma once

//local includes
#include <utils.hpp>

class Decomposer {
	private:
		int lowest_prime_factor(int num);

	public:
		int_vector decompose(int num);
};
