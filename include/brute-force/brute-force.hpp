/*Use brute force to find factors of small numbers*/

#pragma once

//local includes
#include <utils.hpp>

class BruteForce {
	public:
		 std::vector<std::vector<int>> get_factors(mpz_class num);
};
