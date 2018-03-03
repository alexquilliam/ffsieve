/*Use brute force to find factors of small numbers*/

#include <brute-force.hpp>

std::vector<std::vector<int>> BruteForce::get_factors(mpz_class num) {
	long n = num.get_ui();

	std::vector<std::vector<int>> factors(0, std::vector<int>(0));

	for(int i = 1; i <= ceil(sqrt(n)); i++) {
		if(n % i == 0) {
			std::vector<int> vec(2);
			vec[0] = i;
			vec[1] = n / i;

			factors.push_back(vec);
		}
	}

	return factors;
}
