/*Finds the prime factors of a number.*/

#include <decomposer.hpp>

/*Public Functions*/

int_vector Decomposer::decompose(int num) {
	int_vector prime_factors(0);

	int prime = 0;
	int composite = num;

	while(prime != 1) {
		prime = lowest_prime_factor(composite);
		prime_factors.push_back(prime);

		composite /= prime;
	}

	prime_factors[prime_factors.size() - 1] = composite;

	return prime_factors;
}

/*Private Functions*/

int Decomposer::lowest_prime_factor(int num) {
	int lpf = 1;

	int_vector primes = get_primes(num);

	for(unsigned int i = 0; i < primes.size(); i++) {
		if(num % primes[i] == 0) {
			lpf = primes[i];
			return lpf;
		}
	}

	return lpf;
}
