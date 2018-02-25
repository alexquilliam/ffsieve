/*A collection of useful functions and other random things*/

#include <utils.hpp>

int_vector get_primes(int n) {
	int sieve[n];

	for(int i = 0; i < n; i++) {
		sieve[i] = 1;
	}

	sieve[0] = 0;
	sieve[1] = 0;

	int num_primes = 0;
	for(int i = 2; i < n; i++) {
		if(sieve[i] == 1) {
			for(int j = i; i * j < n; j++) {
				sieve[i * j] = 0;
			}
		}

		if(sieve[i] == 1) {
			num_primes++;
		}
	}

	int_vector primes(num_primes);

	for(int i = 0, j = 0; i < n; i++) {
		if(sieve[i] == 1) {
			primes[j++] = i;
		}
	}

	return primes;
}
