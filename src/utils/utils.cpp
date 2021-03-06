/*A collection of useful functions and other random things*/

#include <utils.hpp>

//get the current time in milliseconds
long get_current_time() {
	auto current_time = std::chrono::system_clock::now();
	auto current_time_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(current_time);
	auto epoch = current_time_ms.time_since_epoch();
	long current_time_long = std::chrono::duration_cast<std::chrono::milliseconds>(epoch).count();

	return current_time_long;
}

//check if the factors of a number are actually factors
bool check(mpz_class n, mpz_vector factors) {
	mpz_class product(1);
	for(mpz_class m : factors) {
		product *= m;
	}

	if(n == product) {
		return true;
	}else {
		return false;
	}
}

//get all primes less than n with the Sieve of Eratosthenes
int_vector get_primes(int n) {
	int *sieve = new int[n]; //don't want to overflow the stack

	for(int i = 0; i < n; i++) {
		sieve[i] = 1;
	}

	sieve[0] = 0;
	sieve[1] = 0;

	int num_primes = 0;
	for(int i = 2; i < n; i++) {
		if(sieve[i] == 1) {
			for(long j = i; i * j < n; j++) {
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

	delete[] sieve;

	return primes;
}

//calculate the (approximate) natural logarithim (ln) of a number
double ln(mpz_class n) {
	if(n.fits_uint_p()) {
		return log((double) n.get_ui());
	}

	int n_size = mpz_sizeinbase(n.get_mpz_t(), 10);

	int ppart_size = 5;
	int upart_size = n_size - ppart_size;

	mpz_class ppart(0); //precise part (the first 5 digits)
	mpz_class upart(0); //unprecise part (everything else)

	//10^(n_size - 5)
	mpz_ui_pow_ui(upart.get_mpz_t(), 10, upart_size);

	ppart = n / upart;

	double ppart_log = log(ppart.get_ui());
	double upart_log = upart_size * log(10);

	return ppart_log + upart_log;
}

//an implementation of the modulus operator that works properly on negative numerators
long long modn(long long a, long long b) {
	long long rem = a % b;

	if(rem < 0) {
		rem += b;
	}

	return rem;
}
