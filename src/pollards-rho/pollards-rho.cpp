#include <pollards-rho.hpp>

//modified version of the C++ code from Wikipedia
mpz_vector pollards_rho(mpz_class n, mpz_class cutoff) {
	mpz_class x(2);
	mpz_class y(2);
	mpz_class factor(1);

	mpz_class cutoff_count(0);
	while(factor == 1 && cutoff_count <= cutoff) {
		x = ((x * x) + 1) % n;
		y = ((y * y) + 1) % n;
		y = ((y * y) + 1) % n;
		mpz_abs(factor.get_mpz_t(), mpz_class(x - y).get_mpz_t());
		mpz_gcd(factor.get_mpz_t(), factor.get_mpz_t(), n.get_mpz_t());

		cutoff_count++;
	}

	mpz_vector factors(2);
	if(cutoff_count == cutoff || factor == 1 || factor == n) {
		factors[0] = -1;
		factors[1] = -1;

		return factors;
	}

	factors[0] = factor;
	factors[1] = n / factor;
	return factors;
}
