//local includes
#include <utils.hpp>
#include <decomposer.hpp>
#include <brute-force.hpp>
#include <qsieve.hpp>

int main(int argc, char **argv) {
	mpz_class num(argv[1]);

	//don't bother factoring if the number is prime
	int is_prime = mpz_probab_prime_p(num.get_mpz_t(), 35);
	if(is_prime == 2 || is_prime == 1) {
		std::cout << "Number is prime. Its factors are 1 and " << num << "." << std::endl;
		return 0;
	}

	//just use trial division if the number is small
	if(num < 75000) {
		BruteForce bf;
		std::vector<std::vector<int>> factors = bf.get_factors(num);

		std::cout << "The factors of " << num << " are: " << std::endl << factors << std::endl;

		return 0;
	}

	QuadraticSieve qs;
	qs.get_factors(num);

	return 0;
}
