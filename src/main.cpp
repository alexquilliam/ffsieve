//local includes
#include <utils.hpp>
#include <decomposer.hpp>
#include <brute-force.hpp>
#include <qsieve.hpp>

int main(int argc, char **argv) {
	/*mpz_class num(argv[1]);

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
	}*/

	/*mpz_class coprime_one(90281);
	for(mpz_class coprime_two(90281); coprime_two < 90971; coprime_two += 2) {
		if(!mpz_probab_prime_p(coprime_two.get_mpz_t(), 50)) continue;

		QuadraticSieve qs;
		qs.get_factors(coprime_one * coprime_two);
	}*/

	mpz_class coprime_one("100001");
	mpz_class cp_one_stop(coprime_one + 1000);
	for(; coprime_one < cp_one_stop; coprime_one += 2) {
		if(mpz_probab_prime_p(coprime_one.get_mpz_t(), 50) != 2) continue;
		break;
	}

	mpz_class coprime_two(coprime_one + 2);
	mpz_class cp_two_stop(coprime_two + 1000);
	for(; coprime_two < cp_two_stop; coprime_two += 2) {
		if(mpz_probab_prime_p(coprime_two.get_mpz_t(), 50) != 2) continue;
		break;
	}

	std::cout << "Factoring " << coprime_one * coprime_two << " (" << coprime_one << " * " << coprime_two << ")." << std::endl;

	QuadraticSieve qs;
	qs.get_factors(coprime_one * coprime_two);

	return 0;
}
