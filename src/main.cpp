//local includes
#include <utils.hpp>
#include <decomposer.hpp>
#include <brute-force.hpp>
#include <qsieve.hpp>

int main(int argc, char **argv) {
	mpz_class semiprime_one("200000000001");
	mpz_class sp_one_stop(semiprime_one + 1000);
	for(; semiprime_one < sp_one_stop; semiprime_one += 2) {
		if(mpz_probab_prime_p(semiprime_one.get_mpz_t(), 50) != 2) continue;
		break;
	}

	mpz_class semiprime_two(semiprime_one + 2);
	mpz_class sp_two_stop(semiprime_two + 1000);
	for(; semiprime_two < sp_two_stop; semiprime_two += 2) {
		if(mpz_probab_prime_p(semiprime_two.get_mpz_t(), 50) != 2) continue;
		break;
	}

	mpz_class target_num(semiprime_one * semiprime_two);

	std::cout << "\nFactoring " << target_num << " (" << semiprime_one << " * " << semiprime_two << ")." << std::endl;
	std::cout << "The target number is " << mpz_sizeinbase(target_num.get_mpz_t(), 10) << " digits long, and is " << mpz_sizeinbase(target_num.get_mpz_t(), 2) << " bits." << std::endl;

	long start = get_current_time();
	QuadraticSieve qs;
	qs.get_factors(semiprime_one * semiprime_two);
	long elapsed = get_current_time() - start;

	std::cout << "The quadratic sieve took " << elapsed << "ms to factor " << target_num << "." << std::endl;

	return 0;
}
