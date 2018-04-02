//local includes
#include <utils.hpp>
#include <brute-force.hpp>
#include <pollards-rho.hpp>
#include <qsieve.hpp>

mpz_class generate_semiprime(mpz_class seed) {
	mpz_class semiprime_one(seed);
	mpz_class sp_one_stop(semiprime_one + 1000);
	for(; semiprime_one < sp_one_stop; semiprime_one += 2) {
		if(mpz_probab_prime_p(semiprime_one.get_mpz_t(), 50) == 0) {
			continue;
		}else {
			break;
		}
	}

	mpz_class semiprime_two(semiprime_one + 2);
	mpz_class sp_two_stop(semiprime_two + 1000);
	for(; semiprime_two < sp_two_stop; semiprime_two += 2) {
		if(mpz_probab_prime_p(semiprime_two.get_mpz_t(), 50) == 0) {
			continue;
		}else {
			break;
		}
	}

	return semiprime_one * semiprime_two;
}

int main(int argc, char **argv) {
	mpz_class target_num = generate_semiprime(mpz_class("2000000000001"));

	if(mpz_probab_prime_p(target_num.get_mpz_t(), 50) == 2 || mpz_probab_prime_p(target_num.get_mpz_t(), 50) == 1) {
		std::cout << target_num << " is prime; its factors are 1 and " << target_num << "." << std::endl;
		exit(0);
	}

	std::cout << "\nFactoring " << target_num << std::endl;
	std::cout << "The target number is " << mpz_sizeinbase(target_num.get_mpz_t(), 10) << " digits long, and is " << mpz_sizeinbase(target_num.get_mpz_t(), 2) << " bits." << std::endl << std::endl;

	mpz_vector factors;
	long start = get_current_time();
	factors = pollards_rho(target_num, mpz_class(100000));
	long elapsed = get_current_time() - start;

	if(factors[0] != -1) {
		std::cout << "The factors of " << target_num << " are " << factors << ", found in " << elapsed << "ms with Pollard's rho algorithm." << std::endl;
		exit(0);
	}else {
		std::cout << target_num << " was not factorable with Pollard's rho; trying the quadratic sieve." << std::endl;
	}

	start = get_current_time();
	QuadraticSieve qs;
	qs.get_factors(target_num);
	elapsed = get_current_time() - start;

	std::cout << "The factors of " << target_num << " are " << factors << ", found in " << elapsed << "ms with the quadratic sieve." << std::endl;

	return 0;
}
























/*//local includes
#include <utils.hpp>
#include <decomposer.hpp>
#include <brute-force.hpp>
#include <qsieve.hpp>

int main(int argc, char **argv) {
	std::cout.setstate(std::ios_base::failbit);

	mpz_class value(90200);
	mpz_class end_value(90300);
	int total = 0, sucessful = 0, unsucessful = 0, trivial = 0, prime = 0, loop = 0, other = 0;
	for(; value < end_value; value++) {
		if(mpz_probab_prime_p(value.get_mpz_t(), 50) == 2 || mpz_probab_prime_p(value.get_mpz_t(), 50) == 1) {
			prime++;
			total++;
			continue;
		}

		QuadraticSieve qs;
		std::future<mpz_vector>* future = new std::future<mpz_vector>();
		*future = std::async(&QuadraticSieve::get_factors, &qs, value);

		long start = get_current_time();
		bool fail = false;
		while((*future).wait_for(std::chrono::milliseconds(0)) != std::future_status::ready) {
			if(get_current_time() - start > 500) {
				fail = true;
				break;
			}
		}

		mpz_vector factors;

		if(fail) {
			loop++;
			total++;
			unsucessful++;
			continue;
		}else {
			factors = (*future).get();
		}

		if(std::find(factors.begin(), factors.end(), 1) != factors.end() && std::find(factors.begin(), factors.end(), value) != factors.end()) {
			trivial++;
			total++;
			unsucessful++;
			continue;
		}else if(!check(value, factors)) {
			other++;
			total++;
			unsucessful++;
			continue;
		}else {
			sucessful++;
			total++;
			continue;
		}
	}

	std::cout.clear();
	std::cout << "Out of " << total << " factorizations, " << sucessful << " succeeded, " << unsucessful << " failed, and " << prime << " were prime." << std::endl;
	std::cout << "Of the " << unsucessful << " factorizations that failed, " << trivial << " were trivial, " << loop << " caused an infinite loop, and " << other << " failed completely." << std::endl;
	std::cout.setstate(std::ios_base::failbit);

	return 0;
}*/
