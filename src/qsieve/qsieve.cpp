/*Factor a number with the quadratic sieve*/

#include <qsieve.hpp>

mpz_vector QuadraticSieve::get_factors(mpz_class value) {
	mpz_vector factors(0);

	Collector collector;
	mpz_vector bsmooth_values = collector.collect(value);

	Processor processor;
	factors = processor.process(value, bsmooth_values, collector.get_matrix());

	if(factors[0] == 1 || factors[1] == 1) {
		std::cout << "\nThe factorization of " << value << " is correct, but they are trivial factors (1 and " << value << "." << std::endl;
		return factors;
	}

	bool is_correct = check(value, factors);

	if(is_correct) {
		std::cout << "\nThe factorization of " << value << " (" << factors[0] << " * " << factors[1] << ") is correct." << std::endl;
	}else {
		std::cout << "\nThe factorization of " << value << " (" << factors[0] << " * " << factors[1] << ") is incorrect." << std::endl;
	}

	return factors;
}
