/*Factor a number with the quadratic sieve*/

#include <qsieve.hpp>

mpz_vector QuadraticSieve::get_factors(mpz_class value) {
	mpz_vector factors(0);

	Collector collector;
	mpz_vector bsmooth_values = collector.collect(value);

	Processor processor;
	factors = processor.process(value, bsmooth_values, collector.get_matrix());

	return factors;
}
