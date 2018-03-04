/*Process the data gathered in the collection phase*/

#include <processor.hpp>

mpz_vector Processor::process(mpz_class num, mpz_vector bsmooth_values, Matrix mat) {
	mat.mod();
	mat.transpose();
	mat.to_rref();

	int_vector null_space(mat.size());
	null_space = mat.get_null_space();

	#if DEBUG
		std::cout << "Final RREF matrix: " << std::endl;
		std::cout << mat << std::endl;

		std::cout << "Null space: " << std::endl;
		std::cout << null_space << std::endl;
	#endif

	mpz_vector multiples = get_multiples(num, bsmooth_values, null_space);

	#if DEBUG
		std::cout << "\nMultiples: " << std::endl;
		std::cout << multiples[0] << " and " << multiples[1] << std::endl;
	#endif

	mpz_vector factors(2);
	factors[0] = gcd(num, multiples[0]);
	factors[1] = gcd(num, multiples[1]);

	#if DEBUG
		std::cout << "\nFactors before correction: " << std::endl;
		std::cout << factors[0] << " and " << factors[1] << std::endl;
	#endif

	//try dividing one of the factors by 2
	if(!check(num, factors)) {
		factors[1] /= 2;
	}

	#if DEBUG
		std::cout << "\nFactors: " << std::endl;
		std::cout << factors[0] << " and " << factors[1] << std::endl;
	#endif

	return factors;
}

mpz_vector Processor::get_multiples(mpz_class num, mpz_vector bsmooth_values, int_vector null_space) {
	mpz_vector square_factors(0);
	for(size_t i = 0; i < bsmooth_values.size(); i++) {
		if(null_space[i] == 1) {
			square_factors.push_back(bsmooth_values[i]);
		}
	}

	mpz_vector originals(0);
	for(size_t i = 0; i < square_factors.size(); i++) {
		mpz_class original(square_factors[i]);
		original += num;
		original = sqrt(original);

		originals.push_back(original);
	}

	mpz_class b_square(1);
	for(mpz_class m : square_factors) {
		b_square *= m;
	}

	b_square = sqrt(b_square);

	mpz_class o_square(1);
	for(mpz_class m : originals) {
		o_square *= m;
	}

	mpz_vector multiples(2);
	multiples[0] = o_square - b_square;
	multiples[1] = o_square + b_square;

	return multiples;
}
