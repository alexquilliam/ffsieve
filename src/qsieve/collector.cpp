/*Collect the information needed for sieveing*/

#include <collector.hpp>
#include <chrono>

/*Public methods*/

mpz_vector Collector::collect(mpz_class n) {
	//setup
	long bound = get_bound(n);

	int_vector factor_base = get_factor_base(n, bound);

	std::cout << "\nFactor base:" << std::endl;
	std::cout << factor_base << std::endl;

	mpz_class mpz_n_sqrt(sqrt(n));
	long n_sqrt = mpz_n_sqrt.get_ui();
	long block_size = n_sqrt / 2;
	int num_blocks = 0;
	mpz_vector poly_values = generate_poly_values(n, n_sqrt + (num_blocks++ * block_size), block_size);
	poly_values.erase(poly_values.begin()); //remove negative value

	std::cout << "\nPolynomial values:" << std::endl;
	std::cout << poly_values << std::endl;

	mpz_vector bsmooth_values(0);
	std::vector<std::vector<int>> exp_vec(0, std::vector<int>(0));

	size_t i = 0;
	while(bsmooth_values.size() < factor_base.size()) {
		int_vector exp_vector = is_bsmooth(poly_values[i], factor_base);

		if(exp_vector.size() != 0) {
			bsmooth_values.push_back(poly_values[i]);
			exp_vec.push_back(exp_vector);
		}

		i++;

		if(i == poly_values.size()) {
			mpz_vector new_poly_batch = generate_poly_values(n, n_sqrt + (num_blocks++ * block_size), block_size);
			poly_values.insert(poly_values.end(), new_poly_batch.begin(), new_poly_batch.end());
		}
	}

	std::cout << "\nb-smooth values:" << std::endl;
	std::cout << bsmooth_values << std::endl;

	std::cout << "\nb-smooth values size: " << bsmooth_values.size() << ", factor base size: " << factor_base.size() << std::endl;

	exp_matrix = Matrix(exp_vec);

	std::cout << "\nInput exponent matrix: " << std::endl;
	std::cout << get_matrix() << std::endl;

	return bsmooth_values;
}

Matrix Collector::get_matrix() {
	return exp_matrix;
}

/*Private methods*/

//determine if a number is bsmooth up to a bound; return the factor exponent vector if it is smooth
//otherwise, return an empty vector of size 0
int_vector Collector::is_bsmooth(mpz_class n, int_vector smooth_primes) {
	int_vector factor_exp_vector(smooth_primes.size());

	mpz_class quotient(n);
	mpz_class remainder(0);

	for(size_t i = 0; i < smooth_primes.size(); i++) {
		mpz_tdiv_r_ui(remainder.get_mpz_t(), quotient.get_mpz_t(), smooth_primes[i]);

		if(remainder != 0) {
			factor_exp_vector[i] = 0;
			continue;
		}

		//divide out the relevant powers of the smooth primes
		while(true) {
			mpz_tdiv_r_ui(remainder.get_mpz_t(), quotient.get_mpz_t(), smooth_primes[i]);

			if(remainder != 0) {
				break;
			}else {
				mpz_divexact_ui(quotient.get_mpz_t(), quotient.get_mpz_t(), smooth_primes[i]);
				factor_exp_vector[i]++;
			}
		}
	}

	if(quotient == 1) {
		return factor_exp_vector;
	}else {
		return int_vector(0);
	}
}

//generate the values of the polynomial x^2 - n over a range of x
mpz_vector Collector::generate_poly_values(mpz_class n, long poly_start, long block_size) {
	long poly_end = poly_start + block_size;

	mpz_vector poly_values(block_size);
	for(long x = poly_start; x < poly_end; x++) {
		poly_values[x - poly_start] = x * x;
		poly_values[x - poly_start] -= n;
	}

	return poly_values;
}

//get all of the primes less than b (bound) where n/p (the legendre symbol) is 1
int_vector Collector::get_factor_base(mpz_class n, long bound) {
	int_vector primes = get_primes((int) bound);

	int_vector factor_base(0);

	//calculate the legendre symbol n/p, where p is a prime
	for(size_t i = 0; i < primes.size(); i++) {
		mpz_class mpzprime(primes[i]);
		if(mpz_legendre(n.get_mpz_t(), mpzprime.get_mpz_t()) == 1) {
			factor_base.push_back(primes[i]);
		}
	}

	//the legendre symbol is not defined at 2, but all numbers have a quadratic residue mod 2
	if(factor_base[0] != 2) {
		factor_base.insert(factor_base.begin(), 2);
	}

	return factor_base;
}

//claculate the upper bound of the factor base with the formula
//ceil(e^((.5 + o(1))sqrt(ln(N)ln(ln(N))))). Here, o(1) will equal .22.
//it can also be calculated by 2*e^(sqrt(ln(n)ln(ln(n))/4))
long Collector::get_bound(mpz_class n) {
	long double exponent = 0;
	double o = .22;

	exponent = (.5 + o) * sqrt(ln(n) * ln(ln(n)));

	return ceil(exp(exponent));
}
