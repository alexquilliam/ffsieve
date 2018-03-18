#include <collector.hpp>

mpz_vector Collector::collect(mpz_class n) {
	long bound = get_bound(n);

	#if VERBOSE
		std::cout << "\nGot smoothness bound." << std::endl;
	#endif

	int_vector factor_base = get_factor_base(n, bound);

	#if VERBOSE
		std::cout << "Got factor base." << std::endl;
	#endif

	mpz_class n_sqrt(sqrt(n));
	long block_size = 100;//get_cache_size(1) / sizeof(mpz_class);

	#if VERBOSE
		std::cout << "Got polynomial block size." << std::endl;
	#endif

	mpz_vector bsmooth_values = sieve(n, factor_base, block_size);

	#if VERBOSE
		std::cout << "Got bsmooth values." << std::endl;
	#endif

	return bsmooth_values;
}

//sieve, returning the bsmooth values
mpz_vector Collector::sieve(mpz_class n, int_vector factor_base, int block_size) {
	const size_t factor_base_size = factor_base.size();

	mpz_vector bsmooth_values(0);
	std::vector<std::vector<int>> exp_mat(0, std::vector<int>(0));

	mpz_class n_sqrt(sqrt(n));

	//generate the first batch of polynomials
	int num_blocks = 0;
	mpz_vector poly_chunk = generate_poly_values(n, n_sqrt + (num_blocks++ * block_size) + 1, block_size);

	//get the start location of the sieve (the modular square roots) for each smooth factor
	int_vector sieve_start(factor_base_size * 2);
	long long ceil_n_sqrt(n_sqrt.get_ui() + 1);
	for(size_t i = 0, j = 0; i < factor_base_size; i++, j += 2) {
		int_vector solutions = msqrt(n, factor_base[i]);

		if(solutions[0] == solutions[1]) {
			sieve_start[j] = modn(solutions[0] - ceil_n_sqrt, factor_base[i]);
			sieve_start[j + 1] = -1;
		}else {
			sieve_start[j] = modn(solutions[0] - ceil_n_sqrt, factor_base[i]);
			sieve_start[j + 1] = modn(solutions[1] - ceil_n_sqrt, factor_base[i]);
		}
	}

	std::vector<double> factor_base_log(factor_base_size);
	for(size_t i = 0; i < factor_base_size; i++) {
		factor_base_log[i] = log(factor_base[i]);
	}

	//begin sieving until a sufficent number of bsmooth values have been found
	while(bsmooth_values.size() < factor_base_size) {
		int_vector sieve(poly_chunk.size());

		//start sieving
		for(size_t i = 0, j = 0; i < factor_base_size; i++, j += 2) {
			for(size_t k = j; k < j + 2; k++) {
				if(sieve_start[k] == -1) {
					continue;
				}

				size_t l;
				for(l = sieve_start[k]; l < poly_chunk.size(); l += factor_base[i]) {
					sieve[l] += factor_base_log[i];
				}

				sieve_start[k] = (l + factor_base[i]) % block_size;
			}
		}

		double x = .2;
		double cutoff = log((2 * x) * n_sqrt.get_ui());
		for(size_t i = 0; i < poly_chunk.size(); i++) {
			if(sieve[i] >= cutoff) {
				int_vector exp_vec = is_bsmooth(poly_chunk[i], factor_base);
				if(exp_vec.size() != 0) {
					bsmooth_values.push_back(poly_chunk[i]);
					exp_mat.push_back(exp_vec);

					if(bsmooth_values.size() >= factor_base_size) {
						exp_matrix = Matrix(exp_mat);
						return bsmooth_values;
					}
				}
			}
		}

		poly_chunk = generate_poly_values(n, n_sqrt + (num_blocks++ * block_size) + 1, block_size);
	}

	exp_matrix = Matrix(exp_mat);
	std::cout << exp_matrix << std::endl;

	return bsmooth_values;
}

//find the modular square root (the solution to r^2 === a (mod p)) via brute force
int_vector Collector::msqrt(mpz_class a, int p) {
	int_vector sqrts(2);

	mpz_class mpz_p(p);
	mpz_class mod_a(a % mpz_p);
	for(mpz_class r(0); r < a; r++) {
		if((r * r) % mpz_p == mod_a) {
			sqrts[0] = r.get_ui();
			sqrts[1] = p - r.get_ui();
			return sqrts;
		}
	}

	return sqrts;
}

Matrix Collector::get_matrix() {
	return exp_matrix;
}

//determine if a number is bsmooth up to a bound; return the factor exponent vector if it is smooth
//otherwise, return an empty vector of size 0
int_vector Collector::is_bsmooth(mpz_class n, int_vector smooth_primes) {
	int_vector factor_exp_vector(smooth_primes.size());

	if(n.fits_ulong_p()) {
		unsigned long quotient = n.get_ui();
		unsigned long remainder = 0;

		for(size_t i = 0, smooth_prime_size = smooth_primes.size(); i < smooth_prime_size; i++) {
			ldiv_t divide_results = ldiv(quotient, smooth_primes[i]);
			remainder = divide_results.rem;

			if(remainder != 0) {
				factor_exp_vector[i] = 0;
				continue;
			}

			//divide out the relevant powers of the smooth primes
			while(true) {
				ldiv_t power_divide_results = ldiv(quotient, smooth_primes[i]);
				remainder = power_divide_results.rem;

				if(remainder != 0) {
					break;
				}else {
					quotient = power_divide_results.quot;
					factor_exp_vector[i]++;
				}
			}
		}

		if(quotient == 1) {
			return factor_exp_vector;
		}else {
			return int_vector(0);
		}
	}else {
		mpz_class quotient(n);
		mpz_class remainder(0);

		for(size_t i = 0, smooth_prime_size = smooth_primes.size(); i < smooth_prime_size; i++) {
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
}

//generate the values of the polynomial x^2 - n over a range of x
mpz_vector Collector::generate_poly_values(mpz_class n, mpz_class poly_start, long block_size) {
	mpz_vector poly_values(block_size);
	mpz_class poly_end(poly_start + block_size);

	//avoid excessive reallocation
	mpz_class x_size(poly_start + block_size);

	mpz_t x;
	mpz_init2(x, mpz_sizeinbase(x_size.get_mpz_t(), 2));
	mpz_set(x, poly_start.get_mpz_t());

	int i = 0;
	while(i < block_size) {
		mpz_mul(poly_values[i].get_mpz_t(), x, x);
		poly_values[i] -= n;

		i++;
		mpz_add_ui(x, x, 1);
	}

	return poly_values;
}

//get all of the primes less than b (bound) where n/p (the legendre symbol) is 1
int_vector Collector::get_factor_base(mpz_class n, long bound) {
	int_vector primes = get_primes((int) bound);

	int_vector factor_base(0);

	//calculate the legendre symbol n/p, where p is a prime
	for(size_t i = 0, primes_size = primes.size(); i < primes_size; i++) {
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

//get the cache size, in bytes, of the L1, L2, or L3 caches, according to cache_index
//if it returns -1, something went wrong
int Collector::get_cache_size(int cache_index) {
	std::ifstream file;

	if(cache_index == 1) {
		file.open("/sys/devices/system/cpu/cpu0/cache/index1/size");
	}else if(cache_index == 2) {
		file.open("/sys/devices/system/cpu/cpu0/cache/index2/size");
	}else if(cache_index == 3) {
		file.open("/sys/devices/system/cpu/cpu0/cache/index3/size");
	}else {
		return -1;
	}

	if(!file) {
		return -1;
	}

	char data[20];
	file >> data;

	file.close();

	for(int i = 0; i < 20; i++) {
		if(data[i] == '\0') {
			data[i - 1] = '\0';
			break;
		}
	}

	int size = std::stoi(data);

	return size * 1000;
}
