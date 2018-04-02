#include <collector.hpp>

size_t BSMOOTH_VALUE_CUTOFF = 0;

mpz_vector Collector::collect(mpz_class n) {
	long bound = get_bound(n);

	#if VERBOSE
		std::cout << "\nGot smoothness bound." << std::endl;
	#endif

	int_vector factor_base = get_factor_base(n, bound);

	BSMOOTH_VALUE_CUTOFF = factor_base.size() + 1;

	#if VERBOSE
		std::cout << "Got factor base." << std::endl;
	#endif

	mpz_class n_sqrt(sqrt(n));
	long block_size = 1000; //100000

	#if VERBOSE
		std::cout << "Got polynomial block size." << std::endl;
	#endif

	mpz_vector bsmooth_values = sieve(n, factor_base, block_size);

	#if VERBOSE
		std::cout << "Got bsmooth values." << std::endl;
	#endif

	return bsmooth_values;
}

int canidates = 0, confirmed = 0, total = 0;
double sieve_average = 0;

//sieve, returning the bsmooth values
mpz_vector Collector::sieve(mpz_class n, int_vector factor_base, int block_size) {
	const size_t factor_base_size = factor_base.size();

	mpz_vector bsmooth_values(0);
	std::vector<std::vector<int>> exp_mat(0, std::vector<int>(0));

	mpz_class n_sqrt(sqrt(n));

	//generate the first batch of polynomials
	int num_blocks = 0;
	mpz_vector poly_chunk = generate_poly_values(n, (num_blocks++ * block_size) - (2 * block_size), block_size);

	//get the start location of the sieve (the modular square roots) for each smooth factor
	std::vector<std::vector<int>> sieve_start = generate_solutions(n, factor_base, 1);

	int max_exp = 2;
	int_vector ext_factor_base(0);
	const size_t ext_factor_base_size = factor_base_size * (max_exp - 1);

	std::vector<std::vector<int>> ext_sieve_start(0);

	int_vector small_primes(0);
	for(size_t i = 0; i < factor_base_size / 2; i++) {
		small_primes.push_back(factor_base[i]);
	}

	for(int e = 2; e <= max_exp; e++) {
		int_vector tmp_factor_base(small_primes.size());
		for(size_t i = 0; i < small_primes.size(); i++) {
			tmp_factor_base[i] = pow(small_primes[i], e);
		}

		std::vector<std::vector<int>> tmp_solutions = generate_solutions(n, small_primes, e);

		ext_factor_base.insert(ext_factor_base.end(), tmp_factor_base.begin(), tmp_factor_base.end());
		ext_sieve_start.insert(ext_sieve_start.end(), tmp_solutions.begin(), tmp_solutions.end());
	}

	std::vector<double> factor_base_log(factor_base_size);
	for(size_t i = 0; i < factor_base_size; i++) {
		factor_base_log[i] = log(factor_base[i]);
	}

	std::vector<double> ext_factor_base_log(ext_factor_base_size);
	for(size_t i = 0; i < ext_factor_base_size; i++) {
		ext_factor_base_log[i] = log(ext_factor_base[i]);
	}

	//begin sieving until a sufficent number of bsmooth values have been found
	while(bsmooth_values.size() <= BSMOOTH_VALUE_CUTOFF) {
		std::vector<double> sieve(poly_chunk.size());

		//first sieve pass
		for(size_t i = 0; i < factor_base_size; i++) {
			for(size_t j = 0; j < sieve_start[i].size(); j++) {
				if(sieve_start[i][j] == 0) {
					continue;
				}

				size_t l;
				for(l = sieve_start[i][j]; l < poly_chunk.size(); l += factor_base[i]) {
					sieve[l] += factor_base_log[i];
				}

				sieve_start[i][j] = (l + factor_base[i]) % block_size;
			}
		}

		//second sieve pass; sieve for prime powers
		for(size_t i = 0; i < ext_sieve_start.size(); i++) {
			for(size_t j = 0; j < ext_sieve_start[i].size(); j++) {
				if(ext_sieve_start[i][j] == 0) {
					continue;
				}

				int l;
				for(l = ext_sieve_start[i][j]; l < block_size; l += ext_factor_base[i]) {
					sieve[l] += ext_factor_base_log[i];
				}

				ext_sieve_start[i][j] = (l + ext_factor_base[i]) % block_size;
			}
		}

		//filter_smooth_values(n, poly_chunk, sieve, factor_base, &bsmooth_values, &exp_mat, cutoff);
		for(int i = 0; i < block_size; i++) {
			total++;

			//int cutoff = log((2 * (i + 1)) * n_sqrt.get_ui());
			if(sieve[i] >= 0) {
				canidates++;
				int_vector exp_vec = is_bsmooth(poly_chunk[i], factor_base);
				if(exp_vec.size() != 0) {
					//std::cout << poly_chunk[i] << ", " << sieve[i] << ", " << cutoff << ", " << n_sqrt << std::endl;
					sieve_average += sieve[i];
					confirmed++;
					bsmooth_values.push_back(poly_chunk[i]);
					exp_mat.push_back(exp_vec);

					if(bsmooth_values.size() >= BSMOOTH_VALUE_CUTOFF) {
						exp_matrix = Matrix(exp_mat);
						std::cout.precision(3);
						std::cout << confirmed << " confirmed out of " << canidates << " canidates (" << ((double) confirmed / canidates) * 100 << "% passage rate), out of " << total << " total (" << ((double) canidates) * 100 / total << "% detection rate)." << std::endl;
						std::cout << "Average sieve value: " << sieve_average / confirmed << std::endl;
						return bsmooth_values;
					}
				}
			}
		}

		poly_chunk = generate_poly_values(n, (num_blocks++ * block_size) - (2 * block_size), block_size);
	}

	return bsmooth_values;
}

std::vector<std::vector<int>> Collector::generate_solutions(mpz_class n, int_vector factor_base, int power) {
	const size_t factor_base_size = factor_base.size();
	long long n_sqrt = mpz_class(sqrt(n) + 1).get_ui();
	std::vector<std::vector<int>> all_solutions(0);

	for(size_t i = 0; i < factor_base_size; i++) {
		int_vector solutions = msqrt_pow(n, factor_base[i], power);

		if(std::adjacent_find(solutions.begin(), solutions.end(), std::not_equal_to<int>()) == solutions.end()) {
			solutions[0] = modn(solutions[0] - n_sqrt, factor_base[i]);

			solutions.resize(1);
			all_solutions.push_back(solutions);
			continue;
		}

		for(size_t j = 0; j < solutions.size(); j++) {
			solutions[j] = modn(solutions[j] - n_sqrt, factor_base[i]);
		}

		all_solutions.push_back(solutions);
	}

	return all_solutions;
}

int_vector Collector::msqrt_pow(mpz_class n, int p, int e) {
	if(n % pow(p, e) == 0) {
		return {0, 0};
	}

	int_vector new_solutions(2);

	//deal with powers of 2
	if(p == 2) {
		long p_exp = pow(2, e);
		int residue = mpz_class(n % p_exp).get_ui();

		if(e == 2) {
			if(residue == 1) {
				new_solutions = {1, 3};
			}else {
				new_solutions.resize(0);
			}
		}else if(e == 3) {
			if(residue == 1) {
				new_solutions = {1, 3, 5, 7};
			}else {
				new_solutions.resize(0);
			}
		}else {
			if(n % 8 == 1) {
				new_solutions.resize(4);

				int inc = 0;
				for(int i = 1; i < p_exp; i += 2) {
					if((i * i) % p_exp == residue) {
						new_solutions[inc++] = i;
					}
				}
			}else {
				new_solutions.resize(0);
			}
		}

		return new_solutions;
	}

	int_vector orig_solutions = msqrt(n, p);
	if(e == 1) {
		return orig_solutions;
	}

	//the rest of the powers
	mpz_class p_exp(pow(p, e));
	for(int i = 0; i < 2; i++) {
		mpz_class x(orig_solutions[i]);
		mpz_class c(n);

		mpz_powm_ui(x.get_mpz_t(), x.get_mpz_t(), pow(p, e - 1), p_exp.get_mpz_t());

		int c_exp = ((pow(p, e) - (2 * pow(p, e - 1)) + 1)) / 2;
		mpz_powm_ui(c.get_mpz_t(), c.get_mpz_t(), c_exp, p_exp.get_mpz_t());

		mpz_class res((x * c) % p_exp);
		if(res.fits_sint_p()) {
			new_solutions[i] = res.get_ui();
		}else {
			throw std::overflow_error("msqrt_pow: value is too large to be converted to int.");
		}
	}

	return new_solutions;
}

//find the modular square root (the solution to r^2 === a (mod p)) with the tonelli-shanks algorithm
int_vector Collector::msqrt(mpz_class n, int p) {
	if(n % p == 0) {
		return {0, 0};
	}

	mpz_class mpz_p(p);
	mpz_class product(0);

	int s = p - 1, e = 0;
	while(s % 2 == 0) {
		s /= 2;
		e++;
	}

	mpz_class q(2);
	while(true) {
		mpz_powm_ui(product.get_mpz_t(), q.get_mpz_t(), (p - 1) / 2, mpz_p.get_mpz_t());

		if(product == p - 1) {
			break;
		}

		q++;
	}

	mpz_class x(0);
	mpz_class b(0);
	mpz_class g(0);

	mpz_powm_ui(x.get_mpz_t(), n.get_mpz_t(), (s + 1) / 2, mpz_p.get_mpz_t());
	mpz_powm_ui(b.get_mpz_t(), n.get_mpz_t(), s, mpz_p.get_mpz_t());
	mpz_powm_ui(g.get_mpz_t(), q.get_mpz_t(), s, mpz_p.get_mpz_t());

	int r = e;

	while(true) {
		int m = 0;
		while(true) {
			mpz_powm_ui(product.get_mpz_t(), b.get_mpz_t(), pow(2, m), mpz_p.get_mpz_t());
			if(product == 1) {
				break;
			}

			m++;
		}

		if(m == 0) {
			int_vector sqrts(2);
			if(x.fits_sint_p()) {
				sqrts[0] = x.get_ui();
				sqrts[1] = p - x.get_ui();
			}else {
				throw std::overflow_error("msqrt: value is too large to be converted to int.");
			}

			return sqrts;
		}

		mpz_powm_ui(product.get_mpz_t(), g.get_mpz_t(), pow(2, r - m - 1), mpz_p.get_mpz_t());
		x = (x * product) % mpz_p;

		mpz_powm_ui(g.get_mpz_t(), g.get_mpz_t(), pow(2, r - m), mpz_p.get_mpz_t());

		b = (b * g) % mpz_p;

		if(b == 1) {
			int_vector sqrts(2);
			if(x.fits_sint_p()) {
				sqrts[0] = x.get_ui();
				sqrts[1] = p - x.get_ui();
			}else {
				throw std::overflow_error("msqrt: value is too large to be converted to int.");
			}

			return sqrts;
		}

		r = m;
	}
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

	mpz_class x(poly_start);

	int i = 0;
	while(i < block_size) {
		poly_values[i] = (x * x) - n;

		if(poly_values[i] < 0) {
			poly_values[i] *= -1;
		}

		i++;
		x++;
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

	//deal with the special case of 2
	int residue = mpz_class(n % 8).get_ui();
	if(residue == 1 || residue == 7) {
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
