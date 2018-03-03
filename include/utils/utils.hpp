/*A collection of useful functions and other random things*/

#pragma once

//std includes
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

//gmp includes
#include <gmp.h>
#include <gmpxx.h>

//typedefs
typedef std::vector<mpz_class> mpz_vector;
typedef std::vector<int> int_vector;

/*Print Functions*/

//overload the << operator to print a vector
template<typename T> std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    out << "[";

    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];

		if (i != last) {
            out << ", ";
		}
    }

    out << "]";

	fflush(stdout);

    return out;
}

bool check(mpz_class n, mpz_vector factors);

/*Math Functions*/

//get all primes less than n with the Sieve of Eratosthenes
int_vector get_primes(int n);

//calculate the natural logarithim (ln) of a number
double ln(mpz_class n);
