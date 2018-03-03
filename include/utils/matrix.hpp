/*A matrix class with relevant utility functions*/

#pragma once

//std includes
#include <bitset>

//local includes
#include <utils.hpp>

//overload the << operator to print a matrix
template<typename T> std::ostream& operator<< (std::ostream& out, const std::vector<std::vector<T>>& mat) {
	for(size_t i = 0; i < mat.size(); i++) {
		std::cout << mat[i] << std::endl;
	}

	fflush(stdout);

    return out;
}

class Matrix : public std::vector<std::vector<int>> {
	private:
		std::vector<int> null_space;
	public:
		Matrix(std::vector<std::vector<int>> vec);
		Matrix();
		void mod();
		void transpose();
		void to_rref();
		std::vector<int> get_null_space();
};
