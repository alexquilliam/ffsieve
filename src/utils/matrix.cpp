/*A matrix class with relevant utility functions*/

#include <matrix.hpp>

/*Constructors*/

Matrix::Matrix() {};

Matrix::Matrix(std::vector<std::vector<int>> vec) {
	for(size_t i = 0; i < vec.size(); i++) {
		this -> push_back(vec[i]);
	}
}

/*Public methods*/

void Matrix::mod() {
	for(size_t i = 0; i < this -> size(); i++) {
		for(size_t j = 0; j < this -> at(0).size(); j++) {
			this -> at(i).at(j) %= 2;
		}
	}
}

void Matrix::transpose() {
	std::vector<std::vector<int>> vec(this -> at(0).size(), std::vector<int>(this -> size()));

	for(size_t i = 0; i < this -> size(); i++) {
		for(size_t j = 0; j < this -> at(0).size(); j++) {
			vec[j][i] = this -> at(i).at(j);
		}
	}

	Matrix mat(vec);

	this -> clear();
	Matrix().swap(*this);
	this -> swap(mat);
}

void Matrix::to_rref() {
	std::vector<std::vector<int>> vec = *this;
	int m = vec.size(); //org: m=vals, trans: m=exps
	int n = vec.at(0).size(); //org: n=exps, trans: n=vals

	if(m < n) {
		std::cerr << "\nWARNING: There are more columns than rows in the input matrix.\n" << std::endl;
	}

	//reduction
	for(int pivot = 0; pivot < n; pivot++) {
		if(vec.at(pivot).at(pivot) == 0) {
			for(int i = pivot; i < n; i++) {
				if(vec.at(i).at(pivot) == 1) {
					vec.at(i).swap(vec.at(pivot));
				}
			}
		}

		for(int i = pivot + 1; i < n; i++) {
			if(vec.at(i).at(pivot) == 1) {
				for(int j = 0; j < m; j++) {
					vec.at(i).at(j) ^= vec.at(pivot).at(j);
				}
			}
		}
	}

	//back substitution
	int line;
	for(line = 0; line < m; line++) {
		if(vec.at(line).at(line) != 1) {
			break;
		}
	}

	for(int i = line - 1; i >= 0; i--) {
		for(int j = i - 1; j >= 0; j--) {
			if(vec.at(j).at(i)) {
				for(int k = 0; k < m; k++) {
					vec.at(j).at(k) ^= vec.at(i).at(k);
				}
			}
		}
	}

	std::vector<int> null_space(m);
	null_space.at(line) = 1;

	for(int i = 0; i < line; i++) {
		null_space.at(i) = vec.at(i).at(line);
	}

	this -> null_space = null_space;

	//put results back into the matrix
	for(int i = 0; i < m; i++) {
		this -> at(i).swap(vec.at(i));
	}
}

std::vector<int> Matrix::get_null_space() {
	return this -> null_space;
}
