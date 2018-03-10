/*A matrix class with relevant utility functions*/

#include <matrix.hpp>

/*Constructors*/

#define vthis (*this)

Matrix::Matrix() {};

Matrix::Matrix(std::vector<std::vector<int>> vec) {
	for(size_t i = 0, vec_size = vec.size(); i < vec_size; i++) {
		vthis.push_back(vec[i]);
	}
}

/*Public methods*/

void Matrix::mod() {
	for(size_t i = 0, n = vthis[0].size(); i < n; i++) {
		for(size_t j = 0, m = vthis.size(); j < m; j++) {
			vthis[i][j] %= 2;
		}
	}
}

void Matrix::transpose() {
	std::vector<std::vector<int>> vec(vthis[0].size(), std::vector<int>(vthis.size()));

	for(size_t i = 0, m = vthis.size(); i < m; i++) {
		for(size_t j = 0, n = vthis[0].size(); j < n; j++) {
			vec[j][i] = vthis[i][j];
		}
	}

	Matrix(vec).swap(vthis);
}

void Matrix::to_rref() {
	std::vector<std::vector<int>> vec = vthis;
	int m = vec.size();
	int n = vec[0].size();

	//reduction
	for(int pivot = 0; pivot < n; pivot++) {
		if(vec[pivot][pivot] == 0) {
			for(int i = pivot; i < n; i++) {
				if(vec[i][pivot] == 1) {
					vec[i].swap(vec[pivot]);
				}
			}
		}

		for(int i = pivot + 1; i < n; i++) {
			if(vec[i][pivot] == 1) {
				for(int j = 0; j < m; j++) {
					vec[i][j] ^= vec[pivot][j];
				}
			}
		}
	}

	//back substitution
	int line;
	for(line = 0; line < m; line++) {
		if(vec[line][line] != 1) {
			break;
		}
	}

	for(int i = line - 1; i >= 0; i--) {
		for(int j = i - 1; j >= 0; j--) {
			if(vec[j][i]) {
				for(int k = 0; k < m; k++) {
					vec[j][k] ^= vec[i][k];
				}
			}
		}
	}

	std::vector<int> null_space(m);
	null_space[line] = 1;

	for(int i = 0; i < line; i++) {
		null_space[i] = vec[i][line];
	}

	vthis.null_space = null_space;

	//put results back into the matrix
	for(int i = 0; i < m; i++) {
		vthis[i].swap(vec[i]);
	}
}

std::vector<int> Matrix::get_null_space() {
	return vthis.null_space;
}
