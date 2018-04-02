/*A matrix class with relevant utility functions*/

#include <matrix.hpp>

/*Constructors*/

#define vthis (*this)

Matrix::Matrix() {};

Matrix::Matrix(std::vector<std::vector<int>> vec) {
	vthis.resize(vec.size());
	for(size_t i = 0, vec_size = vec.size(); i < vec_size; i++) {
		vthis[i] = vec[i];
	}
}

/*Public methods*/

void Matrix::mod() {
	for(size_t i = 0, n = vthis[0].size(); i < n; i++) {
		for(size_t j = 0, m = vthis.size(); j < m; j++) {
			vthis[j][i] %= 2;
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
	size_t m = vthis.size();
	size_t n = vthis[0].size();

	std::vector<std::vector<int>> mat = vthis;

	int current_row = 0;
	int pivot = 0, next_pivot = 0;
	for(size_t c = 0; c < n; c++) {
		next_pivot = -1;

		for(size_t i = current_row; i < m; i++) {
			if(mat[i][c]) {
				next_pivot = i;
				break;
			}
		}

		if(next_pivot == -1) {
			continue;
		}

		if(next_pivot != current_row) {
			mat[next_pivot].swap(mat[current_row]);
		}

		pivot = current_row;
		current_row++;

		for(size_t i = current_row; i < m; i++) {
			if(mat[i][c]) {
				for(size_t j = 0; j < n; j++) {
					mat[i][j] ^= mat[pivot][j];
				}
			}
		}
	}

	size_t line;
	for(line = 0; line < m; line++) {
		if(mat[line][line] != 1) {
			break;
		}
	}

	for(int i = line - 1; i >= 0; i--) {
		for(int j = i - 1; j >= 0; j--) {
			if(mat[j][i]) {
				for(size_t k = 0; k < m; k++) {
					mat[j][k] ^= mat[i][k];
				}
			}
		}
	}

	std::vector<int> null_space(m);
	null_space[line] = 1;

	for(size_t i = 0; i < line; i++) {
		null_space[i] = mat[i][line];
	}

	for(size_t i = 0; i < null_space.size(); i++) {
		if(null_space[i] != 0) {
			null_space[i] = 1;
		}
	}

	vthis.null_space = null_space;

	for(size_t i = 0; i < m; i++) {
		vthis[i].swap(mat[i]);
	}
}

std::vector<int> Matrix::get_null_space() {
	return vthis.null_space;
}
