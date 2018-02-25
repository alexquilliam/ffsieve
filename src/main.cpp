//local includes
#include <utils.hpp>
#include <decomposer.hpp>

int main(int argc, char **argv) {
	Decomposer decomposer;
	int_vector prime_factors = decomposer.decompose(100);
	std::cout << prime_factors << std::endl;

	return 0;
}
