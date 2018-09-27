#include "Header.h"


double LinPotential(const permutation &prm, u64 a, u64 b, unsigned n)
{
	const int size(1 << n);
	const int mod(size - 1);
	const double pi = std::acos(-1);
	std::complex<double> sum(0, 0);
	const std::complex<double> i(0, 1);
	const std::complex<double> e = std::exp(i*pi * 2.0 / (double)size); 
	for (unsigned x = 1; x < size; x++)
	{
		sum += std::pow(e, a*x & mod)*std::conj(std::pow(e, b*prm[x] & mod));
	}
	return std::abs(sum);
}