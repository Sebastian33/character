#include "Header.h"
#include <random>
#include <fstream>

double LinPotential(const permutation &prm, u64 a, u64 b, unsigned n)
{
	const int size(1 << n);
	const int mod(size - 1);
	const double pi = std::acos(-1);
	cmplx sum(0, 0);
	const cmplx e = std::exp(cmplx(0, 1) * pi * 2.0 / (double)size); 
	for (unsigned x = 0; x < size; x++)
	{
		sum += std::pow(e, a*x & mod) * std::conj(std::pow(e, b*prm[x] & mod));
	}
	return pow(std::abs(sum), 2);
}

permutation generatePermutation(unsigned n)
{
	const int size(1 << n);
	std::mt19937_64 engine;
	std::uniform_int_distribution<int> dstrb;
	permutation prm(size);

	for (int i = 0; i < size; i++)
	{
		prm[i] = i;
	}

	for (int i = size - 1; i > 1; i--)
	{
		std::swap(prm[i], prm[dstrb(engine) % i]);
	}
	return prm;
}

void experiment(const permutation &prm, unsigned n)
{
	const int size(1 << n);
	const int mod(size - 1);
	const double pi = std::acos(-1);
	cmplx sum(0, 0);
	const cmplx e = std::exp(cmplx(0, 1) * pi * 2.0 / (double)size);
	std::vector<std::vector<double>> data(size, std::vector<double>(size, 0));

	for (int a = 1; a < size; a++)
	{
		for (int b = 1; b < size; b++)
		{
			for (unsigned x = 0; x < size; x++)
			{
				//sum += std::pow(e, a*x & mod) * std::conj(std::pow(e, b*prm[x] & mod));
				sum += std::pow(e, (a*x - b * prm[x]) & mod); //better, faster, stronger ...
			}
			data[a - 1][b - 1] = std::norm(sum);
			sum = 0;
		}
	}

	std::ofstream out("data.txt");
	for (int i = 1; i < size; i++)
	{
		out << ' ' << i;
	}
	out << std::endl;

	for (int a = 1; a < size; a++)
	{
		out << a << ' ';
		for (int b = 1; b < size; b++)
		{
			out << data[a - 1][b - 1] << ' ';
		}
		out << std::endl;
	}
	out.close();
}