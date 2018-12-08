#include "Header.h"
#include <random>
#include <fstream>

const std::vector<int> ord{ 
	0, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256,
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256,
	8, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	4, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	8, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	2, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	8, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	4, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	8, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256, 
	16, 256, 128, 256, 64, 256, 128, 256, 32, 256, 128, 256, 64, 256, 128, 256 };

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

table experiment(const permutation &prm, unsigned n)
{
	const int size(1 << n);
	const int mod(size - 1);
	const double pi = std::acos(-1);
	cmplx sum(0, 0);
	const cmplx e = std::exp(cmplx(0, 1) * pi * 2.0 / (double)size);
	table data(size, std::vector<double>(size, 0));

	for (unsigned a = 1; a < size; a++)
	{
		for (unsigned b = 1; b < size; b++)
		{
			for (unsigned x = 0; x < size; x++)
			{
				//sum += std::pow(e, a*x & mod) * std::conj(std::pow(e, b*prm[x] & mod));
				sum += std::pow(e, (a*x - b * prm[x]) & mod); //better, faster, stronger ...
			}
			data[a][b] = std::norm(sum);
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
			out << data[a][b] << ' ';
		}
		out << std::endl;
	}
	out.close();
	return data;
}

void distribution(const table &data, unsigned n)
{
	const int absval(1 << 2*n);
	const int size(1 << n);
	int precision(20);
	std::vector<int> dstrb(absval / precision, 0);

	int j;
	int end(3400);
	dstrb[0] = 0;
	for (int i = precision; i < end; i+=precision)
	{
		j = i / precision;
		for (int r = 1; r < size; r++)
		{
			for (int c = 1; c < size; c++)
			{
				dstrb[j] += static_cast<int>(data[r][c] < i);
			}
		}
	}

	std::ofstream out("distribution.txt");
	for (auto d : dstrb)
		out << d << ' ';
	out.close();
}

double parkTheorem(const table &potentials, unsigned n, unsigned brnchIndx)
{
	const int size(1 << n);
	const int denum(size*size);
	double sum1(0), tmp;
	for (int a = 1; a < size; a++)
	{
		tmp = 0;
		for (int b = 1; b < size; b++)
		{
			tmp += (ord[b] - 1)*pow(potentials[a][b], brnchIndx);
		}
		if (tmp > sum1)
			sum1 = tmp;
	}

	double sum2(0);
	for (int b = 1; b < size; b++)
	{
		tmp = 0;
		for (int a = 1; a < size; a++)
		{
			tmp += (ord[a] - 1)*pow(potentials[a][b], brnchIndx);
		}
		if (tmp > sum2)
			sum2 = tmp;
	}
	return ((sum1 > sum2) ? sum1 : sum2) / pow(denum, brnchIndx);
}

table fromFile(unsigned n)
{
	int size(1 << n);
	std::ifstream in("data.txt");
	int tmp;

	for (int i = 1; i < size; i++)
		in >> tmp;

	table data(size, std::vector<double>(size, 0));
	for (int a = 1; a < size; a++)
	{
		in >> tmp;
		for (int b = 1; b < size; b++)
		{
			in >> data[a][b];
		}
	}
	in.close();
	return data;
}