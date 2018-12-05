#pragma once
#include <complex>
#include <vector>

using u64 = unsigned long long;
using permutation = std::vector<unsigned long long>;
using cmplx = std::complex<double>;
using table = std::vector<std::vector<double>>;

double LinPotential(const permutation &prm, u64 a, u64 b, unsigned n);
permutation generatePermutation(unsigned n);
table experiment(const permutation &prm, unsigned n);
void distribution(const table &data, unsigned n);