#pragma once
#include <complex>
#include <vector>
//TODO: make as a class
using u64 = unsigned long long;
using permutation = std::vector<unsigned long long>;
using cmplx = std::complex<double>;

double LinPotential(const permutation &prm, u64 a, u64 b, unsigned n);
permutation generatePermutation(unsigned n);
void experiment(const permutation &prm, unsigned n);