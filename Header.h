#pragma once
#include <complex>
#include <vector>

using u64 = unsigned long long;
using permutation = std::vector<unsigned long long>;

double LinPotential(const permutation &prm, u64 a, u64 b, unsigned n);