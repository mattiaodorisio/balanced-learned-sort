
#include <iostream>
#include <vector>
#include <random>

#include "utils.h"
#include "bls.h"

int main(int argc, char ** argv) {
  
  constexpr size_t sz = 1'000'000;
  constexpr uint_fast32_t seed = 42;

  std::vector<double> v(sz);
  
  std::mt19937 rnd(seed);
  std::exponential_distribution<> dist(2.);
  
  // Populate the input
  for (size_t i = 0; i < v.size(); ++i) {
    v[i] = dist(rnd);
  }

  ls_framework::bls<std::vector<double>::iterator>(v.begin(), v.end());

  std::cout << (std::is_sorted(v.begin(), v.end()) ? "Sorted" : "Not Sorted") << std::endl;

  return 0;
}
