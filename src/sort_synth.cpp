#include <iostream>
#include <chrono>
#include <thread>

#include <cstring>
#include <unistd.h>

#include "utils.h"
#include "distr.h"

//#define SORT_CHECKS

#if defined ULS600 || defined ULS1000
#include "bls.h"
#include "indexes/linear.h"

#elif defined BLS
#include "bls.h"
#include "indexes/linear.h"
#include "indexes/linear_balanced.h"

#elif defined LS21
#include "bls.h"
#include "indexes/linear.h"
#include "indexes/my_rmi.h"

#elif defined IPS
#include "ips4o.hpp"

#elif defined LS20
#include "learned_sort.h"

#elif defined PDQ
#include "pdqsort.h"

#elif defined RS
#include "radix_sort.h"

#elif defined SKA
#include "ska_sort.hpp"

#elif defined STDSORT
#include <algorithm>

#else
#error Unspecified Algorithm
#endif

typedef double data_t;
using vec_iter = std::vector<data_t>::iterator;

template <typename F>
void benchmark(std::string benchmark_name, const size_t size, const distr_t distribution, F && fun) {  
  constexpr long runs = 25;
  long times[runs];

  for (size_t i = 0; i < runs; ++i)
  {
    const uint seed = i;
    std::vector<data_t> v(size);

    // Cache the zipf distribution (slow generation)
    if (distribution == distr_t::ZIPF) {
      std::ifstream ifp("zipf_" + std::to_string(size) + ".cache", std::ios::in | std::ios::binary);
      if (ifp.good()) {
        ifp.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(v[0]));
        ifp.close();
      } else {
        generate_data(distribution, v.begin(), v.end(), seed);
        std::ofstream ofp("zipf_" + std::to_string(size) + ".cache", std::ios::out | std::ios::binary);
        ofp.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(v[0]));
        ofp.close();
      }
    }
    else {
      generate_data(distribution, v.begin(), v.end(), seed);
    }

#ifdef SORT_CHECKS
    auto checksum = get_checksum(v);
#endif
    
    // Start the timer and execute the function
    const auto start = std::chrono::steady_clock::now();
    fun(v.begin(), v.end());
    const auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    times[i] = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    
#if SORT_CHECKS
    if (get_checksum(v) != checksum) {
      std::cout << "Checksum ERROR" << std::endl;
      std::cerr << "Checksum ERROR" << std::endl;
      return;
    }

    if (!std::is_sorted(v.begin(), v.end())) {
      std::cout << "Not sorted ERROR" << std::endl;
      std::cerr << "Not sorted ERROR" << std::endl;
      return;
    }
#endif
  }

  // Compute avg, median and stddev
  double mean = std::accumulate(times, times + runs, 0.0) / runs;

  auto cond_round = [mean] (double val) -> double {return mean > 100'000 ? round(val / 1000.) : val / 1000.; };
  benchmark_name += "," + std::to_string(distribution);
  for (size_t i = 0; i < runs; ++i)
    std::cout << benchmark_name << "," << (i+1) << "," << cond_round(times[i]) << ",ms" << std::endl;

  double std_dev = 0;
  for(size_t i = 0; i < runs; ++i)
    std_dev += std::pow(times[i] - mean, 2);
  std_dev = std::sqrt(std_dev / runs);
  std::sort(times, times + runs);
  const double median = runs % 2 != 0 ? times[runs/2] : (times[runs/2] + times[runs/2 - 1]) / 2;

  std::cout << benchmark_name << ",mean," << cond_round(mean) << ",ms" << std::endl;
  std::cout << benchmark_name << ",median," << cond_round(median) << ",ms" << std::endl;
  std::cout << benchmark_name << ",stddev," << cond_round(std_dev) << ",ms" << std::endl;
}

int stick_this_thread_to_core(int core_id) {
   int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
   if (core_id < 0 || core_id >= num_cores)
      return EINVAL;

   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(core_id, &cpuset);

   pthread_t current_thread = pthread_self();    
   return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

int main(int argc, char ** argv) {
  
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " input_size_mln dup/nodup" << std::endl;
    return -1;
  }

  const size_t size = std::atol(argv[1]) * 1'000'000;
  const bool dup = (std::strcmp(argv[2], "dup") == 0);

  if (size > 2000 * 1'000'000) {
    std::cerr << "WARNING: Error the input is too big: " << size << std::endl;
    return -1;
  }

  if (dup == false && std::strcmp(argv[2], "nodup") != 0) {
    std::cerr << "ERROR: dup/nodup unspecified" << std::endl;
    return -1;
  }

  // Before starting the benchmark pin this thread to a core
  bool pinned = true;
  if (stick_this_thread_to_core(34) != 0) {
    std::cerr << "WARNING: cannot pin this thread to a core" << std::endl;
    pinned = false;
  }

  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  std::cout << "Synth benchmark " << (dup ? "dup" : "nodup") << ", Size: " << size << ", Pinning: " << (pinned ? "yes" : "no") << ", host: " << hostname << std::endl;

  std::vector<distr_t> distributions;
  if (dup) distributions = {EIGHT_DUPS, MODULO, ROOT_DUPS, TWO_DUPS};
  else distributions = {CHI_SQUARED, EXPONENTIAL, LOGNORMAL, MIX_GAUSS, NORMAL, UNIFORM, ZIPF};
  
  std::cout << "round,algo,distribution,run,time,unit\n";

   const int round = 0;
   {
    std::string str_r = std::to_string(round) + ",";
    for (int d = 0; d < distributions.size(); ++d) {
#if defined ULS1000
      benchmark(str_r + "ULS1000", size, distributions[d], [](auto first, auto last) { learned_sort_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 1000>(first, last); });
#elif defined ULS600
      benchmark(str_r + "ULS600", size, distributions[d], [](auto first, auto last) { learned_sort_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 600>(first, last); });
#elif defined BLS
      benchmark(str_r + "BLS", size, distributions[d], [](auto first, auto last) { learned_sort_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 600>(first, last); });
#elif defined IPS
      benchmark(str_r + "IPS4O", size, distributions[d], [](auto first, auto last) { ips4o::sort(first, last); });
#elif defined STDSORT
      benchmark(str_r + "std::sort", size, distributions[d], [](auto first, auto last) { std::sort(first, last); });
#elif defined SKA
      benchmark(str_r + "SkaSort", size, distributions[d], [](auto first, auto last) { ska_sort(first, last); });
#elif defined LS20
      benchmark(str_r + "LearnedSort_2_0", size, distributions[d], [](auto first, auto last) { learned_sort::sort(first, last); });
#elif defined LS21
      benchmark(str_r + "LearnedSort_2_1", size, distributions[d], [](auto first, auto last) { learned_sort_framework::sort<indexes::MyRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 500>(first, last); });
#elif defined RS
      benchmark(str_r + "RadixSort", size, distributions[d], [](auto first, auto last) { radix_sort(first, last); });
#elif defined PDQ
      benchmark(str_r + "PDQS", size, distributions[d], [](auto first, auto last) { pdqsort(first, last); });
#endif

      std::this_thread::sleep_for(std::chrono::seconds(10));
    }
  }

  return 0;
}
