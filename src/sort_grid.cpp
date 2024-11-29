
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <vector>
#include <chrono>

#include "utils.h"
#include "distr.h"

#include "bls/bls.h"

#if defined ULS
#include "bls/indexes/linear.h"

#elif defined BLS
#include "bls/indexes/linear.h"
#include "bls/indexes/linear_balanced.h"

#elif defined LS21
#include "bls/indexes/linear.h"
#include "bls/indexes/new_rmi.h"

#else
#error Unspecified learned algo

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
    else
      generate_data(distribution, v.begin(), v.end(), seed);
    
    // Start the timer and execute the function
    auto start = std::chrono::steady_clock::now();
    fun(std::move(v));
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    times[i] = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  }

  // Compute avg, median and stddev
  double mean = std::accumulate(times, times + runs, 0.0) / runs;
  benchmark_name += "," + std::to_string(distribution);

  auto cond_round = [mean] (double val) -> double {return mean > 100'000 ? round(val / 1000.) : val / 1000.; };
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
  
  std::vector<distr_t> distributions = {CHI_SQUARED, EXPONENTIAL, LOGNORMAL, MIX_GAUSS, NORMAL, UNIFORM, ZIPF};
  
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " input_size_mln" << std::endl;
    return -1;
  }

  const size_t size = std::atol(argv[1]) * 1'000'000;
  
  if (size > 2000 * 1'000'000) {
    std::cerr << "WARNING: Error the input is too big: " << size << std::endl;
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
  std::cout << "Grid benchmark -> Size: " << size << ", Pinning: " << (pinned ? "yes" : "no") << ", host: " << hostname << std::endl;
  
  std::cout << "round,algo,distribution,run,time,unit\n";

  const int round = 0;
  {
    std::string str_r = std::to_string(round) + ",";
    for (int d = 0; d < distributions.size(); ++d) {

#if defined LS21
      benchmark(str_r + "rmi_" + std::to_string(300), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,   300>(v.begin(), v.end());});
      benchmark(str_r + "rmi_" + std::to_string(400), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,   400>(v.begin(), v.end());});
      benchmark(str_r + "rmi_" + std::to_string(500), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,   500>(v.begin(), v.end());});
      benchmark(str_r + "rmi_" + std::to_string(600), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,   600>(v.begin(), v.end());});
      benchmark(str_r + "rmi_" + std::to_string(700), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,   700>(v.begin(), v.end());});
      benchmark(str_r + "rmi_" + std::to_string(800), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,   800>(v.begin(), v.end());});
      benchmark(str_r + "rmi_" + std::to_string(900), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,   900>(v.begin(), v.end());});
      benchmark(str_r + "rmi_" + std::to_string(1000), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 1000>(v.begin(), v.end());});
      
#elif defined BLS
      benchmark(str_r + "lbm_" + std::to_string( 300), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  300>(v.begin(), v.end());});
      benchmark(str_r + "lbm_" + std::to_string( 400), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  400>(v.begin(), v.end());});
      benchmark(str_r + "lbm_" + std::to_string( 500), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  500>(v.begin(), v.end());});
      benchmark(str_r + "lbm_" + std::to_string( 600), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  600>(v.begin(), v.end());});
      benchmark(str_r + "lbm_" + std::to_string( 700), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  700>(v.begin(), v.end());});
      benchmark(str_r + "lbm_" + std::to_string( 800), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  800>(v.begin(), v.end());});
      benchmark(str_r + "lbm_" + std::to_string( 900), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  900>(v.begin(), v.end());});
      benchmark(str_r + "lbm_" + std::to_string(1000), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 1000>(v.begin(), v.end());});
      
#elif defined ULS
      benchmark(str_r + "linear_" + std::to_string( 300), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  300>(v.begin(), v.end());});
      benchmark(str_r + "linear_" + std::to_string( 400), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  400>(v.begin(), v.end());});
      benchmark(str_r + "linear_" + std::to_string( 500), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  500>(v.begin(), v.end());});
      benchmark(str_r + "linear_" + std::to_string( 600), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  600>(v.begin(), v.end());});
      benchmark(str_r + "linear_" + std::to_string( 700), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  700>(v.begin(), v.end());});
      benchmark(str_r + "linear_" + std::to_string( 800), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  800>(v.begin(), v.end());});
      benchmark(str_r + "linear_" + std::to_string( 900), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>,  900>(v.begin(), v.end());});
      benchmark(str_r + "linear_" + std::to_string(1000), size, distributions[d], [](std::vector<data_t> && v) {ls_framework::sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 1000>(v.begin(), v.end());});

#endif

    }
  }

  return 0;
}
