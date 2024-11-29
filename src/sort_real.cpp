
#include <iostream>
#include <cstring>
#include <cmath>
#include <chrono>
#include <thread>
#include <numeric>
#include <random>
#include <unistd.h>

#include "utils.h"

#if defined ULS600 || defined ULS1000
#include "bls/bls.h"
#include "bls/indexes/linear.h"

#elif defined BLS
#include "bls/bls.h"
#include "bls/indexes/linear.h"
#include "bls/indexes/linear_balanced.h"

#elif defined LS21
#include "bls/bls.h"
#include "bls/indexes/linear.h"
#include "bls/indexes/new_rmi.h"

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

std::vector<std::string> datasets_unsigned = {
  "../data/Chic/Start.txt",
  "../data/NYC/Pickup.txt",
  "../data/Stks/Date.txt",
  "../data/books_200M_uint64",
  "../data/books_400M_uint64",
  "../data/books_600M_uint64",
  "../data/books_800M_uint64",
  "../data/fb_200M_uint64_cleaned",
  "../data/osm_cellids_200M_uint64",
  "../data/osm_cellids_400M_uint64",
  "../data/osm_cellids_600M_uint64",
  "../data/osm_cellids_800M_uint64",
  "../data/wiki_ts_200M_uint64"
};
const size_t start_unsigned = 0;

std::vector<std::string> datasets_double = {
  "../data/Chic/Tot.txt",
  "../data/NYC/Dist.txt",
  "../data/NYC/Tot.txt",
  "../data/Sof/Hum.txt",
  "../data/Sof/Press.txt",
  "../data/Sof/Temp.txt",
  "../data/Stks/Low.txt",
  "../data/Stks/Open.txt",
  "../data/Stks/Vol.txt"
};
const size_t start_double = start_unsigned + datasets_unsigned.size();

template <typename data_t, typename F>
void benchmark(std::string benchmark_name, const size_t dataset_index, F && fun) {  
  constexpr long runs = 25;
  long times[runs];

  std::vector<data_t> input;

  // Read from file
  std::vector<std::string>& dataset = typeid(data_t) == typeid(double) ? datasets_double : datasets_unsigned;
  const size_t start_point = typeid(data_t) == typeid(double) ? start_double : start_unsigned;
  std::string data_path = dataset[dataset_index - start_point];
  std::ios_base::openmode mode = data_path.ends_with(".txt") ? std::ios::in : std::ios::in | std::ios::binary;
  std::ifstream is(data_path, mode);

  if (!is.is_open()) {
    std::cerr << "File error " << dataset_index << std::endl; 
    return;
  }

  if (!data_path.ends_with(".txt")) {
    data_t size; 
    is.read((char *) &size, sizeof(data_t));
    input.resize(size);
    is.read((char *) input.data(), size * sizeof(data_t));
  } else {
    std::istream_iterator<data_t> start(is), end;
    input.clear();
    input.insert(input.begin(), start, end);
  }

  if (std::is_sorted(input.begin(), input.end())) {
    constexpr const uint seed = 42;
    std::shuffle(input.begin(), input.end(), std::mt19937(seed));
  }

  for (size_t i = 0; i < runs; ++i)
  {
    // Copy the input vector in v
    std::vector<data_t> v = input;
    
    // Start the timer and execute the function
    auto start = std::chrono::steady_clock::now();
    fun(std::move(v));
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    times[i] = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  }

  // Compute avg, median and stddev
  double mean = std::accumulate(times, times + runs, 0.0) / runs;
  benchmark_name += "," + std::to_string(dataset_index);

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
  
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << std::endl;
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
  std::cout << "Real benchmark -> Pinning: " << (pinned ? "yes" : "no") << ", host: " << hostname << std::endl;
  
  using vec_iter_u = std::vector<uint64_t>::iterator;
  using vec_iter_d = std::vector<double>::iterator;

  std::cout << "round,algo,distribution,run,time,unit\n";

  const int round = 0;
  {
    std::string str_r = std::to_string(round) + ",";
    for (int i = start_unsigned; i < start_unsigned + datasets_unsigned.size(); ++i) {
#if defined IPS
      benchmark<uint64_t>(str_r + "IPS4O", i, [](std::vector<uint64_t> && v) { ips4o::sort(v.begin(), v.end()); });
#elif defined SKA
      benchmark<uint64_t>(str_r + "SkaSort", i, [](std::vector<uint64_t> && v) { ska_sort(v.begin(), v.end()); });
#elif defined LS20 
      benchmark<uint64_t>(str_r + "LearnedSort_2_0", i, [](std::vector<uint64_t> && v) { learned_sort::sort(v.begin(), v.end()); });
#elif defined ULS1000
      benchmark<uint64_t>(str_r + "ULS1000", i, [](std::vector<uint64_t> && v) { ls_framework::sort<indexes::MinMaxIndex<vec_iter_u>, indexes::MinMaxIndex<vec_iter_u>, 1000>(v.begin(), v.end()); });
#elif defined ULS600
      benchmark<uint64_t>(str_r + "ULS600", i, [](std::vector<uint64_t> && v) { ls_framework::sort<indexes::MinMaxIndex<vec_iter_u>, indexes::MinMaxIndex<vec_iter_u>, 600>(v.begin(), v.end()); });
#elif defined BLS
      benchmark<uint64_t>(str_r + "BLS", i, [](std::vector<uint64_t> && v) { ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter_u>, indexes::MinMaxIndex<vec_iter_u>, 600>(v.begin(), v.end()); });
#elif defined LS21
      benchmark<uint64_t>(str_r + "LearnedSort_2_1", i, [](std::vector<uint64_t> && v) { ls_framework::sort<indexes::NewRmiIndex<vec_iter_u>, indexes::MinMaxIndex<vec_iter_u>, 600>(v.begin(), v.end()); });
#elif defined RS
      benchmark<uint64_t>(str_r + "RadixSort", i, [](std::vector<uint64_t> && v) { radix_sort(v.begin(), v.end()); });
#elif defined PDQ
      benchmark<uint64_t>(str_r + "PDQS", i, [](std::vector<uint64_t> && v) { pdqsort(v.begin(), v.end()); });
#elif defined STDSORT
      benchmark<uint64_t>(str_r + "std::sort", i, [](std::vector<uint64_t> && v) { std::sort(v.begin(), v.end()); });
#endif

      std::this_thread::sleep_for(std::chrono::seconds(10));
    }

    for (int i = start_double; i < start_double + datasets_double.size(); ++i) {
#if defined IPS
      benchmark<double>(str_r + "IPS4O", i, [](std::vector<double> && v) { ips4o::sort(v.begin(), v.end()); });
#elif defined SKA
      benchmark<double>(str_r + "SkaSort", i, [](std::vector<double> && v) { ska_sort(v.begin(), v.end()); });
#elif defined LS20
      benchmark<double>(str_r + "LearnedSort_2_0", i, [](std::vector<double> && v) { learned_sort::sort(v.begin(), v.end()); });
#elif defined ULS1000
      benchmark<double>(str_r + "ULS1000", i, [](std::vector<double> && v) { ls_framework::sort<indexes::MinMaxIndex<vec_iter_d>, indexes::MinMaxIndex<vec_iter_d>, 1000>(v.begin(), v.end()); });
#elif defined ULS600
      benchmark<double>(str_r + "ULS600", i, [](std::vector<double> && v) { ls_framework::sort<indexes::MinMaxIndex<vec_iter_d>, indexes::MinMaxIndex<vec_iter_d>, 600>(v.begin(), v.end()); });
#elif defined BLS
      benchmark<double>(str_r + "BLS", i, [](std::vector<double> && v) { ls_framework::sort<indexes::BalancedMinMaxIndex<vec_iter_d>, indexes::MinMaxIndex<vec_iter_d>, 600>(v.begin(), v.end()); });
#elif defined LS21
      benchmark<double>(str_r + "LearnedSort_2_1", i, [](std::vector<double> && v) { ls_framework::sort<indexes::NewRmiIndex<vec_iter_d>, indexes::MinMaxIndex<vec_iter_d>, 600>(v.begin(), v.end()); });
#elif defined RS
      benchmark<double>(str_r + "RadixSort", i, [](std::vector<double> && v) { radix_sort(v.begin(), v.end()); });
#elif defined PDQ
      benchmark<double>(str_r + "PDQS", i, [](std::vector<double> && v) { pdqsort(v.begin(), v.end()); });
#elif defined STDSORT
      benchmark<double>(str_r + "std::sort", i, [](std::vector<double> && v) { std::sort(v.begin(), v.end()); });
#endif

      std::this_thread::sleep_for(std::chrono::seconds(10));
    }
  }

  return 0;
}
