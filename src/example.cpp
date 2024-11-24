
#include <iostream>
#include <cstring>
#include <unistd.h>

#include "utils.h"

#include "learned_sort.h"
#include "bls.h"
#include "ips4o.hpp"

#include "indexes/linear.h"
#include "indexes/linear_balanced.h"
#include "indexes/my_rmi.h"

#include "pdqsort.h"
#include "radix_sort.h"
#include "ska_sort.hpp"

template <typename data_t>
struct learned_cfg {
  static constexpr const size_t kMaxBuckets = 1000;
  using iterator = std::vector<data_t>::iterator;
  using bucket_type = uint;
  using less = std::less<>;
  static constexpr const int logBuckets(const std::ptrdiff_t n) {
    return 10;
  }
};

template <typename data_t, size_t m>
using lbm_t = indexes::BalancedMinMaxIndex<learned_cfg<data_t>, m>;

std::vector<std::string> datasets_unsigned = {
  "../../data/books_200M_uint32"
};
const size_t start_unsigned = 0;

template <typename data_t, typename F>
void benchmark(std::string benchmark_name, const size_t dataset_index, F && fun) {
  constexpr const uint seed = 42;
  
  constexpr long runs = 5;
  long times[runs];

  for (size_t i = 0; i < runs; ++i)
  {
    std::vector<data_t> v;

    // Read from file
    std::vector<std::string>& dataset = datasets_unsigned;
    const size_t start_point = start_unsigned;
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
    	v.resize(size);
   	  is.read((char *) v.data(), size * sizeof(data_t));
    } else {
    	std::istream_iterator<data_t> start(is), end;
    	v.clear();
    	v.insert(v.begin(), start, end);
    }

    if (std::is_sorted(v.begin(), v.end()))
      std::shuffle(v.begin(), v.end(), std::mt19937(seed));
    
    // Start the timer and execute the function
    const auto start = std::chrono::steady_clock::now();
    fun(std::move(v));
    const auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    times[i] = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  }

  // Compute avg, median and stddev
  double mean = std::accumulate(times, times + runs, 0.0) / runs;
  double std_dev = 0;
  for(size_t i = 0; i < runs; ++i)
    std_dev += std::pow(times[i] - mean, 2);
  std_dev = std::sqrt(std_dev / runs);
  std::nth_element(times, times + runs/2, times + runs);
  double median = times[runs/2];
  if (runs % 2 == 0) {
    std::nth_element(times, times + runs/2+1, times + runs);
    median = (median + times[runs/2+1]) / 2;
  }

  benchmark_name += "_" + std::to_string(dataset_index);

  auto cond_round = [&] (double val) -> double {return mean > 100'000 ? round(val / 1000.) : val / 1000.; };
  for (size_t i = 0; i < runs; ++i)
    std::cout << benchmark_name << "," << (i+1) << "," << cond_round(times[i]) << ",ms" << std::endl;

  std::cout << benchmark_name << "_mean," << cond_round(mean) << ",ms" << std::endl;
  std::cout << benchmark_name << "_median," << cond_round(median) << ",ms" << std::endl;
  std::cout << benchmark_name << "_stddev," << cond_round(std_dev) << ",ms" << std::endl;
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
  
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n_rounds" << std::endl;
    return -1;
  }

  const int n_rounds = std::atoi(argv[1]);

  if (n_rounds < 0 || n_rounds > 100) {
    std::cerr << "WARNING: The number of rounds seems very high" << std::endl;
    return -1;
  }

  // Before starting the benchmark pin this thread to a core
  bool pinned = true;
  if (stick_this_thread_to_core(0) != 0) {
    std::cerr << "WARNING: cannot pin this thread to a core" << std::endl;
    pinned = false;
  }

  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  std::cout << "Real benchmark -> Rounds: " << n_rounds << ", Pinning: " << (pinned ? "yes" : "no") << ", host: " << hostname << std::endl;
  
  for (int round = 0; round < n_rounds; ++round) {
    std::string str_r = std::to_string(round) + "_";
    for (int i = start_unsigned; i < start_unsigned + datasets_unsigned.size(); ++i) {
      //benchmark<uint32_t>(str_r + "std::sort", i, [&](std::vector<uint32_t> && v) { std::sort(v.begin(), v.end()); });
      benchmark<uint32_t>(str_r + "IPS4O", i, [&](std::vector<uint32_t> && v) { ips4o::sort(v.begin(), v.end()); });
      //benchmark<uint32_t>(str_r + "SkaSort", i, [&](std::vector<uint32_t> && v) { ska_sort(v.begin(), v.end()); });
      benchmark<uint32_t>(str_r + "LearnedSort_2_0", i, [&](std::vector<uint32_t> && v) { learned_sort::sort(v.begin(), v.end()); });
      benchmark<uint32_t>(str_r + "ULS", i, [&](std::vector<uint32_t> && v) { learned_sort_framework::sort<indexes::MinMaxIndex<learned_cfg<uint32_t>>, indexes::MinMaxIndex<learned_cfg<uint32_t>>, 600, false, false>(v.begin(), v.end()); });
      benchmark<uint32_t>(str_r + "BLS", i, [&](std::vector<uint32_t> && v) { learned_sort_framework::sort<lbm_t<uint32_t, 4>, indexes::MinMaxIndex<learned_cfg<uint32_t>>, 600, false, false>(v.begin(), v.end()); });
      benchmark<uint32_t>(str_r + "LearnedSort_2_1", i, [&](std::vector<uint32_t> && v) { learned_sort_framework::sort<indexes::MyRmiIndex<learned_cfg<uint32_t>>, indexes::MinMaxIndex<learned_cfg<uint32_t>>, 600, false, false>(v.begin(), v.end()); });
      benchmark<uint32_t>(str_r + "RadixSort", i, [&](std::vector<uint32_t> && v) { radix_sort(v.begin(), v.end()); });
      benchmark<uint32_t>(str_r + "PDQS", i, [&](std::vector<uint32_t> && v) { pdqsort(v.begin(), v.end()); });
    }
  }

  return 0;
}
