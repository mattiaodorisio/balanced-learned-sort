#pragma once

/**
 * This file has been adapted from: 
 * https://github.com/anikristo/LearnedSort/blob/master/src/generators.h
 * https://github.com/anikristo/LearnedSort/blob/master/src/utils.h
 * 
 * @file generators.h
 * @author Ani Kristo (anikristo@gmail.com)
 * @brief Synthetic data generator functions
 *
 * @copyright Copyright (c) 2020 Ani Kristo (anikristo@gmail.com)
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
#include <iterator>

using namespace std;

// An enumeration of data distribution types to generate for the benchmarks.
enum distr_t {
  CHI_SQUARED,
  EIGHT_DUPS,
  EXPONENTIAL,
  LOGNORMAL,
  MIX_GAUSS,
  MODULO,
  NORMAL,
  ROOT_DUPS,
  TWO_DUPS,
  UNIFORM,
  ZIPF,

  NONE
};
static inline std::string distr_to_string(distr_t d) {
  switch (d)
  {
  case CHI_SQUARED: return "CHI_SQUARED";                 
  case EIGHT_DUPS: return "EIGHT_DUPS";  
  case EXPONENTIAL: return "EXPONENTIAL";  
  case LOGNORMAL: return "LOGNORMAL";  
  case MIX_GAUSS: return "MIX_GAUSS";  
  case MODULO: return "MODULO";  
  case NORMAL: return "NORMAL";  
  case ROOT_DUPS: return "ROOT_DUPS";  
  case TWO_DUPS: return "TWO_DUPS";  
  case UNIFORM: return "UNIFORM";  
  case ZIPF: return "ZIPF";
  default: return "UNKNOWN";
  }
  return "UNKNOWN";
}

template <class ForwardIt>
void exponential_distr(ForwardIt first, const ForwardIt last, double lambda = 2, std::mt19937&& generator = std::mt19937(std::random_device{}())) {
  // Initialize random engine with normal distribution
  const size_t size = std::distance(first, last);
  exponential_distribution<> distribution(lambda);

  // Populate the input
  for (; first < last; ++first) {
    auto val = distribution(generator);
    if constexpr (std::is_integral_v<typename std::iterator_traits<ForwardIt>::value_type>) {
      val *= size; // Amplify to have some meaning with integers
    }
    *first = val;
  }
}

template <class ForwardIt>
void lognormal_distr(ForwardIt first, const ForwardIt last, double mean = 0, double stddev = 0.5,
                          double scale = 0, std::mt19937&& generator = std::mt19937(std::random_device{}())) {
  typedef typename std::iterator_traits<ForwardIt>::value_type T;
  const size_t size = std::distance(first, last);
  // Adjust the default scale parameter w.r.t. the numerical type
  if (!(is_same<float, T>() || is_same<double, T>()) && scale <= 0)
    scale = size;
  else if (scale <= 0)
    scale = 1;

  // Initialize random engine with normal distribution
  lognormal_distribution<> distribution(mean, stddev);

  // Populate the input
  for (; first < last; ++first) {
    *first = distribution(generator) * scale;
  }
}

template <class ForwardIt>
void modulo_distr(ForwardIt first, const ForwardIt last, size_t mod = 16) {
  // Populate the input
  const size_t size = std::distance(first, last);
  for (size_t i = 0; i < size; i++) {
    first[i] = i % mod;
  }
}

template <class ForwardIt>
void normal_distr(ForwardIt first, const ForwardIt last, double mean = 0, double stddev = 1, std::mt19937&& generator = std::mt19937(std::random_device{}())) {
  typedef typename std::iterator_traits<ForwardIt>::value_type T;
  // Initialize random engine with normal distribution
  normal_distribution<> distribution(mean, stddev);
  const size_t size = std::distance(first, last);

  // Populate the input
  for (; first < last; ++first) {
    auto val = distribution(generator);
    if constexpr (std::is_integral_v<T>) {
      val *= size; // Amplify to have some meaning with integers
      val += std::numeric_limits<T>::max() / 2; // To avoid underflow with unsigned types
    }
    *first = val;
  }
}

template <class ForwardIt>
void uniform_distr(ForwardIt first, const ForwardIt last, double a = 0, double b = -1, std::mt19937&& generator = std::mt19937(std::random_device{}())) {
  // Adjust the default parameters
  typedef typename std::iterator_traits<ForwardIt>::value_type T;
  const size_t size = std::distance(first, last);
  if (a == 0 && b == -1) {
    b = size;
    if (std::is_signed<T>::value) a = -1. * size;
  }

  // Initialize random engine with normal distribution
  uniform_real_distribution<> distribution(a, b);

  // Populate the input
  for (; first < last; ++first) {
    *first = distribution(generator);
  }
}

template <class ForwardIt>
void mix_of_gauss_distr(ForwardIt first, const ForwardIt last, size_t num_gauss = 5, std::mt19937&& generator = std::mt19937(std::random_device{}())) {
  // Generate the means
  const size_t size = std::distance(first, last);
  typedef typename std::iterator_traits<ForwardIt>::value_type T;

  vector<double> means(num_gauss);
  uniform_distr(means.begin(), means.end(), -500, 500, std::forward<std::mt19937>(generator));

  // Generate the stdevs
  vector<double> stdevs(num_gauss);
  uniform_distr(stdevs.begin(), stdevs.end(), 0, 100, std::forward<std::mt19937>(generator));

  // Generate the weights
  vector<double> weights(num_gauss);
  uniform_distr(weights.begin(), weights.end(), 0, 1, std::forward<std::mt19937>(generator));

  // Normalize the weights
  double sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.);
  std::for_each(weights.begin(), weights.end(),
                [&](auto x) { x /= sum_of_weights; });

  // Initialize random distribution selector
  discrete_distribution<int> index_selector(weights.begin(), weights.end());

  // Start generating random numbers from normal distributions
  for (; first < last; ++first) {
    auto random_idx = index_selector(generator);
    normal_distribution<> distribution(means[random_idx], stdevs[random_idx]);

    auto val = distribution(generator);
    if constexpr (std::is_integral_v<T>) {
      val *= size; // Amplify to have some meaning with integers
      val += std::numeric_limits<T>::max() / 2; // To avoid underflow with unsigned types
    }
    *first = val;
  }
}

template <class ForwardIt>
void chi_squared_distr(ForwardIt first, const ForwardIt last, double k = 4, std::mt19937&& generator = std::mt19937(std::random_device{}())) {
  typedef typename std::iterator_traits<ForwardIt>::value_type T;
  const size_t size = std::distance(first, last);
  // Initialize random number generator
  chi_squared_distribution<> distribution(k);

  // Populate the input
  for (; first < last; ++first) {
    auto val = distribution(generator);
    if constexpr (std::is_integral_v<T>) {
      val *= size; // Amplify to have some meaning with integers
    }
    *first = val;
  }
}

/**
 *
 * Adapted from
 * https://stackoverflow.com/questions/9983239/how-to-generate-zipf-distributed-numbers-efficiently
 *
 * NOTE: This generator is relatively slow, and the generation cannot be
 * parallelized due to static variables.
 */
template <class ForwardIt>
void zipf_distr(ForwardIt first_, const ForwardIt last_, double skew = 0.75,
                     size_t cardinality = 1e8) {
  typedef typename std::iterator_traits<ForwardIt>::value_type T;
  const size_t size = std::distance(first_, last_);

  // Start generating numbers
  for (size_t i = 0; i < size; ++i) {
    static bool first = true;  // Static first time flag
    static double c = 0;       // Normalization constant
    static double *sum_probs;  // Pre-calculated sum of probabilities

    // Compute normalization constant on first call only
    if (first) {
      for (size_t i = 1; i <= cardinality; ++i)
        c = c + (1.0 / pow((double)i, skew));
      c = 1.0 / c;

      sum_probs = (double *)malloc((cardinality + 1) * sizeof(*sum_probs));
      sum_probs[0] = 0;
      for (size_t i = 1; i <= cardinality; ++i) {
        sum_probs[i] = sum_probs[i - 1] + c / pow((double)i, skew);
      }
      first = false;
    }

    // Pull a uniform random number (0 < z < 1)
    double z;
    do {
      z = 1. * rand() / RAND_MAX;
    } while ((z == 0) || (z == 1));

    // Map z to the value
    size_t low = 1, high = cardinality, mid;
    T zipf_value = 0;  // Computed exponential value to be returned
    do {
      mid = floor((low + high) / 2);

      if (sum_probs[mid] >= z && sum_probs[mid - 1] < z) {
        zipf_value = mid;
        break;
      } else if (sum_probs[mid] >= z) {
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    } while (low <= high);

    first_[i] = zipf_value;
  }
}

/**
 * A[i] = i mod sqrt(N)
 */
template <class ForwardIt>
void root_dups_distr(ForwardIt first, const ForwardIt last) {
  // Populate the input
  const size_t size = std::distance(first, last);
  const size_t root = std::sqrt(size);

  for (size_t i = 0; i < size; i++) {
    first[i] = i % root;
  }
}

/**
 * Adapted from:
 * https://github.com/ips4o/ips4o-benchmark-suite/blob/4d11b0311833f59a31858314b5927101f6077389/src/generator/generator.hpp#L449
 *
 * A[i] = i^2 + N/2 mod N
 */
template <class ForwardIt>
void two_dups_distr(ForwardIt first, const ForwardIt last) {
  typedef typename std::iterator_traits<ForwardIt>::value_type T;
  const size_t size = std::distance(first, last);

  // Find the largest power of 2 that is less than size
  unsigned long largest_power_of_two = 1;
  while (2 * largest_power_of_two <= size) {
    largest_power_of_two *= 2;
  }

  // Limit maximum size of values for 32-bit keys
  if constexpr (std::is_same_v<T, unsigned>) {
    unsigned max = std::numeric_limits<T>::max();
    largest_power_of_two = std::min<unsigned long>(largest_power_of_two, max);
  }

  // Populate the input
  for (size_t i = 0; i < size; ++i) {
    first[i] = static_cast<T>((i * i + largest_power_of_two / 2) %
                            largest_power_of_two);
  }
}

/**
 * Adapted from:
 * https://github.com/ips4o/ips4o-benchmark-suite/blob/4d11b0311833f59a31858314b5927101f6077389/src/generator/generator.hpp#L488
 *
 * A[i] = i^8 + N/2 mod N
 */
template <class ForwardIt>
void eight_dups_distr(ForwardIt first, const ForwardIt last) {
  typedef typename std::iterator_traits<ForwardIt>::value_type T;
  const size_t size = std::distance(first, last);

  // Find the largest power of two
  unsigned long largest_power_of_two = 1;
  while (2 * largest_power_of_two <= size) {
    largest_power_of_two *= 2;
  }

  if (std::is_same<T, unsigned>())
    largest_power_of_two = std::min<unsigned long>(
        largest_power_of_two, std::numeric_limits<T>::max());

  // Populate the input
  for (size_t i = 0; i < size; ++i) {
    unsigned long temp = (i * i) % largest_power_of_two;
    temp = (temp * temp) % largest_power_of_two;
    first[i] = static_cast<T>((temp * temp + largest_power_of_two / 2) %
                            largest_power_of_two);
  }
}

// Utility to return the right distribution of the array
template <class ForwardIt>
void generate_data(distr_t data_distr, ForwardIt first, const ForwardIt last, uint seed = std::random_device{}()) {
  std::mt19937 generator(seed);
  switch (data_distr) {
    case CHI_SQUARED:
      return chi_squared_distr(first, last, 4, std::move(generator));
      break;

    case EIGHT_DUPS:
      return eight_dups_distr(first, last);
      break;

    case EXPONENTIAL:
      exponential_distr(first, last, 2, std::move(generator));
      break;

    case LOGNORMAL:
      return lognormal_distr(first, last, 0, 0.5, 0, std::move(generator));
      break;

    case MIX_GAUSS:
      return mix_of_gauss_distr(first, last, 5, std::move(generator));
      break;

    case MODULO:
      return modulo_distr(first, last);
      break;

    case NORMAL:
      return normal_distr(first, last, 0, 1, std::move(generator));
      break;

    case ROOT_DUPS:
      return root_dups_distr(first, last);
      break;

    case TWO_DUPS:
      return two_dups_distr(first, last);
      break;

    case UNIFORM:
      return uniform_distr(first, last, 0, 1, std::move(generator));
      break;

    case ZIPF:
      return zipf_distr(first, last);
      break;

    default:
      return normal_distr(first, last, 0, 1, std::move(generator));
      break;
  }

}
