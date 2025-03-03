#pragma once

#include <iterator>
#include <random>

#define RANDOM_SAMPLING

namespace indexes {

    static constexpr bool is64bit = sizeof(size_t) == 8;

    static std::linear_congruential_engine<
            std::uintptr_t, is64bit ? 6364136223846793005u : 1664525u,
            is64bit ? 1442695040888963407u : 1013904223u, 0u>
            random_generator;

    /**
     * Selects a random sample in-place.
     * Adapted from https://github.com/ips4o/ips4o/blob/master/include/ips4o/sampling.hpp
     */
    template <class It>
    inline static void selectSample(It begin, const It end,
                    typename std::iterator_traits<It>::difference_type num_samples) {
        auto n = end - begin;
        while (num_samples--) {
            const auto i = std::uniform_int_distribution<
                    typename std::iterator_traits<It>::difference_type>(0, --n)(random_generator);
            std::swap(*begin, begin[i]);
            ++begin;
        }
    }

    template <typename RandomIt, typename count_type = uint>
    class LinearModel {
        using value_type = typename std::iterator_traits<RandomIt>::value_type;

    protected:
        const int sampling_offset;
        size_t max_index;
        double slope, intercept;

    public:
        LinearModel(const float sampling_rate = 0.01) :
            sampling_offset(1. / sampling_rate) {}

    protected:
        inline void setMaxIndex(size_t n_buckets) {
            // Set max_index
            max_index = static_cast<size_t>(n_buckets);
        }

        inline size_t sample(RandomIt begin, RandomIt end) const {
            auto n = end - begin;
            const int offset = this->sampling_offset * 4 > n ? 1 : this->sampling_offset;
            const auto n_samples = n / offset;

#ifndef RANDOM_SAMPLING
            if (offset > 1) {
                // Swap the samples at the begining of the sequence
                for (RandomIt current = std::next(begin), it = begin + offset; it < end; it+=offset, ++current) {
                    std::swap(*it, *current);
                }
            }
#else
            random_generator.seed(42); // helpful for the benchmarks repeatability
            selectSample(begin, end, n_samples);
#endif

            return n_samples;
        }
    
    public:
        inline bool check_homogeneity() const { return true; }

        inline count_type classify(value_type val) const {
            return static_cast<count_type>(std::max(0.,std::min(max_index - 1., slope * val + intercept)));
        }

        inline count_type classify_expanded(const value_type value, const count_type expansion_size, const count_type prev_index) const {
            return static_cast<count_type>(std::max(0., std::min((value * slope + intercept - prev_index) * expansion_size, expansion_size - 1.)));
        }
    };

    template <typename RandomIt, typename count_type = uint>
    class MinMaxIndex : public LinearModel<RandomIt, count_type> {
        using value_type = typename std::iterator_traits<RandomIt>::value_type;

    public:
        MinMaxIndex(const float sampling_rate = 0.01) : 
            LinearModel<RandomIt, count_type>(sampling_rate) {}

        bool build(const RandomIt begin, const RandomIt end, const count_type n_buckets) {
            
            LinearModel<RandomIt, count_type>::setMaxIndex(n_buckets);
            size_t n_samples = LinearModel<RandomIt, count_type>::sample(begin, end);

            RandomIt min_iter, max_iter;
            std::tie(min_iter, max_iter) = std::minmax_element(begin, begin + n_samples);
            value_type min = *min_iter, max = *max_iter;
            // Check edge condition
            if (max == min) {
                this->slope = 1;
                this->intercept = 0;
                return false;
            } else {
                this->slope = 1. / (max - min);
                this->intercept = -this->slope * min;

                // Rescale -> without this we should compute the prediction as [(slope * x + intercept) * scale_factor]
                const double scale_factor = static_cast<double>(this->max_index);
                this->slope *= scale_factor;
                this->intercept *= scale_factor;
                                
                return true;
            }
        }
    };
}
