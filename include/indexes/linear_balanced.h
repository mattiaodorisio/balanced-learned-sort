/*
This functions compute the optimal balanced P partitions for an array of weights.
It's inspired from:
https://cs.stackexchange.com/questions/125821/how-can-i-improve-my-algorithm-for-finding-optimally-balanced-p-way-partitioning
and the related python code:
https://replit.com/@Apass_Jack/optimallybalancedpartitioncs125821#main.py
*/

#pragma once

#include "linear.h"

namespace indexes {

    template <class RandomIt, size_t oversampling_factor = 4, typename count_type = uint>
    class BalancedMinMaxIndex : public MinMaxIndex<RandomIt, count_type> {
        using value_type = typename std::iterator_traits<RandomIt>::value_type;
    
    public:
        BalancedMinMaxIndex(const float sampling_rate = 0.01) : 
            MinMaxIndex<RandomIt, count_type>(sampling_rate) {}

        bool build(const RandomIt begin, const RandomIt end, count_type num_buckets) {
            // Set num_buckets
            this->max_index = static_cast<size_t>(num_buckets) * oversampling_factor;
            mapping_table.resize(num_buckets * oversampling_factor);

            auto n = end - begin;
            const int offset = this->sampling_offset * 4 > n ? 1 : this->sampling_offset;
            const auto n_samples = n / offset;

            // Swap the samples at the begining of the sequence
#ifndef RANDOM_SAMPLING
            for (RandomIt current = std::next(begin), it = begin + offset; it < end; it+=offset, ++current) {
                std::swap(*it, *current);
            }
#else
            random_generator.seed(42); // helpful for the benchmarks
            selectSample(begin, end, n_samples);
#endif

            RandomIt sample_begin = begin;
            RandomIt sample_end = begin + n_samples;

            RandomIt min_iter, max_iter;
            std::tie(min_iter, max_iter) = std::minmax_element(sample_begin, sample_end);
            value_type min = *min_iter, max = *max_iter;

            // Check edge condition
            if (min == max) {
                // Pick another samples by scanning a sequence linearly
                sample_begin = begin + n_samples;
                sample_end = begin + 2 * n_samples;

                std::tie(min_iter, max_iter) = std::minmax_element(sample_begin, sample_end);
                min = *min_iter, max = *max_iter;
            }

            // Check the condition again
            if (max == min) {
                this->slope = 1;
                this->intercept = 0;
                return false;
            } else {
                this->slope = 1. / (max - min);
                this->intercept = -this->slope * min;

                // Rescale -> without this we should compute the prediction as [(slope * x + intercept) * scale_factor]
                const double scale_factor = static_cast<double>(num_buckets * oversampling_factor);
                this->slope *= scale_factor;
                this->intercept *= scale_factor;
            }

            // Build the mapping table
            compute_balanced_partition(sample_begin, sample_end, mapping_table.data());

            return true;
        }

        inline bool check_homogeneity() const { return true; }

        std::vector<count_type> mapping_table;

        inline count_type classify(value_type value) const {
            return mapping_table[static_cast<count_type>(std::min(std::max(0., value * this->slope + this->intercept), this->max_index - 1.))];
        }

    private:
        static count_type count_least_partition(const size_t * weights, const size_t n, const size_t w, const count_type threshold) {
            /* Count the number of partitions up to a threshold, returns threshold + 1 
            if the number of partition is greater than the threshold */

            size_t i = 0;
            count_type n_partitions = 0;

            while (i < n && n_partitions <= threshold) {
                size_t seg_total = weights[i];
                i += 1;
                while (i < n) {
                    if (seg_total + weights[i] <= w) {
                        seg_total += weights[i];
                        i += 1;
                    }
                    else
                        break;
                }

                ++n_partitions;
            }

            return n_partitions;
        }

        static count_type least_partition(const size_t * weights, const size_t n, const size_t w, size_t * partitions) {
            size_t i = 0;
            count_type part_idx = 0;
            while (i < n) {
                size_t seg_total = weights[i];
                i += 1;
                while (i < n) {
                    if (seg_total + weights[i] <= w) {
                        seg_total += weights[i];
                        i += 1;
                    }
                    else
                        break;
                }   

                partitions[part_idx++] = i;
            }

            return part_idx;
        }

        count_type optimally_balanced_partition(const size_t * weights, const size_t n, const count_type P, size_t * partitions, size_t sum = 0) const {
            // Edge conditions
            if (P <= 1) {
                partitions[0] = n;
                return 1;
            } else if (P >= n) {
                for (size_t i = 0; i < n; ++i)
                    partitions[i] = i + 1;
                return n;
            }

            if (sum <= 0)
                sum = std::accumulate(weights, weights + n, 0);
            size_t low = 1;
            size_t high = sum;

            // First iteartion, choose mid:
            {
                size_t mid = 4 * high / P;
                count_type n_part = count_least_partition(weights, n, mid, P);
                if (n_part <= P)
                    high = mid;
                else if (n_part > P)
                    low = mid;
            }

            while (low < high - 1) {
                size_t mid = (low + high) / 2;
                count_type n_part = count_least_partition(weights, n, mid, P);
                if (n_part < P)
                    high = mid;
                else if (n_part > P)
                    low = mid;
                else {
                    high = mid;
                    break;
                }
            }

            count_type n_part = least_partition(weights, n, high, partitions);
            return n_part;
        }

        inline void compute_balanced_partition(RandomIt begin, RandomIt end, count_type *map_result) {
            const size_t expanded_size = this->max_index;
            const size_t reduced_size = expanded_size / oversampling_factor;

            const size_t n_samples = std::distance(begin, end);
            if (n_samples <= reduced_size) {
                for (size_t i = 0; i < reduced_size; ++i)
                    for (size_t j = 0; j < oversampling_factor; ++j)
                        map_result[i*oversampling_factor+j] = i;
                return;
            }
            
            size_t cnt_hist[expanded_size]{0};
            for (auto it = begin; it < begin + n_samples; ++it) {
                count_type pred = MinMaxIndex<RandomIt, count_type>::classify(*it);
                ++cnt_hist[pred];
            }

            size_t partitions[reduced_size];
            count_type n_part = optimally_balanced_partition(cnt_hist, expanded_size, reduced_size, partitions, n_samples);

            size_t prev_partition = 0;
            for (count_type i = 0; i < n_part; ++i) {
                for (count_type j = prev_partition; j < partitions[i]; ++j)
                    map_result[j] = i;
                prev_partition = partitions[i];
            }
        }
    };
}
