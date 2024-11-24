/*
Two-layer RMI
*/

#pragma once

#include "linear.h"

#include <vector>

namespace indexes {

    template <typename RandomIt, typename count_type = uint, size_t n_leaves = 600>
    class MyRmiIndex {
        using value_type = typename std::iterator_traits<RandomIt>::value_type;

    protected:
        struct leaf_model {
            double slope;
            double intercept;
        };

        double root_slope, root_intercept;
        std::vector<leaf_model> leaf_models;
        const float sampling_offset;
        count_type max_index;

    public:
        MyRmiIndex(const float sampling_rate = 0.01) : 
            sampling_offset(1. / sampling_rate) {}
        
        inline bool build(const RandomIt begin, const RandomIt end, count_type num_buckets) {
            struct leaf_obj {
                size_t counter = 0;
                value_type min;
                value_type max;
            };
            
            // Set num_buckets
            max_index = num_buckets;
            leaf_models.resize(n_leaves);
            std::vector<leaf_obj> leafs(n_leaves);
            
            auto n = end - begin;
            const int offset = this->sampling_offset * 4 > n ? 1 : this->sampling_offset;
            auto n_samples = n / offset;
            
#ifndef RANDOM_SAMPLING
            if (offset > 1) {
                // Swap the samples at the begining of the sequence
                for (RandomIt current = std::next(begin), it = begin + offset; it < end; it+=offset, ++current) {
                    std::swap(*it, *current);
                }
            }
#else
            random_generator.seed(42); // helpful for the benchmarks
            selectSample(begin, end, n_samples);
#endif

            RandomIt min_iter, max_iter;
            std::tie(min_iter, max_iter) = std::minmax_element(begin, begin + n_samples);
            value_type min = *min_iter, max = *max_iter;

            double slope_, intercept_;
            // Check edge condition
            if (max == min) {
                root_slope = 1;
                root_intercept = 0;
                return false;
            } else {
                slope_ = 1. / (max - min);
                intercept_ = -slope_ * min;

                // Rescale w.r.t. the number of buckets
                slope_ *= n_leaves;
                intercept_ *= n_leaves;
            }

            root_slope = slope_;
            root_intercept = intercept_;

            // Init leafs
            for (int i = 0; i < n_leaves; ++i) {
                leafs[i].min = std::numeric_limits<value_type>::max();
                leafs[i].max = std::numeric_limits<value_type>::lowest();
            }
            
            // Scan the samples
            for (RandomIt it = begin; it < begin + n_samples; ++it) {
                // Cpmpute the root equation
                const value_type val = it[0];
                const count_type leaf_idx = static_cast<count_type>(std::max(0., std::min(slope_ * val + intercept_, n_leaves - 1.)));
                leafs[leaf_idx].min = std::min(leafs[leaf_idx].min, val);
                leafs[leaf_idx].max = std::max(leafs[leaf_idx].max, val);
                ++leafs[leaf_idx].counter;
            }

            double current_y = 0;
            value_type current_x;

            int i = 0;
            while(leafs[i].counter == 0) {
                leaf_models[i].slope = leaf_models[i].intercept = 0;
                ++i;
            }
            current_x = leafs[i].min;
            for (; i < n_leaves; ++i) {
                const double y_increment = 1. * leafs[i].counter / n_samples;
                if (leafs[i].counter < 1) {
                    // Cannot train the model with 0 points
                    leaf_models[i].slope = 0;
                    leaf_models[i].intercept = current_y;
                } else {
                    leaf_models[i].slope = y_increment / (leafs[i].max - current_x);
                    leaf_models[i].intercept = current_y - leaf_models[i].slope * current_x;
                }

                // Rescale
                leaf_models[i].slope *= num_buckets;
                leaf_models[i].intercept *= num_buckets;
                
                current_y += y_increment;
                current_x = leafs[i].max;
            }

            return true;
        }

        inline bool check_homogeneity() const { return true; }

        inline count_type classify(value_type value) const {
            // Compute the leaf index using the root model
            count_type leaf_idx = static_cast<count_type>(std::min(std::max(0., value * root_slope + root_intercept), n_leaves - 1.));

            // Compute the bucket prediction using the leaf model
            return static_cast<count_type>(std::min(std::max(0., value * leaf_models[leaf_idx].slope + leaf_models[leaf_idx].intercept), this->max_index - 1.));
        }
    };
}
