#pragma once

#include "bls/bls.h"
#include "bls/indexes/linear.h"
#include "bls/indexes/linear_balanced.h"
#include "bls/indexes/new_rmi.h"

namespace ls_framework {

    template <typename vec_iter>
    auto bls = sort<indexes::BalancedMinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 600, vec_iter>;

    template <typename vec_iter>
    auto uls = sort<indexes::MinMaxIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 1000, vec_iter>;

    template <typename vec_iter>
    auto ls21 = sort<indexes::NewRmiIndex<vec_iter>, indexes::MinMaxIndex<vec_iter>, 500, vec_iter>;

}