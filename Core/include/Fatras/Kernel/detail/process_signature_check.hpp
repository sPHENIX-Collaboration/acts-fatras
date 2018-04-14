// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PROCESS_SIGNATURE_CHECK_HPP
#define FATRAS_PROCESS_SIGNATURE_CHECK_HPP

#include <type_traits>
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Fatras {

/// The following operator has to be inplemented in order to satisfy
/// as an sampler for fast simulation
///
/// @code
///  bool
///  operator()(cache_t& cache, 
///             generator_t& generator,
///             const detector_t& detector,
///             const particle_t& in,             
///             std::vector<particle_t>& out) const { return false; }
///
/// @endcode
namespace detail {

  namespace {    
    template <typename T,
              typename cache_t,
              typename generator_t,
              typename detector_t,
              typename particle_t,
              typename = decltype(std::declval<T>().
                                  operator()(std::declval<cache_t&>(),
                                             std::declval<generator_t&>(),
                                             std::declval<const detector_t&>(),
                                             std::declval<const particle_t&>(),
                                             std::declval<std::vector<particle_t>&>()))>
                
    std::true_type
    test_physics_list(int);

    template <typename, typename, typename, typename, typename>
    std::false_type
    test_physics_list(...);
    
    // clang-format on
    template <typename T, 
              typename cache_t, 
              typename generator_t, 
              typename detector_t,
              typename particle_t>
    struct physics_list_signature_check
        : decltype(test_physics_list<T, cache_t, 
                                        generator_t,
                                        detector_t,
                                        particle_t>(0))
    {
    };

    // clang-format on
  }  // end of anonymous namespace

  template <typename T, 
            typename cache_t, 
            typename generator_t,
            typename detector_t,
            typename particle_t>
  constexpr bool physics_list_signature_check_v
      = physics_list_signature_check<T, 
                                     cache_t, 
                                     generator_t,
                                     detector_t,
                                     particle_t>::value;
}  // namespace detail

}  // namespace Fatras

#endif  // FATRAS_PROCESS_SIGNATURE_CHECK_HPP