// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SELECTOR_SIGNATURE_CHECK_HPP
#define FATRAS_SELECTOR_SIGNATURE_CHECK_HPP

#include <type_traits>
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Fatras {

/// The following operator has to be inplemented in order to satisfy
/// as an sampler for fast simulation
///
/// @code
///  bool
///  operator()(const particle_t& particle) const { return true; }
///
/// @endcode
namespace detail {

  namespace {    
    template <typename T,
              typename particle_t,
              typename = decltype(std::declval<T>().
                                  operator()(std::declval<const particle_t&>()))>
                
    std::true_type
    test_selector_list(int);

    template <typename, typename>
    std::false_type
    test_selector_list(...);
    
    template <typename T, typename particle_t>
    struct selector_list_signature_check
        : decltype(test_selector_list<T,particle_t>(0))
    {
    };

  }  // end of anonymous namespace

  template <typename T, typename particle_t>
  constexpr bool selector_list_signature_check_v
      = selector_list_signature_check<T, particle_t>::value;
}  // namespace detail

}  // namespace Fatras

#endif  // FATRAS_SELECTOR_SIGNATURE_CHECK_HPP