// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SELECTOR_LIST_IMPLEMENTATION_HPP
#define FATRAS_SELECTOR_LIST_IMPLEMENTATION_HPP

namespace Fatras {

namespace detail {

namespace {

  template <typename... selectors>
  struct selector_list_impl;

  /// Recursive call pattern
  /// - make sure it bails out 
  template <typename first, typename... others>
  struct selector_list_impl<first, others...>
  {
    template <typename T, typename particle_t>
    static bool
    select(const T& slector_tuple, 
           const particle_t& particle)
    {
      // pick the first select
      const auto&  this_selector = std::get<first>(slector_tuple);
      bool selected = this_selector(particle);
      // recursive call on the remaining ones, none of the selectors
      // is allowed to fail - a single failed selector rejects the particle
      // @todo check if rhs is actually evaluated if lhs fails (should not!)
      return (selected 
            && selector_list_impl<others...>::select(slector_tuple, particle));
    }
  };

  /// Final call pattern
  template <typename last>
  struct selector_list_impl<last>
  {
    template <typename T, typename particle_t>
    static bool
    select(const T& slector_tuple, const particle_t& particle)
    {
      // this is the last select in the tuple
      const auto& this_selector = std::get<last>(slector_tuple);
      return this_selector(particle);
    }
  };

  /// Empty call pattern
  template <>
  struct selector_list_impl<>
  {
    template <typename T, typename particle_t>
    static bool
    select(const T&, const particle_t& particle)
    {
      return true;
    }
  };

} // namespace

} // namespace detail

} // namespace Acts
#endif  // FATRAS_SELECTOR_LIST_IMPLEMENTATION_HPP
