// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SELECTOR_LIST_HPP
#define FATRAS_SELECTOR_LIST_HPP

#include "ACTS/Utilities/detail/Extendable.hpp"
#include "ACTS/Utilities/detail/MPL/all_of.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"
#include "Fatras/Kernel/detail/selector_list_implementation.hpp"
#include "Fatras/Kernel/detail/selector_signature_check.hpp"

namespace Fatras {

/// This is the PhysicsList struct that is used for fast simulation
/// Users can add a variable list of selectors in order to drive the
/// physics simulation
template <typename... selectors>
struct SelectorList : private Acts::detail::Extendable<selectors...>
{
private:
  static_assert(not Acts::detail::has_duplicates_v<selectors...>,
                "same action type specified several times");

  using Acts::detail::Extendable<selectors...>::tuple;

public:
  
  using Acts::detail::Extendable<selectors...>::get;

  /// Call operator that is that broadcasts the call to the tuple()
  ///
  /// @tparam particle_t is the particle type used in simulation
  ///
  /// @param[in] particle to be checked for further processing
  ///
  /// @return inticator if the particle is accepted
  template <typename particle_t>
  bool
  operator()(const particle_t& particle) const
  {
    // clang-format off
    static_assert(Acts::detail::all_of_v<detail::selector_list_signature_check_v<selectors, particle_t>...>,
                  "not all particle selectors support the specified interface");
    // clang-format on
    
    // create an emtpy particle vector            
    typedef detail::selector_list_impl<selectors...> impl;
    return impl::select(tuple(),particle);
  }
  
};

}  // namespace Fatras

#endif  // FATRAS_SELECTOR_LIST_HPP