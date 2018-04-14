// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PHYSICS_LIST_HPP
#define FATRAS_PHYSICS_LIST_HPP

#include "ACTS/Utilities/detail/Extendable.hpp"
#include "ACTS/Utilities/detail/MPL/all_of.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"
#include "Fatras/detail/physics_list_implementation.hpp"
#include "Fatras/detail/physics_list_signature_check.hpp"

namespace Fatras {

/// This is the PhysicsList struct that is used for fast simulation
/// Users can add a variable list of processes in order to drive the
/// physics simulation
template <typename... processes>
struct PhysicsList : private Acts::detail::Extendable<processes...>
{
private:
  static_assert(not Acts::detail::has_duplicates_v<processes...>,
                "same action type specified several times");

  using Acts::detail::Extendable<processes...>::tuple;

public:
  
  using Acts::detail::Extendable<processes...>::get;

  /// Call operator that is that broadcasts the call to the tuple()
  /// members of the list
  ///
  /// @tparam cache_t is the cache type from the propagator 
  /// @tparam generator_t is the random number generator type
  /// @tparam detector_t is the detector information type used
  /// @tparam particle_t is the particle type used in simulation
  ///
  /// @param[in] cache is the propgator/stepper cache
  /// @param[in] cache is the generator
  /// @param[in] detector is the necessary detector information
  /// @param[in] ingoig is the ingoing particle
  /// @param[in,out] outgoing are the (eventually) outgoing particles
  ///
  /// @return inticator which would trigger an abort
  template <typename cache_t, 
            typename generator_t,
            typename detector_t, 
            typename particle_t>
  bool
  operator()(cache_t& cache, 
             generator_t& gen,
             const detector_t& det,
             const particle_t& in,             
             std::vector<particle_t>& out) const
  {
    // clang-format off
    static_assert(Acts::detail::all_of_v<detail::physics_list_signature_check_v<processes, cache_t, generator_t, detector_t, particle_t>...>,
                  "not all process processes support the specified interface");
    // clang-format on
    
    // create an emtpy particle vector            
    typedef detail::physics_list_impl<processes...> impl;
    return impl::process(tuple(),cache,gen,det,in,out);
  }
  
};

}  // namespace Fatras

#endif  // FATRAS_PHYSICS_LIST_HPP
