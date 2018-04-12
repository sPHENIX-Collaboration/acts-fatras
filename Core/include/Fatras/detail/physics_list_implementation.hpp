// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PHYSICS_LIST_IMPLEMENTATION_HPP
#define FATRAS_PHYSICS_LIST_IMPLEMENTATION_HPP

namespace Fatras {

namespace detail {

namespace {

  template <typename... processes>
  struct physics_list_impl;

  /// Recursive Call pattern
  /// - make sure it bails out 
  template <typename first, typename... others>
  struct physics_list_impl<first, others...>
  {
    template <typename T,
              typename cache_t,
              typename generator_t,
              typename detector_t,
              typename particle_t>
    static bool
    process(const T& process_tuple, 
           cache_t& cache,
           generator_t& generator,
           const detector_t& detector,
           const particle_t& ingoing,
           std::vector<particle_t>& outgoing)
    {
      // pick the first process
      const auto&  this_process = std::get<first>(process_tuple);
      bool this_process_kills = this_process(cache,
                                             generator,
                                             detector,
                                             ingoing,
                                             outgoing);
      // recursive call on the remaining ones
      return (this_process_kills 
              || physics_list_impl<others...>::process(process_tuple, 
                                                       cache, 
                                                       generator,
                                                       detector, 
                                                       ingoing, 
                                                       outgoing));
    }
  };

  /// Final call pattern
  template <typename last>
  struct physics_list_impl<last>
  {
    template <typename T,
              typename cache_t,
              typename generator_t,
              typename detector_t,
              typename particle_t>
    static bool
    process(const T& process_tuple, 
            cache_t& cache,
            generator_t& generator,
            const detector_t& detector,
            const particle_t& ingoing,
            std::vector<particle_t>& outgoing)
    {
      // this is the last process in the tuple
      const auto& this_process = std::get<last>(process_tuple);
      return this_process(cache,
                          generator,
                          detector,
                          ingoing,
                          outgoing);
    }
  };

  /// Empty call pattern
  template <>
  struct physics_list_impl<>
  {
    template <typename T,
              typename cache_t,
              typename generator_t,
              typename detector_t,
              typename particle_t>
    static bool
    process(const T&, 
            cache_t&,
            generator_t&,
            const detector_t&,
            const particle_t&,
            std::vector<particle_t>&)
    {
      return false;
    }
  };

}  // namespace

}  // namespace Acts
#endif  // FATRAS_PHYSICS_LIST_IMPLEMENTATION_HPP
