// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAbortConditions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Fatras {

struct VoidDetector {};

/// @brief Fatras simulator
///
/// This is called from a Fatras steering algorithm
/// Per call, the generator is provided which
/// is then used to create a Acts propagator plugin.
///
/// @tparam charged_propagator_t Type of the propagator for charged particles
/// @tparam charged_selector_t Type of the slector (list) for charged particles
/// @tparam charged_interactor_t Type of the dresser for chargerd particles
///
/// @tparam neutral_propagator_t Type of the propagator for neutral particles
/// @tparam neutral_selector_t Type of the slector (list) for neutral particles
/// @tparam neutral_interactor_t Type of the dresser for neutral particles
template <typename charged_propagator_t, typename charged_selector_t,
          typename charged_interactor_t, typename neutral_propagator_t,
          typename neutral_selector_t, typename neutral_interactor_t>
struct Simulator {

  Simulator(charged_propagator_t chpropagator, neutral_propagator_t npropagator)
      : chargedPropagator(std::move(chpropagator)),
        neutralPropagator(std::move(npropagator)),
        mlogger(Acts::getDefaultLogger("Simulator", Acts::Logging::INFO)) {}

  charged_propagator_t chargedPropagator;
  charged_selector_t chargedSelector;
  charged_interactor_t chargedInteractor;

  neutral_propagator_t neutralPropagator;
  neutral_selector_t neutralSelector;
  neutral_interactor_t neutralInteractor;

  VoidDetector detector;

  std::shared_ptr<const Acts::Logger> mlogger = nullptr;

  bool debug = false;

  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *mlogger; }

  /// @brief call operator to the simulator
  ///
  /// @tparam generator_t Type of the generator object
  /// @tparam event_collection_t Type of the event collection
  /// @tparam hit_collection_t Type of the hit collection, needs insert()
  ///
  /// @param fatrasGenerator is the event-bound random generator
  /// @param fatrasEvent is the truth event collection
  /// @param fatrasHits is the hit collection
  template <typename generator_t, typename event_collection_t,
            typename hit_collection_t>
  void operator()(generator_t &fatrasGenerator, event_collection_t &fatrasEvent,
                  hit_collection_t &fatrasHits) const {

    // if screen output is required
    typedef Acts::detail::DebugOutputActor DebugOutput;

    // Action list, abort list and options
    typedef Acts::ActionList<charged_interactor_t, DebugOutput>
        ChargedActionList;
    typedef Acts::AbortList<Acts::detail::EndOfWorldReached> ChargedAbortList;
    typedef typename charged_propagator_t::template Options<ChargedActionList,
                                                            ChargedAbortList>
        ChargedOptions;

    // Action list, abort list and
    typedef Acts::ActionList<neutral_interactor_t, DebugOutput>
        NeutralActionList;
    typedef Acts::AbortList<Acts::detail::EndOfWorldReached> NeutralAbortList;
    typedef typename neutral_propagator_t::template Options<NeutralActionList,
                                                            NeutralAbortList>
        NeutralOptions;

    // loop over the input events
    // -> new secondaries will just be attached to that
    for (auto &event : fatrasEvent) {
      // take care here, the simulation can change the
      // particle collection
      for (auto particle = event.out.begin(); particle != event.out.end();
           ++particle) {
        // charged particle detected and selected
        if (chargedSelector(detector, *particle)) {
          // Need to construct them per call to set the particle
          // Options and configuration
          ChargedOptions chargedOptions;
          chargedOptions.debug = debug;
          // get the charged interactor
          auto &chargedInteractor =
              chargedOptions.actionList.template get<charged_interactor_t>();
          // result type typedef
          typedef typename charged_interactor_t::result_type ChargedResult;
          // set the generator to guarantee event consistent entires
          chargedInteractor.generator = &fatrasGenerator;
          // put all the additional information into the interactor
          chargedInteractor.initialParticle = (*particle);
          // create the kinematic start parameters
          Acts::CurvilinearParameters start(nullptr, particle->position,
                                            particle->momentum, particle->q);
          // run the simulation
          const auto &result =
              chargedPropagator.propagate(start, chargedOptions);
          auto &fatrasResult = result.template get<ChargedResult>();
          // a) handle the hits
          // hits go to the hit collection, particle go to the particle
          // collection
          for (auto &fHit : fatrasResult.simulatedHits) {
            fatrasHits.insert(fHit);
          }
          // b) deal with the particles
          const auto &simparticles = fatrasResult.outgoing;
          event.out.insert(particle, simparticles.begin(), simparticles.end());
          // c) screen output if requested
          if (debug) {
            auto &fatrasDebug = result.template get<DebugOutput::result_type>();
            ACTS_INFO(fatrasDebug.debugString);
          }
        } else if (neutralSelector(detector, *particle)) {
          // Options and configuration
          NeutralOptions neutralOptions;
          neutralOptions.debug = debug;
          // get the charged interactor
          auto &neutralInteractor =
              neutralOptions.actionList.template get<neutral_interactor_t>();
          // result type typedef
          typedef typename neutral_interactor_t::result_type NeutralResult;
          // set the generator to guarantee event consistent entires
          neutralInteractor.generator = &fatrasGenerator;
          // put all the additional information into the interactor
          neutralInteractor.initialParticle = (*particle);
          // create the kinematic start parameters
          Acts::NeutralCurvilinearParameters start(nullptr, particle->position,
                                                   particle->momentum);
          const auto &result =
              neutralPropagator.propagate(start, neutralOptions);
          auto &fatrasResult = result.template get<NeutralResult>();
          // a) deal with the particles
          const auto &simparticles = fatrasResult.outgoing;
          event.out.insert(particle, simparticles.begin(), simparticles.end());
          // b) screen output if requested
          if (debug) {
            auto &fatrasDebug = result.template get<DebugOutput::result_type>();
            ACTS_INFO(fatrasDebug.debugString);
          }
        } // neutral processing
      }   // loop over particles
    }     // loop over events
  }
};

} // namespace Fatras