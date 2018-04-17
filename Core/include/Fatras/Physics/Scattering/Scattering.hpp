// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PROCESSES_SCATTERING_HPP
#define FATRAS_PROCESSES_SCATTERING_HPP

#include <cmath>

#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Fatras {

  /// @class Scatterering
  ///
  /// This is the (multiple) scattering plugin to the
  /// Physics list. It needs a scattering formula in order
  /// to provide the the scattering angle in 3D space.
  /// 
  /// There's two options to apply the scattering
  /// - a parametric action that relates phi and theta (default: off)
  /// - an actuall out of direction scattering applying two random numbers
    
  template <typename Formula>  
  struct Scattering
  {
  
    /// Include the log term
    bool     parametric = false;
    double   projectionFactor = 1./std::sqrt(2.);
    
    /// The scattering formula
    Formula         angle;
          
    /// This is the scattering call operator
    template <typename generator_t,
              typename detector_t,
              typename particle_t>
    std::vector<particle_t>
    operator()(generator_t& gen,
               const detector_t& det,
               particle_t& in) const 
    { 
      /// uniform distribution
      UniformDist uniformDist = UniformDist(0., 1.);
      
      // 3D scattering angle
      double angle3D = angle(gen,det,in);
      
      // parametric scattering
      if (parametric){

        // the initial values
        double theta    = in.momentum.theta();
        double phi      = in.momentum.phi();
        double sinTheta = (sin(theta) * sin(theta) > 10e-10) ? sin(theta) : 1.;

        // sample them in an independent way
        double deltaTheta = projectionFactor * angle3D;
        double numDetlaPhi = 0. ; //??
        double deltaPhi   = projectionFactor *  numDetlaPhi/ sinTheta;

        // use bound parameter 
        // (i) phi
        phi += deltaPhi;
        if (phi >= M_PI) phi -= M_PI;
        else if (phi < -M_PI)  phi += M_PI;
        // (ii) theta
        theta += deltaTheta;
        if (theta > M_PI) theta -= M_PI;
        else if (theta < 0.) theta += M_PI;
        
        double sphi   = std::sin(phi);
        double cphi   = std::cos(phi);
        double stheta = std::sin(theta);
        double ctheta = std::cos(theta);

        // assign the new values
        in.momentum = in.p * Acts::Vector3D(cphi*stheta,sphi*stheta,ctheta);
        
      } else {
        // Create a random uniform distribution between in the intervall [0,1]      
        double psi     = 2. * M_PI * uniformDist(gen);
      
        // more complex but "more true"
        Acts::Vector3D pDirection(in.momentum.unit());
        double         x = -pDirection.y();
        double         y =  pDirection.x();
        double         z = 0.;
        
        // if it runs along the z axis - no good ==> take the x axis
        if (pDirection.z() * pDirection.z() > 0.999999) 
          { x = 1.;  y = 0.; }
        // deflector direction
        Acts::Vector3D deflector(x, y, z);
        // rotate the new direction for scattering using theta and  psi
        Acts::RotationMatrix3D rotation;
        rotation = Acts::AngleAxis3D(angle3D, deflector)
            * Acts::AngleAxis3D(psi, pDirection);
        // rotate and set a new direction to the cache 
        in.momentum = in.p * rotation*pDirection.unit();
      }
      // scattering always returns an empty list
      // - it is a non-distructive process
      return {}; 
    }
  
  };
  
}  // end of namespace

#endif  // FATRAS_PROCESSES_SCATTERING_HPP
