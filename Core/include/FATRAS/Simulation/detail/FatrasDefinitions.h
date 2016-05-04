///////////////////////////////////////////////////////////////////
// FatrasDefinitions.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_FATRASDEFINITIONS_H
#define ACTS_FATRAS_FATRASDEFINITIONS_H

#include "ACTS/EventData/ParticleDefinitions.h"
#include "ACTS/Extrapolation/detail/MaterialInteraction.h"
// STD/STL
#include <cmath>  

namespace Fatras {
    
    
    /** the formulas for energy loss evaluation */
    static const Acts::MaterialInteraction s_interactionFormulae;     
    
    /** struct of Particle masses */
    static const Acts::ParticleMasses      s_particleMasses; 
        
    /** KOverA factor in Bethe-Bloch equation [MeV*cm2/gram] */
    static const double                    s_ka_BetheBloch        = 30.7075;
    
    /** Fine structure constant */
    static const double                    s_alpha                = 1./137.;

    /** Multiple scattering paramters */
    static const double                    s_main_RutherfordScott = 13.6;
    static const double                    s_log_RutherfordScott  =  0.038;

    static const double                    s_main_RossiGreisen    = 17.5;
    static const double                    s_log_RossiGreisen     =  0.125;
    
    // statics doubles used for calculations
    static const double                    s_sqrtTwo              =  sqrt(2.);
    static const double                    s_oneOverThree         = 1./3.;
    
    
}

#endif
