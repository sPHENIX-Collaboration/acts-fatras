///////////////////////////////////////////////////////////////////
// FatrasDefinitions.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_FATRASHELPERS_H
#define ACTS_FATRAS_FATRASHELPERS_H 1

namespace Fatras {

    /** Static method: check for momentum cuts - checks on momentum and (optionally) transverse momentum */
    static bool passesMomentumCuts(const Acts::Vector3D& momentum, double cutM, double cutT)
    {
        // explicitely spell the code all cases to avoid unneccessary computations
        if (cutM < 0. && cutT < 0.) return true;
        // only check on transvese cut
        if (cutM < 0.) return momentum.perp2() > cutT*cutT;
        // only check the magnitude cut
        if (cutT < 0.) return momentum.mag2() > cutM*cutM;
        // full check
        return ( momentum.perp2() > cutT*cutT && momentum.mag2() > cutM*cutM );
    }   
    
}

#endif
