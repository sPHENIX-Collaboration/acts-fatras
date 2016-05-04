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
    
    /** Static method : convert to ParticleType from pdg */
    static Acts::ParticleType convertToParticleType(int pdg, bool& stable, bool& exiting, double charge) {
    
        int pdgCode = abs(pdg);
    
        stable       = false;
        exiting      = false;
    
        Acts::ParticleType particleType;
    
      // try to follow number of appearance 
        switch (pdgCode )
        {
        // leptons
            case 11: // e+/e-
            particleType = Acts::electron;
            stable       = true;
            exiting      = false;
            break;
            case 13: // mu+/mu-
            particleType = Acts::muon;
            stable       = false;
            exiting      = false;
            break;
            case 12: // e neutrino
            case 14: // mu neutrino
            case 16: // tau neutrino
            particleType = Acts::nonInteracting;
            stable       = true;
            exiting      = true;
            break;
            case 22: // gamma
            particleType = Acts::photon;
            stable       = true;
            exiting      = false;
            break; 
            case 211: // pi+/pi-
            particleType = Acts::pion;
            stable       = false;
            exiting      = false;
            break;
            case 111: // pi0
            particleType = Acts::pi0;              
            stable       = false;
            exiting      = false;
            break;
            case 2212: // proton
            particleType = Acts::proton;
            stable       = true;
            exiting      = false;
            break;
            case 2112: // neutron               
            particleType = Acts::neutron;
            stable       = true;
            exiting      = true;
            break;
            case 321: // K
            particleType = Acts::kaon;
            stable       = false;
            exiting      = false;
            break;
            case 130: // K_long
            particleType = Acts::k0;
            stable       = false;
            exiting      = false;
            break;
            case 310: // K_short
            particleType = Acts::k0;
            stable       = false;
            exiting      = false;
            break;
            default: // treat mesons as pions
            particleType = charge != 0. ? Acts::pion : Acts::pi0 ;                               
            stable       = false;
            exiting      = false;
        }
    
      // and all baryons as proton hypo
        if (pdgCode > 999 && pdgCode!=2112 )
        {
            particleType = charge != 0. ? Acts::proton : Acts::neutron ;
            stable       = false;
            exiting      = false;
        }
    
      // ignore SUSY particles for now
        if (pdgCode > 1000000)
        {
            particleType = Acts::nonInteracting;
            stable       = true;
            exiting      = true;
        }
    
        return particleType;
    }
    
    
    /** Static method : convert to pdg from ParticleType */
    static int convertToPdg(Acts::ParticleType particleHypo, double charge, bool dist) 
    {
    
        int pdg = 0;
    
        switch (particleHypo) {
            // the electron case
            case Acts::electron   :  {  pdg = 11; pdg *= charge > 0. ? -1 : 1;   } return pdg;  
            // the muon case
            case Acts::muon       :  {  pdg = 13; pdg *= charge > 0. ? -1 : 1;   } return pdg;  
            // the kaon case
            case Acts::kaon       :  {  pdg = 321; pdg *= charge > 0. ? 1 : -1;  } return pdg;
            // the proton case
            case Acts::proton     :  {  pdg = 2212; pdg *= charge > 0. ? 1 : -1; 
                 if (charge*charge < 0.0001)
                     pdg = dist ? 2112 : -2112; } return pdg;
            // the photon case
            case Acts::photon     :  { pdg = 22; } return pdg;
            // the neutron case
            case Acts::neutron     :  { pdg = 2112; } return pdg;
            // the neutral pion case
            case Acts::pi0         :  { pdg = 111; } return pdg;
            // the neutral kaon case
            case Acts::k0          :  { pdg = dist ? 130 : 310;                      } return pdg;
            // the pion case - is the default
            default              :  {  pdg = 211; pdg *= charge > 0. ? 1 : -1; 
                if (charge*charge < 0.0001)
                    pdg = 111; };  
            }
           return pdg;
    }

}

#endif
