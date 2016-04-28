///////////////////////////////////////////////////////////////////
// PdgToParticleHypothesis.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EVENTDATAUTILS_PDGTOPARTICLEHYPOTHESIS_H
#define ACTS_EVENTDATAUTILS_PDGTOPARTICLEHYPOTHESIS_H

#include "ACTS/EventData/ParticleHypothesis.h"
// STD/STL
#include <cmath>  

namespace Fatras {

    /** @class PdgToParticleHypothesis 

    small converter from the (abs) PDG code 
    to the particle hypothsis used in Tracking
    
    @author Andreas.Salzburger@cern.ch
    
    **/
    
    class PdgToParticleHypothesis {
    
      public:
        /** Constructor */
        PdgToParticleHypothesis(){}
    
        /** ~Destructor */
        ~PdgToParticleHypothesis(){}
    
        /** Converter method : PDG -> Particle Hyptothesis */
        Acts::ParticleHypothesis  convert(int pdg, bool& stable, bool& exiting, double charge=1.) const;
    
        /** Converter method : PDG -> Particle Hyptothesis , w/o stable exiting*/
        Acts::ParticleHypothesis  convert(int pdg, double charge=1.) const;
    
    
        /** Converter method : Particle Hyptothesis -> PDG*/
        int convert(Acts::ParticleHypothesis particleHypo, double charge, bool dist=true) const;
    
    
    };
    
    inline Acts::ParticleHypothesis PdgToParticleHypothesis::convert(int pdg, double charge) const {
        bool stable, exiting; 
        return convert(pdg,stable,exiting,charge); 
    }
    
    
    inline Acts::ParticleHypothesis PdgToParticleHypothesis::convert(int pdg, bool& stable, bool& exiting, double charge) const {
    
        int pdgCode = abs(pdg);
    
        stable       = false;
        exiting      = false;
    
        Acts::ParticleHypothesis particleType;
    
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
    
    
    inline int PdgToParticleHypothesis::convert(Acts::ParticleHypothesis particleHypo, double charge, bool dist) const 
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
                     pdg = dist ? 2112 : -2112;             } return pdg;
            // the photon case
            case Acts::photon     :  { pdg = 22;                                   } return pdg;
            // the neutron case
            case Acts::neutron     :  { pdg = 2112;                                } return pdg;
            // the neutral pion case
            case Acts::pi0         :  { pdg = 111;                      } return pdg;
            // the neutral kaon case
            case Acts::k0          :  { pdg = dist ? 130 : 310;                      } return pdg;
            // the pion case - is the default
            default              :  {  pdg = 211; pdg *= charge > 0. ? 1 : -1; 
                if (charge*charge < 0.0001)
                    pdg = 111;                             };  
            }
                 return pdg;
    }

}
#endif
