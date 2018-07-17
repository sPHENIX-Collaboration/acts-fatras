// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Fatras/Kernel/Definitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"

namespace Fatras {

/// The struct for the physics list
///
///
struct ParametricNuclearInt {

struct Config()
{
	bool m_hadronInteractionFromX0;
	double m_hadronInteractionProbabilityScale;
}

/// @brief Constructor with given configuration
/// @param [in] cfg Configuration file
ParametricNuclearInt(Config& cfg);

/// @brief Calculates the hadronic with a given probability given by the properties of the material and the incoming particle
///
/// @tparam generator_t data type of the random number generator
/// @tparam material_t data type of the material
/// @tparam particle_t data type of the particle
/// @param generator random number generator
/// @param material penetrated material
/// @param pathCorrection scaling factor for the passed material
/// @param particle penetrating particle
template <typename generator_t, typename material_t, typename particle_t>
std::vector<particle_t> 
hadronicInteraction(generator_t& generator, const material_t& material, double pathCorrection, particle_t& particle) const;

  /// Call operator
  ///
  /// @tparam generator_t is a random number generator type
  /// @tparam detector_t is the detector information type
  /// @tparam particle_t is the particle information type
  ///
  /// @param[in] generator is the random number generator
  /// @param[in] detector the detector information
  /// @param[in] particle the particle which is being scattered
  ///
  /// @return eventually produced photons
  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &generator,
                                     const detector_t &detector,
                                     particle_t &particle) const;
                                     
private:

/// @brief Calculates the absorption length for various hadrons
///
/// @tparam material_t data type of the material
/// @tparam particle_t data type of the particle
/// @param [in] material material that is penetrated
/// @param [in] particle particle that penetrates the material
///
/// @return absorption length
template <typename material_t, typename particle_t>
double 
absorptionLength(const material_t* matertial, particle_t& particle) const;

// TODO
template<typename generator_t, typename particle_t>
int
diceNumberOfParticles(generator_t& generator, particle_t& particle) const;

/// @brief Calculates the hadron interactions of a particle
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] time starting time of the evolution
/// @param [in] particle particle that interacts
///
/// @return vector of outgoing particles
template<typename generator_t, typename particle_t>
std::vector<particle_t> 
getHadronState(generator_t& generator, double time, particle_t& particle) const;

// Configuration storage
Config m_cfg;
};

ParametricNuclearInt::ParametricNuclearInt(ParametricNuclearInt::Config& cfg) : m_cfg(cfg)
{
}

template <typename material_t, typename particle_t>
double 
ParametricNuclearInt::absorptionLength(const material_t* matertial, particle_t& particle) const 
{
  double al = material->l0();

  if(particle.pdg == 211 || particle.pdg == -211 || particle.pdg == 321 || particle.pdg == 111 || particle.pdg == 311) // TODO: nur K^+?
    al *= 1. / (1. + exp(-0.5 * (particle.p - 270.) * (particle.p - 270.) / 3600.)); // TODO: da kann man sicherlich noch etwas optimieren

  if(particle.pdg == 2212 || particle.pdg == 2212) al *= 0.7;
  if(particle.pdg == 211 || particle.pdg == -211 || particle.pdg == 111) al *= 0.9;

  return al;
}

template<typename generator_t, typename particle_t>
int
diceNumberOfParticles(generator_t& generator, particle_t& particle) const
{
  // sampling of hadronic interaction
  double E = sqrt(particle.p * particle.p + particle.m * particle.m);
  // get the maximum multiplicity    
  double multiplicity_max = 0.25 * E / 1000. + 18.;
  // multiplicity distribution
  double randx, randy, arg = 0.;
  double p1 = 0.;
  double p2 = 0.;
  if (E > 15000.) 
  {
    p1 = 8.69;
    p2 = 2.34;
  } 
  else 
  {
    p1 = 6.77;
    p2 = 2.02;
  }
  
  for (;;) {
    randx = 30. * generator();
    randy = 1.  * generator();
    arg = exp(-0.5 * ((randx - p1) / p2 + exp(-(randx - p1) / p2)));
    if (randy < arg && randx > 3 && randx < multiplicity_max) break;
  }
  
  randx *= (1.2 - 0.4 * exp(-0.001 * p));     // trying to adjust
  
  return (int)randx;
}

template<typename generator_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::getHadronState(generator_t& generator, particle_t& particle) const
{  
	std::vector<particle_t> chDef; 
  
  int Npart = diceNumberOfParticles(generator, particle);
  
  // protection against Npart < 3
  if (Npart < 3)
    return chDef;
  
  // create the genParticles
  
  // validation if needed 
	// TODO: bookkeeping via vertex oder einfach weglassen?
  
  if (m_cutChain && ( parent->barcode()>100000 || parent->barcode()==0 ) ) {
    if (m_hadIntValidationTree) m_hadIntValidationTree->Fill();
    return chDef;
  }
  
  int gen_part = 0; // TODO: wird erst am ende benoetigt
    
  // new sampling: sample particle type and energy in the CMS frame of outgoing particles
  // creation of shower particles
  double chargedist = 0.;
  std::vector<double> charge(Npart);
  std::vector<Trk::ParticleHypothesis> childType(Npart);
  std::vector<double> newm(Npart);
  std::vector<int> pdgid(Npart);    
  
  // children type sampling  : simplified
  //double pif = 0.19;
  //double nef = 0.20;
  //double prf = 0.20;

  // sample heavy particles (alpha) but don't save  
  double pif = 0.10; 
  double nef = 0.30;
  double prf = 0.30;
  
  if ( particle == Trk::pion || particle == Trk::kaon || particle == Trk::pi0 || particle == Trk::k0 ) {
      pif = 0.15;
      nef = 0.25;
      prf = 0.25;
    }
  if ( particle == Trk::proton ) {
    pif = 0.06;
    nef = 0.25;
    prf = 0.35;
  }
  if ( particle == Trk::neutron ) {
    pif = 0.03;
    nef = 0.35;
    prf = 0.17;
  }
  
  for (int i=0; i<Npart; i++) {
    chargedist  = CLHEP::RandFlat::shoot(m_randomEngine);
    if (chargedist<pif) {
      charge[i]=0.;
      childType[i]=Trk::pi0;
      newm[i]=s_particleMasses.mass[Trk::pi0]; // MeV
      pdgid[i]=111;
      continue;
    }
    if ( chargedist<2*pif) {
      charge[i]=1.;
      childType[i]=Trk::pion;
      newm[i]=s_particleMasses.mass[Trk::pion]; // MeV
      pdgid[i]=211;
      continue;
    }
    if (chargedist<3*pif) {
      charge[i]=-1.;
      childType[i]=Trk::pion;
      newm[i]=s_particleMasses.mass[Trk::pion]; // MeV
      pdgid[i]=-211;
      continue;
    }
    if (chargedist<3*pif+nef) {
      charge[i]=0.;
      childType[i]=Trk::neutron;
      newm[i]=939.565; // MeV
      pdgid[i]=2112; // neutron
      continue;
    }
    if (chargedist<3*pif+nef+prf) {
      charge[i]=1.;
      childType[i]=Trk::proton;
      newm[i]=s_particleMasses.mass[Trk::proton]; // MeV
      pdgid[i]=2212;
      continue;
    }
    charge[i]=2.;
    childType[i]=Trk::proton;
    newm[i]=4000.;
    pdgid[i]=20000;
  }

  // move the incoming particle type forward
  if ( childType[0] != particle ) {
    for (int i=1; i<Npart; i++) {
      if (childType[i]==particle) {
        childType[i]=childType[0];
        childType[0]=particle;
        double cho = charge[i];
        charge[i]=charge[0];
        charge[0]=parent ? parent->charge() : cho;
	newm[i]=s_particleMasses.mass[childType[i]]; // MeV
	newm[0]=s_particleMasses.mass[childType[0]]; // MeV
        break;
      }
    }
  }

//////////////////////////////// hier ist ein semantischer cut

  /*
  // sample momentum of daughters in CMS frame of outgoing particles  [ exp(-par/p) ]
  std::vector<double> mom;
  mom.clear();mom.reserve(Npart);
  double mom_n = 0.;  
  for (int _npart=0; _npart<Npart; _npart++) {
    rand1  = CLHEP::RandFlat::shoot(m_randomEngine);
    mom_n = -log(rand1)/m_childMomParam * p;
    int ipos = _npart;
    while ( ipos>0 && mom_n>mom[ipos-1]) ipos--;
    mom.insert(mom.begin()+ipos,mom_n);
  }  
  
  // check if configuration acceptable - if not, resample hardest mom
  double momR = 0.;
  for (int i=1; i<Npart; i++) momR += mom[i];
  if (momR < mom[0]) mom[0] = mom[1]+rand1*(momR-mom[1]);
  */

  std::vector<double> mom(Npart);
  std::vector<double> th(Npart);
  std::vector<double> ph(Npart);

  // sample first particle energy fraction and random momentum direction
  double eps = 2./Npart;
  double rnd  = CLHEP::RandFlat::shoot(m_randomEngine);
  mom[0] = 0.5*pow(eps,rnd);          
  th[0]  = acos( 2*CLHEP::RandFlat::shoot(m_randomEngine)-1.);
  ph[0]  = 2*M_PI*CLHEP::RandFlat::shoot(m_randomEngine);
  
  // toss particles around in a way which preserves the total momentum (0.,0.,0.) at this point
  // TODO shoot first particle along the impact direction preferentially

  Amg::Vector3D ptemp(mom[0]*sin(th[0])*cos(ph[0]),mom[0]*sin(th[0])*sin(ph[0]),mom[0]*cos(th[0]));
  double ptot = mom[0];
  
  double theta = 0.; double phi = 0.; 
  for (int i=1; i<Npart-2; i++) {
    eps = 1./(Npart-i); 
    //mom[i] = pow(eps,CLHEP::RandFlat::shoot(m_randomEngine))*(1-ptot);
    mom[i] = ( eps + CLHEP::RandFlat::shoot(m_randomEngine)*(1-eps))*(1-ptot); 
    if (ptemp.mag()<1-ptot) {
      while ( fabs(ptemp.mag()-mom[i])>1-ptot-mom[i] ){
	mom[i] =  ( eps + CLHEP::RandFlat::shoot(m_randomEngine)*(1-eps))*(1-ptot);      
      }
    }
    // max p remaining
    double p_rem=1-ptot-mom[i];
    double cthmax = fmin(1.,(-ptemp.mag()*ptemp.mag()-mom[i]*mom[i]+p_rem*p_rem)/2/ptemp.mag()/mom[i]);
    //if (cthmax<-1.) std::cout <<"problem in theta sampling:p_rem:ptot:pcurr:"<<p_rem<<","<<ptemp.mag()<<","<<mom[i]<< std::endl;
    double rnd  = CLHEP::RandFlat::shoot(m_randomEngine);
    theta = acos( (cthmax+1.)*rnd-1.);          
    phi = 2*M_PI*CLHEP::RandFlat::shoot(m_randomEngine);
    HepGeom::Vector3D<double> test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    HepGeom::Vector3D<double> dnewHep = HepGeom::RotateZ3D(ptemp.phi())*HepGeom::RotateY3D(ptemp.theta())*test;
    Amg::Vector3D dnew( dnewHep.x(), dnewHep.y(), dnewHep.z() );
    th[i]=dnew.theta();    
    ph[i]=dnew.phi();          
    ptemp += mom[i]*dnew;
    ptot += mom[i];
  }
  
  eps = 0.5; 
  mom[Npart-2] = pow(eps,CLHEP::RandFlat::shoot(m_randomEngine))*(1-ptot);
  mom[Npart-1] = 1-ptot-mom[Npart-2];
  
  if (ptemp.mag()<1-ptot) {
    while (mom[Npart-1]+mom[Npart-2]<ptemp.mag()) { 
      mom[Npart-2] = pow(eps,CLHEP::RandFlat::shoot(m_randomEngine))*(1-ptot);
      mom[Npart-1] = 1-ptot-mom[Npart-2];
    }
  }
  if (ptemp.mag()<fabs(mom[Npart-1]-mom[Npart-2]) ) {
    double diff = ptemp.mag()*CLHEP::RandFlat::shoot(m_randomEngine);
    double sum = mom[Npart-1]-mom[Npart-2];
    mom[Npart-2]=0.5*(sum+diff);  
    mom[Npart-1]=0.5*(sum-diff);  
  }
  double cth =(-ptemp.mag()*ptemp.mag()-mom[Npart-2]*mom[Npart-2]+mom[Npart-1]*mom[Npart-1])/2/ptemp.mag()/mom[Npart-2];
  if (fabs(cth)>1.) cth = (cth>0.) ? 1. : -1.;
  
  theta = acos(cth);
  phi = 2*M_PI*CLHEP::RandFlat::shoot(m_randomEngine);
  HepGeom::Vector3D<double> test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  HepGeom::Vector3D<double> dnewHep = HepGeom::RotateZ3D(ptemp.phi())*HepGeom::RotateY3D(ptemp.theta())*test;
  Amg::Vector3D dnew( dnewHep.x(), dnewHep.y(), dnewHep.z() );
  
  th[Npart-2]=dnew.theta();    
  ph[Npart-2]=dnew.phi();    
  ptemp += mom[Npart-2]*dnew;
  Amg::Vector3D dlast = -ptemp;
  th[Npart-1]=dlast.theta(); 
  ph[Npart-1]=dlast.phi();    
  
  // particle sampled, rotate, boost and save final state
  //CLHEP::HepLorentzVector bv(p*particleDir.unit().x(),p*particleDir.unit().y(),p*particleDir.unit().z(),s_currentGenParticle->momentum().e()+mtot);  
  double etot = 0.;
  for (int i=0;i<Npart; i++) etot += sqrt(mom[i]*mom[i]+newm[i]*newm[i]);
  double summ = 0.;
  for (int i=0;i<Npart; i++) summ += newm[i];

  // std::cout <<"hadronic interaction: current energy, expected :"<< etot <<","<< sqrt(summ*summ+2*summ*p+m*m)<< std::endl;
  // rescale (roughly) to the expected energy
  float scale = sqrt(summ*summ+2*summ*p+m*m)/etot;
  etot = 0.;
  for (int i=0;i<Npart; i++) {
    mom[i] *= scale;
    etot += sqrt(mom[i]*mom[i]+newm[i]*newm[i]);
  }

  
  CLHEP::HepLorentzVector bv(p*particleDir.unit().x(),p*particleDir.unit().y(),p*particleDir.unit().z(),sqrt(etot*etot+p*p));  
  std::vector<CLHEP::HepLorentzVector> childBoost(Npart);
  
  //std::cout <<"boost vector:"<<p<<","<<bv.e()<<","<<bv.m()<<std::endl;
  //std::cout <<"etot, mother E,m:"<<etot<<","<<E<<","<<m<<std::endl;
  
  Amg::Vector3D in(0.,0.,0.); 
  Amg::Vector3D fin(0.,0.,0.); 
  
  for (int i=0; i<Npart; i++) {
    Amg::Vector3D dirCms(sin(th[i])*cos(ph[i]),sin(th[i])*sin(ph[i]),cos(th[i])); 
    //Amg::Vector3D rotDirCms=HepGeom::RotateZ3D(particleDir.phi())*HepGeom::RotateY3D(particleDir.theta())*dirCms; 
    Amg::Vector3D childP = mom[i]*dirCms;
    in += childP;
    CLHEP::HepLorentzVector newp(childP.x(),childP.y(),childP.z(),sqrt(mom[i]*mom[i]+newm[i]*newm[i]));
    CLHEP::HepLorentzVector childPB = newp.boost(bv.boostVector());
    childBoost[i]=childPB;
    fin += Amg::Vector3D(childPB.x(),childPB.y(),childPB.z());
  } 
  
  double eout = 0.;
  
  // child particle vector for TruthIncident
  //  Reserve space for as many paricles as created due to hadr. int.
  //  However, the number of child particles for TruthIncident might be
  //  smaller due to (momentum) cuts
  ISF::ISFParticleVector           children(Npart);
  ISF::ISFParticleVector::iterator childrenIt = children.begin();
  unsigned short                numChildren = 0;
  
  for (int i=0; i<Npart; i++) {
    if (pdgid[i]<10000) {
      Amg::Vector3D childP = Amg::Vector3D(childBoost[i].x(),childBoost[i].y(),childBoost[i].z());
      Amg::Vector3D chP = Amg::Vector3D(sin(th[i])*cos(ph[i]),sin(th[i])*sin(ph[i]),cos(th[i]));
      
      eout += childBoost[i].e();     
      
      // validation if needed
      if (m_hadIntValidationTree && m_hadIntChildren < MAXHADINTCHILDREN){
	m_hadIntChildPdg[m_hadIntChildren]      = pdgid[i];   
	m_hadIntChildP[m_hadIntChildren]        = childP.mag();
	m_hadIntChildPcms[m_hadIntChildren]     = mom[i];
	m_hadIntChildTh[m_hadIntChildren]        = childP.unit().dot(particleDir);
	m_hadIntChildThc[m_hadIntChildren]       =chP.dot(particleDir);
	m_hadIntChildPhi[m_hadIntChildren]      = childP.phi();
	m_hadIntChildEta[m_hadIntChildren]      = childP.eta();
	double deltaPhi = m_hadIntMotherPhi - m_hadIntChildPhi[m_hadIntChildren];
	// rescale the deltaPhi
	deltaPhi -= deltaPhi > M_PI ? M_PI : 0.;
	deltaPhi += deltaPhi < -M_PI ? M_PI : 0.;		 
	m_hadIntChildDeltaPhi[m_hadIntChildren] = deltaPhi;
	m_hadIntChildDeltaEta[m_hadIntChildren] = m_hadIntMotherEta - m_hadIntChildEta[m_hadIntChildren];
	++m_hadIntChildren;
      }      
      
      if (childP.mag()> m_minimumHadOutEnergy) {
	// get the new particle    
	double mass = s_particleMasses.mass[ childType[i] ];
	
	// create the particle
	ISF::ISFParticle *child = new ISF::ISFParticle ( vertex,
							 childP,
							 mass,
							 charge[i],
							 pdgid[i],
							 time, 
							 *parent );
	// in the validation mode, add process info
	if (m_validationMode) {
	  ISF::ParticleUserInformation* validInfo = new ISF::ParticleUserInformation();
	  validInfo->setProcess(m_processCode);
	  if (parent->getUserInformation()) validInfo->setGeneration(parent->getUserInformation()->generation()+1);
	  else validInfo->setGeneration(1);     // assume parent is a primary track
	  child->setUserInformation(validInfo);
	}
	// record child for TruthIncident
	*childrenIt = child;
	++childrenIt; numChildren++;
      }
      
      gen_part++;
    }
  } // particle loop
  
  children.resize(numChildren);
  ISF::ISFTruthIncident truth( const_cast<ISF::ISFParticle&>(*parent),
			       children,
			       m_processCode,
			       parent->nextGeoID(),
			       ISF::fKillsPrimary );
  m_truthRecordSvc->registerTruthIncident( truth);

  // save info for validation
  if (m_validationMode && m_validationTool) {
    Amg::Vector3D* nMom = 0;
    m_validationTool->saveISFVertexInfo(m_processCode,vertex,*parent,p*particleDir,nMom,children);
    delete nMom;
  }

  
  m_hadIntChildE = eout;
  
  if (m_hadIntValidationTree) m_hadIntValidationTree->Fill();
  
  ATH_MSG_VERBOSE( "[ had ] it was kinematically possible to create " << gen_part << " shower particles " ); 
  
  return children;

}

template <typename generator_t, typename material_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::hadronicInteraction(generator_t& generator, const material_t& material, double pathCorrection, particle_t& particle) const
{
	const material_t* extMprop = dynamic_cast<const material_t*>(&material);
	double prob = 0.;

	// m_hadIntProbScale is used later, not here
	if (extMprop && !m_cfg.m_hadronInteractionFromX0) 
	{	  
		double al = absorptionLength(extMprop, particle);  // in mm
	
	    if (al > 0.) 
			prob = exp(-pathCorrection * extMprop->thickness() / al);
	    else
			prob = exp(-pathCorrection * extMprop->thicknessInL0());
	} 
	else 
		// using approximation lambda = 0.37 * Z * X0 instead -- giving a warning message
		prob = exp(-pathCorrection * mprop.thicknessInX0() / (0.37 * mprop.averageZ()));
  
	// apply a global scalor of the probability
	// (1. - prob) is generally O(0.01), so this is the right way to scale it
	if (generator() < (1. - prob) * m_hadronInteractionProbabilityScale) 
		return getHadronState(generator, particle);
  
	// no hadronic interactions were computed
	return {};  
}

} // namespace Fatras
