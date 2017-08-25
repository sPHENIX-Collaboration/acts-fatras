///////////////////////////////////////////////////////////////////
// IRandomNumbers.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_CORE_IRANDOMNUMBERSVC_H
#define ACTS_CORE_IRANDOMNUMBERSVC_H 1

namespace Fatras {
  
  /** @class IRandomNumberSvc
   * 
   * Random number service for drawing number with a specified distribution
   * 
   *  @author Noemi Calace -at- cern.ch
   */
  
  /** @enum Distribution
   * Distribution, enum for drawing number with a specified distribution
   */

  enum Distribution : unsigned int {
      Flat          = 0,
      Gauss         = 1,
      GaussZiggurat = 2,
      Landau        = 3,
      Gamma         = 4
  };
  
  class IRandomNumbers {
    
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////
  public:
    /** Virtual destructor */
    virtual ~IRandomNumbers() {}
      
    /** draw the random number 
     * with a specified distribution */
    virtual double draw(Fatras::Distribution, double k=1.0, double lambda=1.0) const = 0;
      
  };
}

#endif //> ! ACTS_CORE_IRANDOMNUMBERSVC_H
