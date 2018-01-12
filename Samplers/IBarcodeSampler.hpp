//////////////////////////////////////////////////////////////////
// IBarcodeSampler.h, Acts project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_IBARCODESAMPLER_H
#ifndef ACTS_FATRAS_IBARCODESAMPLER_H

namespace Fatras {

/// @class IBarcodeSampler
///
/// Interface class for Barcode creation in Fatras
///
/// @author Andreas.Salzburger@cern.ch

class IBarcodeSampler
{
public:
  /// Virtual destructor
  virtual ~IBarcodeSampler() = default;

  /// Interface for processing of the presampled conversion on layer
  /// @param
  /// @todo there is something wrong
  virtual barcode_type
  generateChildBarcode(barcode_type mother, ) const;
};
}

#endif /* ACTS_FATRAS_IBARCODESAMPLER_H */
