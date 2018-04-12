//
//  RandomNumbersDistributions.cpp
//  ACTFW
//
//  Created by Hadrien Grasland on 27/06/17.
//
//

#include "Fatras/RandomNumberDistributions.hpp"

Fatras::LandauDist::param_type::param_type(double mean, double scale)
  : mean(mean), scale(scale)
{
}

bool
Fatras::LandauDist::param_type::operator==(const param_type& other) const
{
  return (mean == other.mean) && (scale == other.scale);
}

Fatras::LandauDist::LandauDist(double mean, double scale) : m_cfg(mean, scale)
{
}

Fatras::LandauDist::LandauDist(const param_type& cfg) : m_cfg(cfg)
{
}

Fatras::LandauDist::result_type
Fatras::LandauDist::min() const
{
  return -std::numeric_limits<double>::infinity();
}

Fatras::LandauDist::result_type
Fatras::LandauDist::max() const
{
  return std::numeric_limits<double>::infinity();
}

bool
Fatras::LandauDist::operator==(const LandauDist& other) const
{
  return (m_cfg == other.m_cfg);
}
