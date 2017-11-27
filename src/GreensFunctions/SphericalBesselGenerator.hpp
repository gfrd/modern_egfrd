#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <memory>
#include "DefsGf.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// for a given orders this table stores pre-calculated values (j and jdot) or (y and ydot), from start with step delta, for N pairs
struct sbTable
{
   uint N;
   double start;
   double delta;
   std::unique_ptr<double[]> values;
};

// --------------------------------------------------------------------------------------------------------------------------------

class SphericalBesselGenerator
{
public:
   SphericalBesselGenerator();     // constructor loads table from disk (on first use)!

   GF_EXPORT double j(uint n, double z) const;
   GF_EXPORT double y(uint n, double z) const;

   uint MinNj() const { return minNj_; };
   uint MinNy() const { return minNy_; };
   uint MaxNj() const { return maxNj_; };
   uint MaxNy() const { return maxNy_; };

   GF_EXPORT static const SphericalBesselGenerator& instance();

private:

   uint minNj_;
   uint maxNj_;
   uint minNy_;
   uint maxNy_;

   std::unique_ptr<sbTable[]> tableJ_;
   std::unique_ptr<sbTable[]> tableY_;
};

// --------------------------------------------------------------------------------------------------------------------------------
