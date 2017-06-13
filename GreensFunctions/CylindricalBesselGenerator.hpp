#ifndef CYLINDRICALBESSELGENERATOR_HPP
#define CYLINDRICALBESSELGENERATOR_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <memory>
#include "DefsGf.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// for a given orders this table stores pre-calculated values (J and Jdot) or (Y and Ydot), from start with step delta, for N pairs
struct cbTable
{
   uint N;
   double start;
   double delta;
   std::unique_ptr<double[]> values;
};

// --------------------------------------------------------------------------------------------------------------------------------

class CylindricalBesselGenerator
{
public:
   CylindricalBesselGenerator();

   GF_EXPORT double J(uint n, double z) const;
   GF_EXPORT double Y(uint n, double z) const;

   uint MinNJ() const { return minNJ_; };
   uint MinNY() const { return minNY_; };
   uint MaxNJ() const { return maxNJ_; };
   uint MaxNY() const { return maxNY_; };

   GF_EXPORT static const CylindricalBesselGenerator& instance();

private:

   uint minNJ_;
   uint maxNJ_;
   uint minNY_;
   uint maxNY_;

   std::unique_ptr<cbTable[]> tableJ_;
   std::unique_ptr<cbTable[]> tableY_;

};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* CYLINDRICALBESSELGENERATOR_HPP */