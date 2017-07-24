#ifndef PAIRGREENSFUNCTION_HPP
#define PAIRGREENSFUNCTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT PairGreensFunction : public GreensFunction
{
public:
   PairGreensFunction(double D, double kf, double r0, double sigma) : GreensFunction(D), kf_(kf), r0_(r0), sigma_(sigma) {}

   double getkf() const { return kf_; }
   double getr0() const { return r0_; }
   double getSigma() const { return sigma_; }

   virtual double drawTheta(double d, double r, double t) const { ASSERT(false); UNUSED(d); UNUSED(r); UNUSED(t); /* NOT IMPLEMENTED */ return 0; }

protected:
   const double kf_;
   const double r0_;
   const double sigma_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* PAIRGREENSFUNCTION_HPP */