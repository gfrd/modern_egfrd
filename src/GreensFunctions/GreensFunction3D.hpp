#ifndef GREENSFUNCTION3D_HPP
#define GREENSFUNCTION3D_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "PairGreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// Pair Green's function for the case where the pair never interact. kf == sigma == 0.
class GF_EXPORT GreensFunction3D : public PairGreensFunction
{
public:

   GreensFunction3D(double D, double r0) : PairGreensFunction(D, 0.0, r0, 0.0) { }

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction3D"; }

   double drawTime(double rnd) const override { UNUSED(rnd); return INFINITY; }

   double drawR(double rnd, double t) const override;

   double drawTheta(double rnd, double r, double t) const override;

private:

   friend class UnitTestGreenFunction3D;
   
   double p_r(double r, double t) const;

   double ip_r(double r, double t) const;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif // GREENSFUNCTION3D_HPP