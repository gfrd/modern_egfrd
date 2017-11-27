#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction3DAbsSym : public GreensFunction
{
public:
   GreensFunction3DAbsSym(double D, double a) : GreensFunction(D), a_(a) {}

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction3DAbsSym"; }

   double geta() const { return a_; }

   double drawTime(double rnd) const override;

   double drawR(double rnd, double t) const override;

private:

   friend class UnitTestGreensFunction3DAbsSym; // for UnitTests

   double ellipticTheta4Zero(double q) const;

   double p_survival(double t) const;

   double p_int_r(double r, double t) const;

   double p_int_r_free(double r, double t) const;

   double p_r_fourier(double r, double t) const;

   const double a_;
};

// --------------------------------------------------------------------------------------------------------------------------------
