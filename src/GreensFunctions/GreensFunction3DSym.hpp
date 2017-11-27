#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// Green's Function for a free diffusion particle.
class GF_EXPORT GreensFunction3DSym : public GreensFunction
{
public:

   GreensFunction3DSym(double D) : GreensFunction(D) { }

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction3DSym"; }

   double drawTime(double rnd) const override { UNUSED(rnd); return INFINITY; }

   double drawR(double rnd, double t) const override;

private:

   double p_r(double r, double t) const;

   double ip_r(double r, double t) const;
};

// --------------------------------------------------------------------------------------------------------------------------------
