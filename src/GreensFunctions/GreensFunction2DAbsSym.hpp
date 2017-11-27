#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction2DAbsSym : public GreensFunction
{
public:

   GreensFunction2DAbsSym(double D, double a) : GreensFunction(D), a_(a) {}

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction2DAbsSym"; }

   double geta() const { return a_; }

   double drawTime(double rnd) const override;

   double drawR(double rnd, double t) const override;

   double p_survival(double t) const;
   double p_int_r(double r, double t) const;
   double p_int_r_free(double r, double t) const;

private:

   struct p_survival_params
   {
      const GreensFunction2DAbsSym* const gf;
      double rnd;
   };

   struct p_r_params
   {
      const GreensFunction2DAbsSym* const gf;
      double t;
      double target;
   };

   static double p_survival_F(double t, const p_survival_params* params);

   static double p_r_free_F(double r, const p_r_params* params);

   static double p_r_F(double r, const p_r_params* params);

   const double a_;
};

// --------------------------------------------------------------------------------------------------------------------------------
