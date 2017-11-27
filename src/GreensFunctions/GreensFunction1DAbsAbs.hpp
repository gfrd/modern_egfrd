#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction1DAbsAbs : public GreensFunction
{
public:
   GreensFunction1DAbsAbs(double D, double r0, double sigma, double a, double v = 0.0)
      : GreensFunction(D), v_(v), sigma_(sigma), a_(a), r0_(r0), lscale_(GfCfg.L_TYPICAL_1D), tscale_(GfCfg.T_TYPICAL_1D) { }

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction1DAbsAbs"; }

   double getsigma() const { return sigma_; }

   double geta() const { return a_; }

   double getv() const { return v_; }

   double getr0() const { return r0_; }

   // Draws the first passage time from the propensity function
   double drawTime(double rnd) const override;

   // Draws the position of the particle at a given time
   double drawR(double rnd, double t) const override;

   // Determines based on the flux ratios if the particle left the left or right boundary
   EventKind drawEventType(double rnd, double t) const override;

private:

   // Calculates the amount of flux leaving the left boundary at time t
   double leaves(double t) const;

   // Calculates the amount of flux leaving the right boundary at time t
   double leavea(double t) const;

   // Calculates the probability of finding the particle inside the 
   // domain at time t so, the survival probability
   double p_survival(double t) const;

   // Calculates the probability density of finding the particle at location r at time t.
   // double prob_r(double r, double t) const;

   uint guess_maxi(double t) const;
   
   double p_survival_table(double  t, DoubleVector& psurvTable) const;
   
   double p_survival_i(uint i, double t, const DoubleVector& table) const;

   double p_survival_table_i_v(uint i) const;

   double p_survival_table_i_nov(uint i) const;

   void createPsurvTable(uint maxi, DoubleVector& table) const;

   double p_int_r_table(double r, double t, DoubleVector& table) const;

   double p_int_r_i(uint i, double r, double t, DoubleVector& table) const;

   void create_p_int_r_Table(double t, uint maxi, DoubleVector& table) const;

   double get_p_int_r_Table_i(uint& i, double t, DoubleVector& table) const
   {
      if (i >= table.size())
         create_p_int_r_Table(t, i + 1, table);
      return table[i];
   }

   const double v_;         // The diffusion constant and drift velocity
   const double sigma_;     // These are the dimensions of our domain; L is calculated as a-sigma
   const double a_;
   const double r0_;
   const double lscale_;       // This is the 'length scale' of your system (1e-14 or 1e6), Although rescaling is discontinued, we use it to check whether a is well-chosen
   const double tscale_;       // This is the time scale of the system, used by drawTime_f
};

// --------------------------------------------------------------------------------------------------------------------------------
