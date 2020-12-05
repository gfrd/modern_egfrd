#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction1DRadAbs : public GreensFunction
{
public:
   GreensFunction1DRadAbs(double D, double k, double r0, double sigma, double a, double v = 0.0)
      : GreensFunction(D), v_(v), k_(k), r0_(r0), sigma_(sigma), a_(a), lscale_(GfCfg.L_TYPICAL_1D), tscale_(GfCfg.T_TYPICAL_1D)
   {
      //set first root.
      calculate_n_roots(1);
   }

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction1DRadAbs"; }

   double geta() const { return a_; }

   double getsigma() const { return sigma_; }

   double getr0() const { return r0_; }

   double getk() const { return k_; }

   double getv() const { return v_; }

   // Calculates the probability density of finding the particle at 
   // location z at timepoint t, given that the particle is still in the 
   // domain.
   double calcpcum(double r, double t) const;

   // Determine which event has occurred, an escape or a reaction. Based 
   // on the fluxes through the boundaries at the given time. Beware: if 
   // t is not a first passage time you still get an answer!
   EventKind drawEventType(double rnd, double t) const override;

   // Draws the first passage time from the propensity function
   double drawTime(double rnd) const override;

   // Draws the position of the particle at a given time, assuming that 
   // the particle is still in the domain
   double drawR(double rnd, double t) const override;

   // These methods are both public and private, they are used by public methods 
   // but can also be called from the 'outside'. This is mainly because of 
   // debugging purposes.

   // Calculates the probability of finding the particle inside the 
   // domain at time t -> the survival probability
   double p_survival(double t) const;

   // Calculates the total probability flux leaving the domain at time t
   double flux_tot(double t) const;

   // Calculates the probability flux leaving the domain through the 
   // radiative boundary at time t
   double flux_rad(double t) const;

   // Calculates the flux leaving the domain through the radiative 
   // boundary as a fraction of the total flux. This is the probability 
   // that the particle left the domain through the radiative
   // boundary instead of the absorbing boundary.
   double fluxRatioRadTot(double t) const;

   // Calculates the probability density of finding the particle at 
   // location r at time t.
   double prob_r(double r, double t) const;

private:

   double An(double a_n) const;
   double Bn(double a_n) const;
   double Cn(double a_n, double t) const;

   struct tan_f_params
   {
      double a;
      double h;
   };

   struct drawT_params
   {
      GreensFunction1DRadAbs const* gf;
      DoubleVector& psurvTable;
      double rnd;
   };

   struct drawR_params
   {
      GreensFunction1DRadAbs const* gf;
      const double t;
      DoubleVector table;
      double rnd;
   };

   /* Functions managing the rootList */

   /* return the rootList size */
   uint rootList_size() const { return static_cast<uint>(rootList_.size()); }

   /* return the n + 1'th root */
   double get_root(uint n) const;

   /* Check the rootList for the first n roots. */
   void calculate_n_roots(uint n) const;

   /* Fills the rootList from i to n. */
   void fill_table_to_n(uint i, uint n);

   /* Guess the number of terms needed for convergence, given t. */
   uint guess_maxi(double t) const;

   /* this is the appropriate definition of the function in gsl. */
   static double tan_f(double x, void *p);

   /* functions for drawTime / p_survival */

   static double drawT_f(double t, void *p);

   double p_survival_table(double  t, DoubleVector& psurvTable) const;

   double p_survival_i(uint i, double t, const DoubleVector& table) const;

   double p_survival_table_i_v(uint i) const;

   double p_survival_table_i_nov(uint i) const;

   void createPsurvTable(DoubleVector& table) const;

   /* functions for drawR */

   static double drawR_f(double z, void* p);

   double p_int_r_table(double r, double t, DoubleVector& table) const;

   double p_int_r_i(uint i, double r, double t, DoubleVector& table) const;

   void create_p_int_r_Table(double t, DoubleVector& table) const;

   double get_p_int_r_Table_i(uint& i, double t, DoubleVector& table) const
   {
      if (i >= table.size())
      {
         calculate_n_roots(i + 1);
         create_p_int_r_Table(t, table);
      }

      return table[i];
   }

   const double v_;         // The diffusion constant and drift velocity
   const double k_;         // The reaction constant
   const double r0_;        // Initial distance of particle
   const double sigma_;     // The left and right boundary of the domain (sets the l_scale, see below)
   const double a_;
   const double lscale_;   // This is the length scale of the system
   const double tscale_;   // This is the time scale of the system.

   DoubleVector rootList_;            /* vector containing the roots 0f tan_f. */
};

// --------------------------------------------------------------------------------------------------------------------------------
