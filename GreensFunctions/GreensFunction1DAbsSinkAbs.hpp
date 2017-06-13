#ifndef GREENSFUNCTION1DABSSINKABS_HPP
#define GREENSFUNCTION1DABSSINKABS_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction1DAbsSinkAbs : public GreensFunction
{
public:
   GreensFunction1DAbsSinkAbs(double D, double k, double r0, double rsink, double sigma, double a)
      : GreensFunction(D), k_(k), r0_(r0), sigma_(sigma), a_(a), rsink_(rsink), tscale_(GfCfg.T_TYPICAL_1D),
      Lr_(r0 >= rsink ? a - rsink : rsink - sigma), Ll_(r0 >= rsink ? rsink - sigma : a - rsink), L0_(fabs(r0 - rsink))
   {
      /* Set variables which define a domain with the sink at the origin.
         Furthermore r0 is assumed to be right from the sink. */
      THROW_UNLESS(std::invalid_argument, a > sigma);
      calculate_n_roots(1);
   }

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction1DAbsSinkAbs"; }

   double geta() const { return a_; }

   double getLr() const { return Lr_; }

   double getLl() const { return Ll_; }

   double getL0() const { return L0_; }

   double getrsink() const { return rsink_; }

   double getsigma() const { return sigma_; }

   double getr0() const { return r0_; }

   double getk() const { return k_; }

   /* Determine which event has occurred at time t. Either an escape
      (left or right boundary) or a reaction with the sink.
      Based on the fluxes through the boundaries at the given time. */
   EventKind drawEventType(double rnd, double t) const override;

   /* Draws the first passage time from the propensity function */
   double drawTime(double rnd) const override;

   /* Draws the position of the particle at a given time, assuming that
      the particle is still in the domain. */
   double drawR(double rnd, double t) const override;

private:


   /* Calculates the probability density of finding the particle at
   location z at timepoint t, given that the particle is still in the
   domain. */
   double calcpcum(double r, double t) const;

   /* Calculates the probability flux leaving the domain through the right
      absorbing boundary at time t. */
   double flux_leavea(double t) const;

   /* Calculates the probability flux leaving the domain through the left
      absorbing boundary at time t. */
   double flux_leaves(double t) const;

   /* Calculates the probability flux leaving the domain through the sink
      at time t. */
   double flux_sink(double t) const;

   /* Calculates the probability of finding the particle inside the
      domain at time t -> the survival probability */
   double p_survival(double t) const;

   /* c.d.f. of the greensfunction with respect to position. */
   double p_int_r_table(double r, double t, DoubleVector& table) const;

   double p_int_r(double r, double t) const;

   /* Calculates the total probability flux leaving the domain at time t. */
   double flux_tot(double t) const;

   /* Calculates the probability flux leaving the domain through the
      sub-domain containing r0 via the absorbing boundary and the flux
      leaving the sub-domain not containing r0 via the absorbing boundary. */
   double flux_abs_Lr(double t) const;
   double flux_abs_Ll(double t) const;

   /* Calculates the probability density of finding the particle at
      location r at time t. */
   double prob_r(double r, double t) const;


   //struct root_f_params
   //{
   //   double Lm_L;
   //   double h;
   //};

   struct lower_upper_params
   {
      double h;
      double Lm_L;
      double long_period;
      double short_period;
      double last_long_root;
      double last_short_root;
      bool last_was_long;
   };

   /* Functions managing the rootList */

   /* return the rootList size */
   uint rootList_size() const { return static_cast<uint>(rootList_.size()); };

   /* return the n + 1'th root */
   double get_root(uint n) const;

   /* Check the rootList for the first n roots. */
   void calculate_n_roots(uint n) const;

   /* Fills the rootList from i to n. */
   void fill_table_to_n(uint i, uint n);

   /* Function returns two positions on the x-axis which straddle the next root. */
   DoublePair get_lower_and_upper();

   /* Function of which we need the roots. */
   double root_f(double x, double Lm_L, double h) const;

   /* Guess the number of terms needed for convergence, given t. */
   uint guess_maxi(double t) const;

   /* Function for calculating the survival probability. */
   double p_survival_table(double t, DoubleVector& psurvTable) const;

   double p_survival_i(uint i, double t, const DoubleVector& table) const;

   double p_survival_table_i(double root_i) const;

   void createPsurvTable(DoubleVector& table) const;

   /* Functions for calculating the greensfunction. */
   double prob_r_r0_i(uint i, double rr, double t) const;

   double prob_r_nor0_i(uint i, double rr, double t) const;

   /* Functions for calculating the fluxes. */
   double flux_tot_i(uint i, double t) const;

   double flux_abs_Lr_i(uint i, double t) const;

   double flux_abs_Ll_i(uint i, double t) const;

   /* functions for calculating the c.d.f. */

   /* i'th term of p_int_r(r') for r' in left domain */
   double p_int_r_leftdomain(uint i, double rr, double t, DoubleVector& table) const;

   /* i'th term of p_int_r(r') for r' in right domain, left of r0 */
   double p_int_r_rightdomainA(uint i, double rr, double t, DoubleVector& table) const;

   /* i'th term of p_int_r(r') for r' in right domain, right of r0 */
   double p_int_r_rightdomainB(uint i, double rr, double t, DoubleVector& table) const;

   /* Fills table with r-independent part of p_int_r_i. */
   void create_p_int_r_Table(double t, DoubleVector& table) const;

   /* Returns i'th r-independent term of p_int_r_i.
      Term is created if not in table. */
   double get_p_int_r_Table_i(uint& i, double t, DoubleVector& table) const
   {
      if (i >= table.size())
      {
         calculate_n_roots(i + 1);
         create_p_int_r_Table(t, table);
      }
      return table[i];
   }

   /* Denominator of the Greens function */
   inline double p_denominator_i(double root_n) const;

   /* Standard form of Greens Function: exp( -Dt root_n ** 2 ) / denominator */
   inline double p_exp_den_i(double t, double root_n, double root_n2) const;

   /* Function for drawR */
   static double drawT_f(double t, void *p);

   /* Function for drawTime */
   static double drawR_f(double r, void *p);

private:
   const double k_;           // The reaction constant
   const double r0_;          // starting position
   const double sigma_;       // The left and right boundary of the domain (sets the l_scale, see below)
   const double a_;
   const double rsink_;       // Position of the sink in the domain.
   const double tscale_;           // This is the time scale of the system.

   /* Greensfunction assumes that the sink is at the origin, and
      consists of two sub-domains: one between a boundary and the sink including
      r0, and one between boundary and sink not including r0. */

   const double Lr_;            // Length of sub-domain which does not include r0.
   const double Ll_;            // Length of sub-domain which does include r0.
   const double L0_;            // Distance between the sink and r0.

   DoubleVector rootList_;                        // Stores all the roots.
   lower_upper_params luparams;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif // GREENSFUNCTION1DRADABS_HPP