#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <tuple>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_lambert.h>
#include "helperFunctionsGf.hpp"
#include "GreensFunction2DRadAbs.hpp"
#include "Logger.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

const uint GreensFunction2DRadAbs::MAX_ALPHA_SEQ = 500;				  // The maximum number of n 
const double GreensFunction2DRadAbs::SCAN_START = 0.001;
const double GreensFunction2DRadAbs::FRACTION_SCAN_INTERVAL = .5; // TODO CHANGED THIS FROM .5 to .2
const double GreensFunction2DRadAbs::CONVERGENCE_ASSUMED = 25;
const double GreensFunction2DRadAbs::INTERVAL_MARGIN = .33;

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction2DRadAbs");

// --------------------------------------------------------------------------------------------------------------------------------

// Greens function class for 2d Green's Function for 2d annulus with radial and
// axial dependence. Inner boundary is radiative (rad) (reaction event), outer 
// boundary is absorbing (abs) (escape event). Different "draw" functions 
// provide a way to draw certain values from the Green's Function, e.g. an
// escape angle theta ("drawTheta" function).
// 
// Based upon code from Riken Institute. 
// Written by Laurens Bossen, Adapted by Martijn Wehrens. FOM Institute AMOLF.

GreensFunction2DRadAbs::GreensFunction2DRadAbs(double D, double kf, double r0, double sigma, double a) :
   PairGreensFunction(D, kf, r0, sigma), h_(kf / (D * 2.0 * M_PI * sigma)), a_(a), estimated_alpha_root_distance_(M_PI / (a - sigma))
   // observed convergence of
   // distance roots f_alpha().
   // ^: Here parent "constructors" are specified that are executed 
   // also a constructor initialization list can be (and is) specified. 
   // These variables will be set before the contents of the constructor are 
   // called.
{
   // Check whether input makes sense, outer boundary a should be > inner 
   // boundary sigma.
   THROW_UNLESS(std::invalid_argument, a > sigma);

   // Clear AlphaTables
   clearAlphaTable();
}

//
// Alpha-related methods
//

// Resets the alpha-tables
void GreensFunction2DRadAbs::clearAlphaTable() const
{
   // Clears all vectors in the alphaTables_
   std::for_each(alphaTables_.begin(), alphaTables_.end(), [](DoubleVector rv) { rv.clear(); });

   // Sets all values of the alpha_x_scan_table_ to zero.
   std::fill(alpha_x_scan_table_.begin(), alpha_x_scan_table_.end(), SCAN_START*estimated_alpha_root_distance_);
   // DEBUG TODO; actual number here might
   // be important for functioning of Bessel functions.

   // Sets all values of the alpha_correctly_estimated_ table to zero.
   std::fill(alpha_correctly_estimated_.begin(), alpha_correctly_estimated_.end(), 0);
}

// The method evaluates the equation for finding the alphas for given alpha. This
// is needed to find the alpha's at which the expression is zero -> alpha is the root.
double GreensFunction2DRadAbs::f_alpha0(double alpha) const
{
   double s_An(sigma_ * alpha);
   double a_An(a_ * alpha);
   // Needed? TODO
   //  double h_s( h * sigma);

   double J0_s_An(gsl_sf_bessel_J0(s_An));
   double J1_s_An(gsl_sf_bessel_J1(s_An));
   double J0_a_An(gsl_sf_bessel_J0(a_An));

   double Y0_s_An(gsl_sf_bessel_Y0(s_An));
   double Y1_s_An(gsl_sf_bessel_Y1(s_An));
   double Y0_a_An(gsl_sf_bessel_Y0(a_An));

   //  double rho1 ( ( (h_s * J0_s_An) + (s_An * J1_s_An) ) * Y0_a_An );
   //  double rho2 ( ( (h_s * Y0_s_An) + (s_An * Y1_s_An) ) * J0_a_An );
   //  return (rho1 - rho2); 

   // Sigma can be divided out, roots will remain same:
   // (Note: actually double checked this).
   double rho1(((h_ * J0_s_An) + (alpha * J1_s_An)) * Y0_a_An);
   double rho2(((h_ * Y0_s_An) + (alpha * Y1_s_An)) * J0_a_An);

   return rho1 - rho2;
}

// f_alpha() Calculates the value of the mathematical function f_alpha(). The
// roots (y=0) of this function are constants in the Green's Functions.
double GreensFunction2DRadAbs::f_alpha(double alpha, int n) const
{
   double s_An(sigma_ * alpha);
   double a_An(a_ * alpha);
   double h_sigma(h_ * sigma_);

   double Jn_s_An(gsl_sf_bessel_Jn(n, s_An));
   double Jn1_s_An(gsl_sf_bessel_Jn(n + 1, s_An));
   double Jn_a_An(gsl_sf_bessel_Jn(n, a_An));

   double Yn_s_An(gsl_sf_bessel_Yn(n, s_An));
   double Yn1_s_An(gsl_sf_bessel_Yn(n + 1, s_An));
   double Yn_a_An(gsl_sf_bessel_Yn(n, a_An));

   double rho1(((h_sigma * Jn_s_An) + (s_An * Jn1_s_An) - n*Jn_s_An) * Yn_a_An);
   double rho2(((h_sigma * Yn_s_An) + (s_An * Yn1_s_An) - n*Yn_s_An) * Jn_a_An);
   return (rho1 - rho2);

   //  Or..?
   //  double rho1 ( ( ((h*sigma-n) * Jn_s_An) + (s_An * Jn1_s_An) ) * Yn_a_An );
   //  double rho2 ( ( ((h*sigma-n) * Yn_s_An) + (s_An * Yn1_s_An) ) * Jn_a_An );
   //  return (rho1 - rho2);

}

// calculates the constant part of the i-th term for the survival probability
double GreensFunction2DRadAbs::p_survival_i(double alpha) const
{
   // get the needed parameters
   double s_An(sigma_*alpha);
   double a_An(a_*alpha);
   double alpha_sq(alpha*alpha);

   // calculate all the required Bessel functions
   double J1_sAn(gsl_sf_bessel_J1(s_An));
   double J0_aAn(gsl_sf_bessel_J0(a_An));
   double Y0_aAn(gsl_sf_bessel_Y0(a_An));
   double Y1_sAn(gsl_sf_bessel_Y1(s_An));

   // calculate C0,n
   double C_i_0(calc_A_i_0(alpha));

   // calculate the integral over Bn,0
   double dB_n_0dr(J1_sAn*Y0_aAn - Y1_sAn*J0_aAn);
   // this is only the part without alpha of dB0,n(sigma)/dr
   double B_n_0_int(double(2.0) / (M_PI*alpha_sq) - (sigma_ / alpha)*dB_n_0dr);

   // return the total result
   double result(C_i_0 * B_n_0_int);

   return result;
}

// Calculates the factor An,0 for (for example) determination of the flux through the outer interface
double GreensFunction2DRadAbs::calc_A_i_0(double alpha) const
{
   double s_An(sigma_*alpha);
   double a_An(a_*alpha);
   double r0An(r0_*alpha);

   // Calculate all the required Bessel functions
   double J0_sAn(gsl_sf_bessel_J0(s_An));
   double J1_sAn(gsl_sf_bessel_J1(s_An));
   double J0_aAn(gsl_sf_bessel_J0(a_An));

   double J0_r0An(gsl_sf_bessel_J0(r0An));
   double Y0_aAn(gsl_sf_bessel_Y0(a_An));
   double Y0_r0An(gsl_sf_bessel_Y0(r0An));

   // Calculate return value
   double alpha_sq(alpha*alpha);

   double rho(h_*J0_sAn + alpha*J1_sAn);
   double rho_sq(rho*rho);

   double B_n_0(J0_r0An*Y0_aAn - Y0_r0An*J0_aAn);

   double A_i_0((alpha_sq * rho_sq * B_n_0) / (rho_sq - (J0_aAn*J0_aAn)*(h_*h_ + alpha_sq)));

   return A_i_0;
}

// Calculates the n-th term of the summation for calculating the flux through  the inner interface (reaction)
double GreensFunction2DRadAbs::leaves_i(double alpha) const
{
   double s_An(sigma_*alpha);
   double a_An(a_*alpha);

   // calculate all the required Bessel functions
   double J1_sAn(gsl_sf_bessel_J1(s_An));
   double Y1_sAn(gsl_sf_bessel_Y1(s_An));
   double J0_aAn(gsl_sf_bessel_J0(a_An));
   double Y0_aAn(gsl_sf_bessel_Y0(a_An));

   // calculate An,0
   double A_i_0(calc_A_i_0(alpha)); // calculate the coefficient A0,n

   // calculate dBn,0(sigma)/dr
   double dB_n_0dr(-alpha*(J1_sAn*Y0_aAn - Y1_sAn*J0_aAn));

   // calculate the total result
   double result(A_i_0 * dB_n_0dr);

   return result;
}

// calculates a table with all the constant factors for the survival probability
void GreensFunction2DRadAbs::createPsurvTable(DoubleVector& table) const
{
   const DoubleVector& alphaTable_0(getAlphaTable(0));	// get the roots for the survival probability

   table.clear();                              // empty the table
   table.reserve(alphaTable_0.size());       // and get the necessary memory

   std::transform(alphaTable_0.begin(), alphaTable_0.end(), std::back_inserter(table), std::bind(&GreensFunction2DRadAbs::p_survival_i, this, std::placeholders::_1));
   // This gets all the roots from 'begin' to 'end' passes them as an argument to p_survival_i and the result is passed to back_inserter
}

// Creates the tables with various Bessel functions used in drawR, the table is used to speed things up
void GreensFunction2DRadAbs::createY0J0Tables(DoubleVector& Y0_Table, DoubleVector& J0_Table, DoubleVector& Y0J1J0Y1_Table, double t) const
{
   const DoubleVector& alphaTable_0(getAlphaTable(0));
   // get the roots for the survival probability
   Y0_Table.clear();                           // empty the table
   J0_Table.clear();
   Y0J1J0Y1_Table.clear();

   Y0_Table.reserve(alphaTable_0.size());    // and get the necessary memory
   J0_Table.reserve(alphaTable_0.size());
   Y0J1J0Y1_Table.reserve(alphaTable_0.size());

   std::tuple<double, double, double> result;

   for (uint count = 0; count < alphaTable_0.size(); count++)
   {
      result = Y0J0J1_constants(alphaTable_0[count], t);

      Y0_Table.push_back(std::get<0>(result));
      J0_Table.push_back(std::get<1>(result));
      Y0J1J0Y1_Table.push_back(std::get<2>(result));
   }
}

// Creates the values for in the tables Y0, J0 and Y0J1J0Y1
std::tuple<double, double, double> GreensFunction2DRadAbs::Y0J0J1_constants(double alpha, double t) const
{
   double s_An(sigma_*alpha);
   double a_An(a_*alpha);
   double alpha_sq(alpha*alpha);

   // calculate all the required Bessel functions	
   double J0_aAn(gsl_sf_bessel_J0(a_An));
   double Y0_aAn(gsl_sf_bessel_Y0(a_An));
   double J1_sAn(gsl_sf_bessel_J1(s_An));
   double Y1_sAn(gsl_sf_bessel_Y1(s_An));

   // calculate An,0
   double A_i_0(calc_A_i_0(alpha)); //_sq * rho_sq * B_n_0)/( rho_sq - J0_bAn*J0_bAn*(h*h + alpha_sq)));

   // calculate the exponent with the time
   double expT(std::exp(-D_*alpha_sq*t));
   // and the product
   double Ai0_expT(A_i_0 * expT / alpha);

   // calculate the large constant term in the intergral of Bn,0
   double Y0J1_J0Y1(Y0_aAn*sigma_*J1_sAn - J0_aAn*sigma_*Y1_sAn);

   return std::make_tuple(Ai0_expT*Y0_aAn, Ai0_expT*J0_aAn, Ai0_expT*Y0J1_J0Y1);
}

// =============================================================================
// Root finding algorithm for alpha function. 
// Roots (y=0) are necessary to calculate 
// Modified by wehrens@amolf.nl. nov, 2011
//
// Note mar 2012:
// The modified root finding algorithm can only work if the roots are calculated
// in a "chronological" sequence (i.e. i = 0, i = 1, i = 2 etc). Therefore roots
// are now all calculated immediately when the roots-table is expanded.
// =============================================================================

// Scans for next interval where a sign change is observed, and thus where a root is expected. 
void GreensFunction2DRadAbs::GiveRootInterval(double& low, double& high, int n) const     // Order of Bessel functions
{
   THROW_UNLESS(std::invalid_argument, static_cast<uint>(n) < alpha_x_scan_table_.size());

   // Variables for function values @ resp. left and right boundary interval.

   // # Get/calculate boundaries and determine width of interval to check for sign change.
   double interval = FRACTION_SCAN_INTERVAL * estimated_alpha_root_distance_;

   // If order (n) is zero, the offset is zero, otherwise take n-1 offset value as first estimate.
   // (This can be done because the offsets only get bigger, the roots shift to the right with the order of the Bessel functions n.)
   if (alpha_x_scan_table_[n] == 0)  // which implies i == 0
   {
      if (n > 0)
      {
         alpha_x_scan_table_[n] = (alphaTables_[n - 1][0]);
      }
   }

   // Define new interval as between x1=offset ("low") and x2=(offset+interval)
   // ("high"). Make sure "low" is > 0.
   low = alpha_x_scan_table_[n];
   high = alpha_x_scan_table_[n] + interval;
   if (low <= 0) // TODO this check should be redundant
   {
      _log.error() << "Left alpha search interval boundary < 0.";//low = EPSILON/L_TYPICAL;
      throw std::exception();
   }

   // # Look for the sign change:

   // Get the values of the function at x "low" and x "high".
   // Different for n=0 because this results in a simpler function.
   // (Note: This code could be optimized by duplicating all involved functions and removing this if-statement. 

   double f_low = n == 0 ? f_alpha0(low) : f_alpha(low, n);
   double f_high = n == 0 ? f_alpha0(high) : f_alpha(high, n);

   // Continue shifting the search interval until a sign change is detected.
   while (f_low * f_high > 0)
   {
      low = high;
      f_low = f_high;

      high += interval;
      f_high = f_alpha(high, n);
   }

   // When above loop has finished, low and high have values in-between which
   // a root of the alpha function should lie have been found. Make sure that 
   // scanning for the next root starts at the end of the domain [low, high] found here.
   alpha_x_scan_table_[n] = high;
}

// Simply returns an interval based upon previous root, estimated interval in-between roots and INTERVAL_MARGIN (see .hpp).
void GreensFunction2DRadAbs::GiveRootIntervalSimple(double& low, double& high, int n, uint i) const
{
   // Offset is simply based on previous root, the interval in which the first 
   // root (i=0) lies is never calculated with this function.
   double previous_root = getAlpha(n, i - 1);

   // Calculates interval [low, high] where root is expected based on the 
   // assumption where in a converging regime, where the deviation from this
   // estimate is not more than INTERVAL_MARGIN.
   low = previous_root + estimated_alpha_root_distance_* (1 - INTERVAL_MARGIN);
   high = previous_root + estimated_alpha_root_distance_ * (1 + INTERVAL_MARGIN);
}

// This function calls the GSL root finder, for roots for which n = 0. (This is a special case for which the function simplifies.)
double GreensFunction2DRadAbs::getAlphaRoot0(double low, double high) const
{
   //f_alpha0_aux_params params = { *this };
   //gsl_function F = { reinterpret_cast<double(*)(double, void*)> (&GreensFunction2DRadAbs::f_alpha0_aux_F), &params };
   auto f = [this](double alpha) { return f_alpha0(alpha); };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(GfCfg.EPSILON / GfCfg.L_TYPICAL_2D, GfCfg.EPSILON, "GreensFunction2DRadAbs::getAlphaRoot0");
}

// This function calls the GSL root finder, for roots for which n > 0. (n = 0 is  a special case for which the function simplifies.)
double GreensFunction2DRadAbs::getAlphaRootN(double low, double high, int n) const
{
   // n is the summation index (the order of the Bessel functions used
   auto f = [this, n](double alpha) { return f_alpha(alpha, n); };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(GfCfg.EPSILON / GfCfg.L_TYPICAL_2D, GfCfg.EPSILON, "GreensFunction2DRadAbs::getAlphaRootN");
}

// Simply calls the correct function to call the rootfinder (either getAlphaRoot0 or getAlphaRootN).
double GreensFunction2DRadAbs::getAlphaRoot(double low, double high, int n) const
{
   return n == 0 ? getAlphaRoot0(low, high) : getAlphaRootN(low, high, n);
}

// Subfunction of getAlpha
// 
// The function just counts the number of times the root lies in the interval 
// that can be used to estimate the next root. If this happens enough, then it
// edits alpha_correctly_estimated_ to a negative value to signify that it can 
// be assumed that for all next roots the distance inbetween roots will
// equal the interval +/- margin.
// TODO: This could be more sophisticated.
void GreensFunction2DRadAbs::decideOnMethod2(uint n, uint i) const
{
   // Since the function can only decide with two alpha's already calculated, 
   // it can't do anything if i = 0.
   if (i > 0)
   {
      double dx(getAlpha(n, i) - getAlpha(n, i - 1)); // note the recursiveness!

      // If the relative deviation from the expected difference is smaller 
      // than the expected margin, increase number of would-be correct 
      // guesses.
      if (fabs(1 - dx / estimated_alpha_root_distance_) < INTERVAL_MARGIN)
      {
         ++alpha_correctly_estimated_[n];
      }
      else
      {
         alpha_correctly_estimated_[n] = 0;
      }

      // If guessing would have worked for last CONVERGENCE_ASSUMED roots, 
      // assume it will for all following roots.
      if (alpha_correctly_estimated_[n] > CONVERGENCE_ASSUMED) {
         alpha_x_scan_table_[n] = -2; // permanently switch
      }
   }
}

// This function searches for roots (y=0) of the so-called alpha-function 
// (::f_alpha()). It either moves a small search interval along the x-axis to
// check for sign-change (which would indicate a root), and calls the GSL
// root finder, or directly calls the root finder if the spacing between 
// roots is found to be converged. 
double GreensFunction2DRadAbs::getAlpha(uint n, uint i) const
{
   THROW_UNLESS(std::invalid_argument, n < alphaTables_.size());

   double current_root_, low, high;

   // # "Administration"
   // Equals reading/writing to/from alphaTables_ to reading/writing from/to to 
   // alphaTables_[n]; n being the order of the Bessel function.
   DoubleVector& alphaTable(this->alphaTables_[n]);

   // Gets it's size 
   uint oldSize(static_cast<uint>(alphaTable.size()));

   // # Expansion of root table

   // If doesn't contain requested value, expand table until value
   if (i >= oldSize)
   {
      // Expand the table, temporarily fill with zeroes
      alphaTable.resize(i + 1, 0);

      // # Calculating the necessary roots to expand the table
      for (uint j = oldSize; j <= i; j++)
      {
         if (alphaTable[j] != 0)
         {
            _log.error() << "tried accessing root that's not 0. Didn't search. i=" << i << ", oldSize=" << oldSize << ", j=" << j;
         }
         else
         {
            // Method 1. SCANNING. If the roots are not expected to lie close enough 
            // to the estimate, use the "scanning" procedure to find an interval 
            // that contains a root. (More robust method.)
            //      If it is established that we can use method 2, 
            // alpha_x_scan_table_[n] will contain a value < 0.
            if (alpha_x_scan_table_[n] >= 0)
            {
               // ### Gets estimate of interval by sign-change-searching
               //      high and low are the return-values of GiveRootInterval.
               GiveRootInterval(low, high, n);

               // ### Finds the root using the GSL rootfinder
               current_root_ = getAlphaRoot(low, high, n);

               // ### Puts the found root in the table.
               alphaTable[j] = current_root_;

               // Check if we can use method 2 for next roots
               decideOnMethod2(n, j);
            }
            // Method 2. ASSUMING ROOTS AT ~FIXED INTERVAL. If next root is expected 
            // to lie at distance close enough to estimated distance the root finder 
            // can be called using a simple estimated interval.
            else
            {
               // ### Get interval by simple extrapolation
               GiveRootIntervalSimple(low, high, n, j);

               // ### Finds the root using the GSL rootfinder
               current_root_ = getAlphaRoot(low, high, n);

               // ### Puts the found root in the table.
               alphaTable[j] = current_root_;
            }
         }
      }

   }
   return alphaTable[i];
}

// calculates the ith term with exponent and time for the survival probability
double GreensFunction2DRadAbs::p_survival_i_exp_table(uint i, double t, const DoubleVector& table) const
{
   double alpha(getAlpha(0, i));
   return std::exp(-D_ * t * alpha * alpha) * table[i];
}

// adds the exponential with the time to the sum. Needed for the calculation of the flux through the outer interface
double GreensFunction2DRadAbs::leavea_i_exp(uint i, double t) const
{
   double alpha(getAlpha(0, i));
   return std::exp(-D_ * t * alpha * alpha) * calc_A_i_0(alpha);
}

// adds the exponential with the time to the sum. Needed for the inner interface (reaction)
double GreensFunction2DRadAbs::leaves_i_exp(uint i, double t) const
{
   double alpha(getAlpha(0, i));
   return std::exp(-D_ * t * alpha * alpha) * leaves_i(alpha);
}

// calculates the Bossen function for a given r
double GreensFunction2DRadAbs::p_int_r_i_exp_table(uint i, double r, const DoubleVector& Y0_aAnTable, const DoubleVector& J0_aAnTable, const DoubleVector& Y0J1J0Y1Table) const
{
   double alpha(getAlpha(0, i));	// get the root An
   double r_An(r*alpha);
   double J1_rAn(gsl_sf_bessel_J1(r_An));
   double Y1_rAn(gsl_sf_bessel_Y1(r_An));
   double result(Y0_aAnTable[i] * r*J1_rAn - J0_aAnTable[i] * r*Y1_rAn - Y0J1J0Y1Table[i]);
   return result;
}

// This tries to guess the maximum number of n iterations it needs for calculating the survival probability
uint GreensFunction2DRadAbs::guess_maxi(double t) const
{
   uint safety(2);
   if (!std::isfinite(t)) return safety;

   double alpha0(getAlpha(0, 0));
   double Dt(D_ * t);
   double thr((exp(-Dt * alpha0 * alpha0) / alpha0) * GfCfg.EPSILON * 1e-1);
   double thrsq(thr * thr);

   if (thrsq <= 0.0) return MAX_ALPHA_SEQ;

   double max_alpha(1.0 / (sqrt(exp(gsl_sf_lambert_W0(2 * Dt / thrsq)) * thrsq)));
   uint maxi(safety + static_cast<uint>(max_alpha * (a_ - sigma_) / M_PI));
   return std::min(maxi, MAX_ALPHA_SEQ);
}

// Calculates the survival probability at a given time.
// This is a little wrapper for the p_survival_table so that you can easily calculate the survival probability at a given time
double GreensFunction2DRadAbs::p_survival(double t) const
{
   DoubleVector psurvTable;
   return p_survival_table(t, psurvTable);
}

// This actually calculates the Survival probability at time t given the particle was at r0 at time 0
// It uses the pSurvTable for efficiency (so you don't have to calculate all the constant factors all
// the time)
double GreensFunction2DRadAbs::p_survival_table(double t, DoubleVector& psurvTable) const
{
   double p;
   uint maxi(guess_maxi(t)); // guess the maximum number of iterations required

   // If the estimated # terms needed for convergence is bigger than number
   // of terms summed over (MAX_ALPHA_SEQ), give error.
   if (maxi == MAX_ALPHA_SEQ)
      _log.warn() << "p_survival_table (used by drawTime) couldn't converge; max terms reached: " << maxi;

   if (psurvTable.size() < maxi + 1)           // if the dimensions are good then this means
   {                                            // that the table is filled
      getAlpha(0, maxi);                      // this updates the table of roots
      createPsurvTable(psurvTable);      // then the table is filled with data
   }
   // A sum over terms is performed, where convergence is assumed. It is not clear if this is a just assumption.
   p = funcSum_all(std::bind(&GreensFunction2DRadAbs::p_survival_i_exp_table, this, std::placeholders::_1, t, psurvTable), maxi);
   return p*M_PI*M_PI_2;
}

// calculates the flux leaving through the inner interface at a given moment
// FIXME: This is inaccurate for small t's!!
double GreensFunction2DRadAbs::leaves(double t) const
{
   double p = funcSum(std::bind(&GreensFunction2DRadAbs::leaves_i_exp, this, std::placeholders::_1, t), MAX_ALPHA_SEQ);
   return M_PI_2*M_PI*D_*sigma_*p; // The minus is not there because the flux is in the negative r
   // direction, and the direction is already accounted for in the derivative of B0,n(r)
   // See also leaves_i
}

// calculates the flux leaving through the outer interface at a given moment
double GreensFunction2DRadAbs::leavea(double t) const
{
   double p = funcSum(std::bind(&GreensFunction2DRadAbs::leavea_i_exp, this, std::placeholders::_1, t), MAX_ALPHA_SEQ);
   return M_PI*D_*p;
}

// calculates the sum of the sequence for drawR based upon the values in the tables and r
double GreensFunction2DRadAbs::p_int_r_table(double r, const DoubleVector& Y0_aAnTable, const DoubleVector& J0_aAnTable, const DoubleVector& Y0J1J0Y1Table) const
{
   double p = funcSum(std::bind(&GreensFunction2DRadAbs::p_int_r_i_exp_table, this, std::placeholders::_1, r, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table), static_cast<uint>(Y0_aAnTable.size()));
   return p*M_PI*M_PI_2;
}

// Wrapper for p_survival_table for the iterator to find the root for drawTime
double GreensFunction2DRadAbs::p_survival_table_F(double t, const p_survival_table_params* params)
{
   return params->rnd - params->gf.p_survival_table(t, params->table);
}

// a wrapper to make p_int_r_table available to the iterator calculating the root
double GreensFunction2DRadAbs::p_int_r_F(double r, const p_int_r_params* params)
{
   return params->gf.p_int_r_table(r, params->Y0_aAnTable, params->J0_aAnTable, params->Y0J1J0Y1Table) - params->rnd;
}

// Draws a first passage time, this could be an escape (through the outer boundary) or a reaction (through
// the inner boundary)
double GreensFunction2DRadAbs::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, sigma_ <= r0_ && r0_ <= a_);

   double t_guess;

   if (r0_ == a_ || a_ == sigma_) // when the particle is at the border or if the PD has no real size
   {
      return 0.0;
   }

   // get some initial guess for the time, dr=sqrt(2dDt) with d
   // the dimensionality (2 in this case)
   double t_Abs(gsl_pow_2(a_ - r0_) / (4.0 * D_));
   if (kf_ == 0.0) // if there was only one absorbing boundary
   {
      t_guess = t_Abs;
   }
   else
   {
      double t_Rad(D_ / gsl_pow_2(kf_ / (2 * M_PI*a_)) + gsl_pow_2(r0_ - sigma_) / D_);
      t_guess = std::min(t_Abs, t_Rad); // take the shortest time to a boundary
   }

   t_guess *= .1;

   double minT(std::min(sigma_ * sigma_ / D_ * GfCfg.MIN_T_FACTOR, t_guess * 1e-7)); // something with determining the lowest possible t

   DoubleVector psurvTable; // this is still empty as of now->not used
   p_survival_table_params params = { *this, psurvTable, rnd };
   gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&p_survival_table_F), &params };

   // put in a upper and lower limit (the picked time cannot be infinite!)
   double low(t_guess);
   double high(t_guess);

   // adjust high and low to make sure that f( low ) and f( high ) straddle.
   double value(GSL_FN_EVAL(&F, t_guess));

   if (value < 0.0)                   // if the function is below zero at the guess the upper
   {                                   // boundary should be moved (passed the zero point)
      do
      {
         high *= 10;
         value = GSL_FN_EVAL(&F, high);

         if (fabs(high) >= 1e10)	// if high time is way too high forget about it
         {
            _log.warn() << "Couldn't adjust high. F(" << high << ")=" << GSL_FN_EVAL(&F, high) << ", r0=" << r0_;
            throw std::exception();
         }
      } while (value < 0.0);
   }
   else                                // if the function is over zero (or at zero!) then the lower
   {                                   // boundary should be moved
      double value_prev(value);
      do
      {
         low *= .1;      // keep decreasing the lower boundary until the function straddles
         value = GSL_FN_EVAL(&F, low);     // get the accompanying value

         if (fabs(low) <= minT || fabs(value - value_prev) < GfCfg.EPSILON*GfCfg.T_TYPICAL_2D)
         {
            _log.warn() << "Couldn't adjust low. F(" << low << ")=" << value;
            return low;
         }
         value_prev = value;
      } while (value >= 0.0);
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(GfCfg.EPSILON * GfCfg.T_TYPICAL_2D, GfCfg.EPSILON, "GreensFunction2DRadAbs::drawTime");
}

// --------------------------------------------------------------------------------------------------------------------------------

// This determines based on the flux at a certain time, if the 'escape' was a reaction or a proper escape
GreensFunction::EventKind GreensFunction2DRadAbs::drawEventType(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, 0 <= rnd && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, sigma_ <= r0_ && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t > 0.0);

   if (kf_ == 0.0) // if there cannot be any flow through the radiating boundary it is always an escape
   {
      return EventKind::IV_ESCAPE;
   }

   // First, check if r0 is close only either to a or sigma relative
   // to Dt.  In such cases, the event type is always ESCAPE or REACTION
   // respectively.   This avoids numerical instability in calculating
   // leavea() and/or leaves().

   // Here, use a rather large threshold for safety.
   uint H(6); 				// 6 times the msd travelled as threshold
   double max_dist(H * sqrt(4.0 * D_ * t));
   double a_dist(a_ - r0_);
   double s_dist(r0_ - sigma_);

   if (a_dist > max_dist)
   {
      if (s_dist < max_dist)
      {
         return EventKind::IV_REACTION;
      }
   }
   else // a_dist < max_dist
   {
      if (s_dist > max_dist)
      {
         return EventKind::IV_ESCAPE;
      }
   }

   double reaction(leaves(t));	// flux through rad boundary
   double escape(leavea(t));	// flux through abs boundary
   double value(reaction / (reaction + escape));

   return rnd <= value ? EventKind::IV_REACTION : EventKind::IV_ESCAPE;
}

// This draws a radius R at a given time, provided that the particle was at r0 at t=0
double GreensFunction2DRadAbs::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_ && r0_ < a_);

   if (t == 0.0) return r0_;// if no time has passed

   double psurv(p_survival(t)); // calculate the survival probability at this time
   // this is used as the normalization factor
   // BEWARE!!! This also produces the roots An and therefore
   // SETS THE HIGHEST INDEX -> side effect
   // VERY BAD PROGRAMMING PRACTICE!!

   DoubleVector Y0_aAnTable;
   DoubleVector J0_aAnTable;
   DoubleVector Y0J1J0Y1Table;
   createY0J0Tables(Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, t);

   p_int_r_params params = { *this, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, rnd * psurv };
   gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&p_int_r_F), &params };

   // adjust low and high starting from r0.
   // this is necessary to avoid root finding in the long tails where
   // numerics can be unstable.
   double low(r0_);             // start with the initial position as the first guess
   double high(r0_);
   double value;
   uint H(3);

   double msd(sqrt(4.0 * D_ * t));
   if (GSL_FN_EVAL(&F, r0_) < 0.0)
   {
      do
      {
         high = r0_ + H * msd;
         if (high > a_)
         {
            if (GSL_FN_EVAL(&F, a_) < 0.0)        // something is very wrong, this should never happen
            {
               _log.warn() << "drawR: p_int_r_table( a ) < 0.0. Returning a.";
               return a_;
            }
            high = a_;
            break;
         }
         value = GSL_FN_EVAL(&F, high);
         ++H;
      } while (value < 0.0);

   }
   else
   {
      do
      {
         low = r0_ - H * msd;
         if (low < sigma_)
         {
            if (GSL_FN_EVAL(&F, sigma_) > 0.0)
            {
               _log.warn() << "drawR: p_int_r_table( sigma ) > 0.0. Returning sigma.";
               return sigma_;
            }

            low = sigma_;
            break;
         }

         value = GSL_FN_EVAL(&F, low);
         ++H;
      } while (value > 0.0);
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(GfCfg.L_TYPICAL_2D*GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction2DRadAbs::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

// The calculates constant factor m,n for the drawing of theta. These factors are summed later.
double GreensFunction2DRadAbs::p_m_alpha(uint n, uint m, double r, double t) const
{
   double alpha(getAlpha(m, n)); // Gets the n-th root using the 
   // besselfunctions of order m.

   double alpha_sq(alpha * alpha);
   double msq(m * m);
   double ssq(sigma_ * sigma_);

   double s_Anm(sigma_ * alpha);
   double a_Anm(a_ * alpha);
   double r0Anm(r0_ * alpha);
   double r_Anm(r * alpha);

   //  const CylindricalBesselGenerator& s(CylindricalBesselGenerator::instance());

   // calculate the needed bessel functions
   double Jm_sAnm(gsl_sf_bessel_Jn(m, s_Anm));
   double Jmp1_sAnm(gsl_sf_bessel_Jn(m + 1, s_Anm));	// prime
   //  double Jm_sAnm   (s.J(m, s_Anm));
   //  double Jmp1_sAnm (s.J(m+1, s_Anm));	// prime

   double Jm_aAnm(gsl_sf_bessel_Jn(m, a_Anm));
   double Ym_aAnm(gsl_sf_bessel_Yn(m, a_Anm));
   //  double Jm_aAnm   (s.J(m, a_Anm));
   //  double Ym_aAnm   (s.Y(m, a_Anm));

   double Jm_r0Anm(gsl_sf_bessel_Jn(m, r0Anm));
   double Ym_r0Anm(gsl_sf_bessel_Yn(m, r0Anm));
   //  double Jm_r0Anm  (s.J(m, r0Anm));
   //  double Ym_r0Anm  (s.Y(m, r0Anm));

   double Jm_rAnm(gsl_sf_bessel_Jn(m, r_Anm));
   double Ym_rAnm(gsl_sf_bessel_Yn(m, r_Anm));
   //  double Jm_rAnm   (s.J(m, r_Anm));
   //  double Ym_rAnm   (s.Y(m, r_Anm));

   // calculating An,m
   double h_ma(h_ - m / sigma_);
   double rho(h_ma*Jm_sAnm + alpha*Jmp1_sAnm);
   double rho_sq(rho*rho);
   // calculating Bn,m(r')
   double B_n_m_r0(Jm_r0Anm * Ym_aAnm - Ym_r0Anm * Jm_aAnm);

   double A_n_m(alpha_sq * rho_sq * B_n_m_r0 / (rho_sq - Jm_aAnm*Jm_aAnm*(h_*h_ + alpha_sq - msq / ssq)));

   // calculating Bn,m(r*)
   double B_n_m_r(Jm_rAnm * Ym_aAnm - Ym_rAnm * Jm_aAnm);

   // calculating the result
   double result(A_n_m * B_n_m_r * exp(-D_*alpha_sq*t));

   return result;
}

// This calculates the m-th constant factor for the drawTheta method. 
double GreensFunction2DRadAbs::p_m(int m, double r, double t) const
{
   // The m-th factor is a summation over n
   return funcSum(std::bind(&GreensFunction2DRadAbs::p_m_alpha, this, std::placeholders::_1, m, r, t), MAX_ALPHA_SEQ, GfCfg.EPSILON);
}

// this should make the table of constants used in the iteration for finding the root for drawTheta
// The index of the array is consistent with the index of the summation
void GreensFunction2DRadAbs::makep_mTable(DoubleVector& p_mTable, double r, double t) const
{
   p_mTable.clear();

   double p_0(p_m(0, r, t)); // This is the p_m where m is 0, for the denominator
   p_mTable.push_back(p_0);               // put it in the table

   double p_1(p_m(1, r, t) / p_0);
   p_mTable.push_back(p_1);               // put the first result in the table

   if (p_1 == 0)
   {
      return; // apparently all the terms are zero? We are finished
      // TODO: is this assumption correct??
   }

   double threshold(fabs(GfCfg.EPSILON * p_1)); // get a measure for the allowed error, is this correct?

   double p_m_abs(fabs(p_1));
   double p_m_prev_abs;
   uint m(1);
   do
   {
      m++;
      if (m >= GfCfg.MAX_ORDER()) // If the number of terms is too large
      {
         _log.warn() << "p_m didn't converge (m= " << m << ", t=" << t << ", r0=" << r0_ << ", r=" << r << ", t_est=" << gsl_pow_2(r - r0_) / D_ << ", continuing...";
         break;
      }

      p_m_prev_abs = p_m_abs;                             // store the previous term
      double p_m(this->p_m(m, r, t) / p_0);       // get the next term

      if (!std::isfinite(p_m))                        // if the calculated value is not valid->exit
      {
         _log.warn() << "makep_mTable: invalid value (p_m=" << p_m << ", m=" << m << ")";
         break;
      }

      p_mTable.push_back(p_m);  // put the result in the table
      p_m_abs = fabs(p_m);      // take the absolute value
   } while (p_m_abs >= threshold || p_m_prev_abs >= threshold || p_m_abs >= p_m_prev_abs);

   // truncate when converged enough.
   // if the current term is smaller than threshold
   // AND the previous term is also smaller than threshold
   // AND the current term is smaller than the previous
}

// This method calculates the constants for the drawTheta method when the particle is at the boundary
double GreensFunction2DRadAbs::dp_m_alpha_at_a(uint n, uint m, double t) const
{
   double alpha(getAlpha(m, n)); // get the n-th root using the besselfunctions of order m

   double alpha_sq(alpha * alpha);
   double msq(m * m);
   double ssq(sigma_ * sigma_);

   double s_Anm(sigma_*alpha);
   double a_Anm(a_*alpha);
   double r0Anm(r0_*alpha);

   double Jm_sAnm(gsl_sf_bessel_Jn(m, s_Anm));
   double Jmp1_sAnm(gsl_sf_bessel_Jn(m + 1, s_Anm));
   double Jm_aAnm(gsl_sf_bessel_Jn(m, a_Anm));
   double Ym_aAnm(gsl_sf_bessel_Yn(m, a_Anm));

   double Jm_r0Anm(gsl_sf_bessel_Jn(m, r0Anm));
   double Ym_r0Anm(gsl_sf_bessel_Yn(m, r0Anm));

   // calculating An,m
   double h_ma(h_ - m / sigma_);
   double rho(h_ma*Jm_sAnm + alpha*Jmp1_sAnm);
   double rho_sq(rho*rho);

   // calculating Bn,m(r')
   double B_n_m_r0(Jm_r0Anm * Ym_aAnm - Ym_r0Anm * Jm_aAnm);

   double A_n_m((alpha_sq * rho_sq * B_n_m_r0) / (rho_sq - (Jm_aAnm*Jm_aAnm)*(h_*h_ + alpha_sq - msq / ssq)));

   // calculating the result
   return A_n_m * exp(-D_*alpha_sq*t);
}

// Makes the sum over n for order m for the constants for the drawtheta Method
double GreensFunction2DRadAbs::dp_m_at_a(int m, double t) const
{
   return  funcSum(std::bind(&GreensFunction2DRadAbs::dp_m_alpha_at_a, this, std::placeholders::_1, m, t), MAX_ALPHA_SEQ, GfCfg.EPSILON);
}

// creates a tables of constants for drawTheta when the particle is at the edge of the domain
void GreensFunction2DRadAbs::makedp_m_at_aTable(DoubleVector& p_mTable, double  t) const
{
   p_mTable.clear();

   double p_0(dp_m_at_a(0, t)); // This is the p_m where m is 0, for the denominator
   p_mTable.push_back(p_0);                  // put it in the table

   double p_1(dp_m_at_a(1, t) / p_0);
   p_mTable.push_back(p_1);                  // put the first result in the table

   if (p_1 == 0) return; // apparently all the terms are zero? We are finished

   double threshold(fabs(GfCfg.EPSILON * p_1)); // get a measure for the allowed error

   double p_m_abs(fabs(p_1));
   double p_m_prev_abs;
   uint m(1);
   do
   {
      m++;
      if (m >= GfCfg.MAX_ORDER()) // If the number of terms is too large
      {
         _log.warn() << "dp_m didn't converge (m=" << m << "), continuing...";
         break;
      }

      p_m_prev_abs = p_m_abs;                             // store the previous term
      double p_m(dp_m_at_a(m, t) / p_0);    // get the next term

      // DEBUG (something to check in the future?)
      // if (p_m_abs == 0) {
      //    std::cerr << "Zero valued term found, but convergence is:" <<
      //        p_mTable[p_mTable.size()-1-1]/p_mTable[p_mTable.size()-2-1];
      // }
      // END DEBUG

      if (!std::isfinite(p_m)) // if the calculated value is not valid->exit
      {
         _log.warn() << "makedp_m_at_aTable: invalid value (p_m=" << p_m << ", m=" << m << ", t=" << t << ", p_0=" << p_0 << ")";
         break;
      }

      p_mTable.push_back(p_m);                              // put the result in the table
      p_m_abs = fabs(p_m);                                  // take the absolute value
   } while (p_m_abs >= threshold || p_m_prev_abs >= threshold || p_m_abs >= p_m_prev_abs);
   // truncate when converged enough.
   // if the current term is smaller than threshold
   // AND the previous term is also smaller than threshold
   // AND the current term is smaller than the previous
}

// This calculates the m-th term of the summation for the drawTheta calculation
// Note that m here starts at 0 and in the equations the sum starts at 1!
double GreensFunction2DRadAbs::ip_theta_n(uint m, double theta, const DoubleVector& p_nTable) const
{
   uint m_p1 = m + 1; // artificial increase of m to make sure m starts at 1
   return p_nTable[m_p1] * sin(m_p1*theta) / m_p1;
}

// calculates the cumulative probability of finding the particle at a certain theta
// It is used by the drawTheta method
// It uses the p_nTable for it to speed things up
double GreensFunction2DRadAbs::ip_theta_table(double theta, const DoubleVector& p_nTable) const
{
   uint maxm(static_cast<uint>(p_nTable.size() - 1)); // get the length of the sum
    // it is shifted one because the first entry should
    // be used (m=0)

   return funcSum_all(std::bind(&GreensFunction2DRadAbs::ip_theta_n, this, std::placeholders::_1, theta, p_nTable), maxm);
}

// function to iterate when drawing the theta
double GreensFunction2DRadAbs::ip_theta_F(double theta, const DoubleVector& p_nTable, double value) const
{
   return theta / (M_PI * 2) + (ip_theta_table(theta, p_nTable) / M_PI) - value;
}

// This method draws a theta given a certain r and time (and initial condition of course) 
double GreensFunction2DRadAbs::drawTheta(double rnd, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, sigma_ <= r0_ && r0_ <= a_);
   THROW_UNLESS(std::invalid_argument, sigma_ <= r && r <= a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   // t == 0 means no move.
   if (t <= GfCfg.T_TYPICAL_2D*GfCfg.EPSILON || D_ == 0 || fabs(r0_ - a_) <= GfCfg.EPSILON*GfCfg.L_TYPICAL_2D || rnd <= GfCfg.EPSILON)
      return 0.0;
   if (r == sigma_) // a reaction has occurred, the angle is irrelevant
      return 0.0;

   // making the tables with constants

   DoubleVector p_mTable;                        // a table with constants to make calculation draws much faster
   if (fabs(r - a_) <= GfCfg.EPSILON*GfCfg.L_TYPICAL_2D)      // If the r is at the outer boundary
      makedp_m_at_aTable(p_mTable, t);  // making the table if particle on the outer boundary
   else
      makep_mTable(p_mTable, r, t);     // making the table of constants for the regular case

   double target = rnd*0.5;
   auto f = [target, p_mTable, this](double theta) { return ip_theta_F(theta, p_mTable, target); };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, 0, M_PI);
   return solver.findRoot(GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction2DRadAbs::drawTheta");
}

std::string GreensFunction2DRadAbs::dump() const
{
   std::ostringstream ss;
   ss << "D=" << D_ << ", sigma=" << sigma_ << ", a=" << a_ << ", kf=" << kf_ << ", r0=" << r0_ << ", h=" << h_;
   return ss.str();
}

// Debug functionality
// Directly outputs probability distribution function value of leaving angle
// for given theta, r and t.
double GreensFunction2DRadAbs::givePDFTheta(double theta, double r, double t) const
{
   // input parameter range checks.
   THROW_UNLESS(std::invalid_argument, sigma_ <= r0_ && r0_ <= a_);
   THROW_UNLESS(std::invalid_argument, sigma_ <= r && r <= a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   // making the tables with constants
   DoubleVector p_mTable;                        // a table with constants to make calculations much faster
   if (fabs(r - a_) <= GfCfg.EPSILON*GfCfg.L_TYPICAL_2D)      // If the r is at the outer boundary
   {
      makedp_m_at_aTable(p_mTable, t);  // making the table if particle on the outer boundary
   }
   else
   {
      makep_mTable(p_mTable, r, t);     // making the table of constants for the regular case
   }

   return ip_theta_F(theta, p_mTable, 0.0);

}

// Output the PDF of r, given time t has passed.
double GreensFunction2DRadAbs::givePDFR(double r, double t) const
{
   if (t == 0.0) return r0_; // if no time has passed

   p_survival(t); // calculate the survival probability at this time
   // this is used as the normalization factor
   // BEWARE!!! This also produces the roots An and therefore
   // SETS THE HIGHEST INDEX -> side effect
   // VERY BAD PROGRAMMING PRACTICE!!
   DoubleVector Y0_aAnTable;
   DoubleVector J0_aAnTable;
   DoubleVector Y0J1J0Y1Table;
   createY0J0Tables(Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, t);

   // Create a struct params with the corrects vars.
   p_int_r_params params = { *this, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, 0.0 };

   // Calculate PDF(r) with these vars
   double PDF(p_int_r_F(r, &params));

   return PDF;
}

void GreensFunction2DRadAbs::dumpRoots(int n) const
{
   //    return alphaTables_[n];
   std::cout << "Roots are: {";

   uint size(static_cast<uint>(alphaTables_.size()));
   DoubleVector& alphaTable(this->alphaTables_[n]);

   for (uint i = 0; i < size; i++) {
      std::cout << alphaTable[i] << ",";
   }
   std::cout << "}.\n";
}

// It is used by the drawTheta method
// It uses the p_nTable for it to speed things up
/*
double GreensFunction2DRadAbs::debug_ip_theta_table( double theta) const
{
const DoubleVector& p_nTable( params->p_nTable );	// table with useful constants

uint maxm( p_nTable.size()-1 );	// get the length of the sum
// it is shifted one because the first entry should
// be used (m=0)

double p( funcSum_all( std::bind( &GreensFunction2DRadAbs::ip_theta_n, this, std::placeholders::_1, theta, p_nTable ), maxm ) );
return p;
}
*/

