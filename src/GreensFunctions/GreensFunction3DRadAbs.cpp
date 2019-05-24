#include <cmath>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <iomanip>
#include <gsl/gsl_sf_legendre.h>
#include "helperFunctionsGf.hpp"
#include "freeFunctions.hpp"
#include "SphericalBesselGenerator.hpp"
#include "GreensFunction3DRadAbs.hpp"
#include "Logger.hpp"
#include <gsl/gsl_sf.h>

// --------------------------------------------------------------------------------------------------------------------------------

const uint GreensFunction3DRadAbs::MAX_ALPHA_SEQ = 2000;
const double GreensFunction3DRadAbs::H = 6.0;                   // a fairly strict criterion for safety.

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction3DRadAbs");

// --------------------------------------------------------------------------------------------------------------------------------

GreensFunction3DRadAbs::GreensFunction3DRadAbs(double D, double kf, double r0, double sigma, double a)
   : PairGreensFunction(D, kf, r0, sigma), a_(a), h_(kf / (4.0 * M_PI * sigma * sigma * D)), hsigmap1_(1.0 + h_ * sigma)
{
   THROW_UNLESS(std::invalid_argument, a >= sigma_);
   clearAlphaTable();
}

// --------------------------------------------------------------------------------------------------------------------------------

// Resets the alpha-tables
void GreensFunction3DRadAbs::clearAlphaTable() const
{
   std::for_each(alphaTables_.begin(), alphaTables_.end(), [](DoubleVector rv) { rv.clear(); });
   alphaOffsetTables_[0] = 0;
   std::fill(alphaOffsetTables_.begin() + 1, alphaOffsetTables_.end(), -1);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::getAlpha0(uint i) const
{
   DoubleVector& aTable(alphaTables_[0]);
   uint oldSize = static_cast<uint>(aTable.size());
   if (oldSize <= i)
   {
      aTable.resize(i + 1);
      for (uint m = oldSize; m <= i; ++m)
         aTable[m] = alpha0_i(m);
   }
   return aTable[i];
}

// --------------------------------------------------------------------------------------------------------------------------------

// This function searches for roots (y=0) of the so-called alpha-function
// (::f_alpha()). It either moves a small search interval along the x-axis to
// check for sign-change (which would indicate a root), and calls the GSL
// root finder, or directly calls the root finder if the spacing between
// roots is found to be converged.
double GreensFunction3DRadAbs::getAlpha(uint n, uint i) const
{
   THROW_UNLESS(std::invalid_argument, n < alphaTables_.size());

   DoubleVector& aTable = alphaTables_[n];
   uint oldSize = static_cast<uint>(aTable.size());
   if (oldSize <= i)
   {
      aTable.resize(i + 1);
      uint offset = alphaOffset(n);
      root_fsolver_wrapper solver;
      for (uint m = oldSize; m <= i; ++m)
         aTable[m] = alpha_i(m + offset, n, solver);
   }
   return aTable[i];
}

// --------------------------------------------------------------------------------------------------------------------------------

// The method evaluates the equation for finding the alphas for given alpha. This
// is needed to find the alpha's at which the expression is zero -> alpha is the root.
double GreensFunction3DRadAbs::f_alpha0(double alpha) const
{
   double alpha_a_m_sigma = alpha * (a_ - sigma_);
   double sin_alpha_a_m_sigma = std::sin(alpha_a_m_sigma);
   double cos_alpha_a_m_sigma = std::cos(alpha_a_m_sigma);
   double term1 = alpha * sigma_ * cos_alpha_a_m_sigma;
   double term2 = hsigmap1_ * sin_alpha_a_m_sigma;
   return term1 + term2;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::f_alpha0_aux(double alpha) const
{
   return (a_ - sigma_) * alpha - std::atan(hsigmap1_ / (sigma_ * alpha));
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::alpha0_i(int i) const
{
   THROW_UNLESS(std::invalid_argument, i >= 0);
   const double e = std::numeric_limits<double>::epsilon();

   double target = i * M_PI + M_PI_2;
   auto f = [target, this](double alpha) { return f_alpha0_aux(alpha) - target; };
   gsl_lambda<decltype(f)> F(f);

   // We know the range of the solution from -Pi/2 <= atan <= Pi/2.
   double interval = M_PI / (a_ - sigma_);
   double low = i * interval + e;
   double high = (i + 1 + 20 * e) * interval;

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(0.0, 1e-15, "GreensFunction3DRadAbs::alpha0_i");
}

// --------------------------------------------------------------------------------------------------------------------------------

void GreensFunction3DRadAbs::updateAlphaTable0(double t) const
{
   DoubleVector& alphaTable_0 = getAlphaTable(0);
   alphaTable_0.clear();
   alphaTable_0.reserve(MAX_ALPHA_SEQ);

   double alpha0_0 = alpha0_i(0);
   alphaTable_0.emplace_back(alpha0_0);

   double Dt = D_ * t;
   double alpha_cutoff = std::sqrt((-std::log(GfCfg.TOLERANCE * 1e-3) / Dt));
   for (uint i = 1; i < MAX_ALPHA_SEQ; ++i)
   {
      double alpha0_i = this->alpha0_i(i);
      alphaTable_0.emplace_back(alpha0_i);
      if (alpha0_i > alpha_cutoff && i >= 10) // make at least 10 terms
         break;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

// f_alpha() Calculates the value of the mathematical function f_alpha(). The
// roots (y=0) of this function are constants in the Green's Functions.
double GreensFunction3DRadAbs::f_alpha(double alpha, int n) const
{
   double aAlpha = a_ * alpha;
   double sigmaAlpha = sigma_ * alpha;
   double hSigma_m_n = h_ * sigma_ - n;
   const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());
   double js1 = s.j(n, sigmaAlpha);
   double ys1 = s.y(n, sigmaAlpha);
   double js2 = s.j(n + 1, sigmaAlpha);
   double ys2 = s.y(n + 1, sigmaAlpha);
   double ja = s.j(n, aAlpha);
   double ya = s.y(n, aAlpha);
   double term1 = (hSigma_m_n * js1 + sigmaAlpha * js2) * ya;
   double term2 = (hSigma_m_n * ys1 + sigmaAlpha * ys2) * ja;
   double factor(2.0 * alpha * std::sqrt(a_ * sigma_) * M_1_PI);
   double result = ((term1 - term2) * factor);
   return result;
}

// --------------------------------------------------------------------------------------------------------------------------------

double G(uint n, uint k)
{
   return gsl_sf_fact(n + k) / (gsl_sf_fact(k) * gsl_sf_fact(n - k));
}

// --------------------------------------------------------------------------------------------------------------------------------

double P(int n, double x)
{
   int sign = 1;
   double sx2 = 1.0;
   double x2sq_r = 1.0 / gsl_pow_2(x + x);
   double result = 0.0;
   uint maxm = n / 2;
   for (uint m = 0; m <= maxm; ++m)
   {
      result += sign * sx2 * G(n, 2 * m);
      sign = -sign;
      sx2 *= x2sq_r;
   }
   return result;
}

DoublePair P2(int n, double x)
{
   int sign = 1;
   double sx2 = 1.0;
   double x2sq_r = 1.0 / gsl_pow_2(x + x);
   double result = 0.0;
   double resultp = 0.0;
   uint maxm = n / 2;
   for (uint m = 0; m <= maxm; ++m)
   {
      double sx2p = sign * sx2;
      result += sx2p * G(n, 2 * m);
      resultp += sx2p * G(n + 1, 2 * m);
      sign = -sign;
      sx2 *= x2sq_r;
   }

   if (n % 2) resultp += sign * sx2 * G(n + 1, n + 1);

   return std::make_pair(result, resultp);
}

// --------------------------------------------------------------------------------------------------------------------------------

double Q(int n, double x)
{
   int sign = 1;
   double sx2 = 1.0 / (x + x);
   double x2sq = sx2 * sx2;
   double result = 0.0;
   uint maxm = (n + 1) / 2; // sum_(0)^((n-1)/2)
   for (uint m = 0; m < maxm; ++m)
   {
      result += sign * sx2 * G(n, 2 * m + 1);;
      sign = -sign;  // (-1)^m
      sx2 *= x2sq;
   }
   return result;
}

DoublePair Q2(int n, double x)
{
   int sign = 1;  // (-1)^m
   double sx2 = 1.0 / (x + x);
   double x2sq = sx2 * sx2;
   double result = 0.0;
   double resultp = 0.0;
   uint maxm = (n + 1) / 2; // sum_(0)^((n-1)/2)
   for (uint m = 0; m < maxm; ++m)
   {
      double sx2p = sign * sx2;
      uint m2p1 = 2 * m + 1;
      result += sx2p * G(n, m2p1);
      resultp += sx2p * G(n + 1, m2p1);

      sign = -sign; // (-1)^m
      sx2 *= x2sq;
   }

   if (!(n % 2)) resultp += sign * sx2 * G(n + 1, n + 1);

   return std::make_pair(result, resultp);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::f_alpha_aux(double alpha, int n) const
{
   if (alpha == 0.0) return -1.0;

   double aAlpha = a_ * alpha;
   double sigmaAlpha = sigma_ * alpha;
   double n_m_hSigma = n - h_ * sigma_;

   /*(a - s) u - ArcTan[(P[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) -
     Q[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))/
     (Q[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) +
     P[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))]
     */

   double Pa = P(n, aAlpha);
   double Qa = Q(n, aAlpha);

   double Ps, Psp;
   std::tie(Ps, Psp) = P2(n, sigmaAlpha);
   double Qs, Qsp;
   std::tie(Qs, Qsp) = Q2(n, sigmaAlpha);

   double Qa_Pa = Qa / Pa;
   double A = sigmaAlpha * Qsp - n_m_hSigma * Ps;
   double B = sigmaAlpha * Psp + n_m_hSigma * Qs;

   // this form, dividing all terms by Pa, prevents overflow.
   double angle = (A - Qa_Pa * B) / (Qa_Pa * A + B);
   double result = (a_ - sigma_) * alpha - std::atan(angle);

   //printf("f_alpha_aux: alpha=%.16g, n=%8d => result=%.16g\r\n", alpha, n, result);
   return result;

}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::alpha_i(int i, int n, const root_fsolver_wrapper& solver) const
{
   double target = M_PI * i + M_PI_2;
   double factor = 1.0 / (a_ - sigma_);
   double low = (target - M_PI_2) * factor;
   double high = (target + M_PI_2) * factor;

   auto f = [n, target, this](double alpha) { return f_alpha_aux(alpha, n) - target; };
   gsl_lambda<decltype(f)> F(f);
   solver.set(&F, low, high);
   return solver.findRoot(1e-6, 1e-15, "GreensFunction3DRadAbs::alpha_i");
}

// --------------------------------------------------------------------------------------------------------------------------------

uint GreensFunction3DRadAbs::alphaOffset(uint n) const
{
   THROW_UNLESS(std::invalid_argument, n < alphaOffsetTables_.size());

   if (alphaOffsetTables_[n] >= 0)
      return alphaOffsetTables_[n];

   uint offset = alphaOffsetTables_[n - 1];
   double factor = 1.0 / (a_ - sigma_);

   double target = (offset * M_PI + M_PI_2);
   // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
   double alphaMid = target * factor;
   double alphaHalfRange = M_PI_2 * factor;
   double low = alphaMid - alphaHalfRange * (1.0 - 1e-3); // avoid zero.
   double high = alphaMid + alphaHalfRange;

   // Here we find the interval where the first positive root is in.
   // We find the first pair of alpha (Pi * offset + Pi/2) +- Pi/2 / (a - sigma)
   // where the values of f_alpha() straddle.
   // The assumption is the interval between roots is not much smaller than Pi / (a - sigma).

   double lowvalue = f_alpha(low, n);
   double highvalue = f_alpha(high, n);

   for (;;)
   {
      if (lowvalue * highvalue < 0) // low and high straddle?
         break;

      ++offset;
      target = M_PI * offset + M_PI_2;
      low = (target - M_PI_2) * factor;
      high = (target + M_PI_2) * factor;

      lowvalue = highvalue;
      highvalue = f_alpha(high, n);
   }

   alphaOffsetTables_[n] = offset;
   return offset;
}

// --------------------------------------------------------------------------------------------------------------------------------

void GreensFunction3DRadAbs::updateAlphaTable(uint n, double t) const
{
   THROW_UNLESS(std::invalid_argument, n < GfCfg.MAX_ORDER());

   if (n == 0)
   {
      updateAlphaTable0(t);
      return;
   }

   uint offset = alphaOffset(n);
   DoubleVector& alphaTable_n = getAlphaTable(n);
   alphaTable_n.clear();
   alphaTable_n.reserve(MAX_ALPHA_SEQ);

   root_fsolver_wrapper solver;
   double alphan0 = alpha_i(offset, n, solver);
   alphaTable_n.emplace_back(alphan0);

   double Dt = D_ * t;
   double alphan0sq = alphan0 * alphan0;
   double threshold = GfCfg.TOLERANCE * 1e-2 * alphan0sq * std::exp(-Dt * alphan0sq);

   uint end = offset + MAX_ALPHA_SEQ;
   for (uint i = offset + 1; ; ++i)
   {
      if (i >= end) { _log.info() << "alphaTable (" << n << "): didn't converge. t=" << std::scientific << std::setprecision(16) << t << ", " << dump().c_str(); break; }

      double alphai = alpha_i(i, n, solver);
      alphaTable_n.emplace_back(alphai);

      // cutoff
      double alphaisq = alphai * alphai;
      if (alphaisq * std::exp(-Dt * alphaisq) < threshold) break;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

// calculates the constant part of the i-th term for the survival probability
double GreensFunction3DRadAbs::p_survival_i(double alpha) const
{
   double sigmasq = sigma_ * sigma_;
   double alphasq = alpha * alpha;
   double denominator = (r0_ * hsigmap1_ * alpha * (-hsigmap1_ * (a_ + a_ * h_ * sigma_ - h_ * sigmasq) + (sigma_ - a_) * sigmasq * alphasq));
   return -2.0 * (h_ * sigmasq * hsigmap1_ - a_ * (hsigmap1_ * hsigmap1_ + sigmasq * alphasq) * std::cos(alpha * (a_ - sigma_))) * num_r0(alpha) / denominator;
}

// --------------------------------------------------------------------------------------------------------------------------------

// Calculates the n-th term of the summation for calculating the flux through  the outer interface (escape)
double GreensFunction3DRadAbs::leavea_i(double alpha) const
{
   double sigmasq = sigma_ * sigma_;
   double alphasq = alpha * alpha;
   double numerator = alpha * (hsigmap1_ * hsigmap1_ + sigmasq * alphasq) * std::cos(alpha * (a_ - sigma_));
   double denominator = 2 * a_ * M_PI * r0_ * hsigmap1_ * (hsigmap1_ * (a_ + a_ * h_ * sigma_ - h_ * sigmasq) + (a_ - sigma_) * sigmasq * alphasq);
   return D_ * numerator * num_r0(alpha) / denominator;
}

// --------------------------------------------------------------------------------------------------------------------------------

// Calculates the n-th term of the summation for calculating the flux through  the inner interface (reaction)
double GreensFunction3DRadAbs::leaves_i(double alpha) const
{
   double sigmasq = sigma_ * sigma_;
   double alphasq = alpha * alpha;
   double numerator = h_ * alpha * num_r0(alpha);
   double denominator = 2 * M_PI * r0_ * ((a_ - sigma_) * sigmasq * alphasq + hsigmap1_ * (a_ + a_ * h_ * sigma_ - h_ * sigmasq));
   return -D_ * numerator / denominator;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::num_r0(double alpha) const
{
   double angle_r0 = alpha * (r0_ - sigma_);
   return alpha * sigma_ * std::cos(angle_r0) + hsigmap1_ * std::sin(angle_r0);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::p_int_r_i(double r, double alpha, double num_r0) const
{
   double angle_r = alpha * (r - sigma_);
   double sigmasq = sigma_ * sigma_;
   double alphasq = alpha * alpha;
   double hsigma = h_ * sigma_;
   double numerator = alpha * (hsigma * sigma_ - hsigma * r * std::cos(angle_r) - (r - sigma_) * std::cos(angle_r)) + (hsigmap1_ + r * sigma_ * alphasq) * std::sin(angle_r);
   double denominator = r0_ * alphasq * ((a_ - sigma_) * sigmasq * alphasq + hsigmap1_ * (a_ + a_ * h_ * sigma_ - h_ * sigmasq));
   return 2.0 * numerator * num_r0 / denominator;
}

// --------------------------------------------------------------------------------------------------------------------------------

// calculates a table with all the constant factors for the survival probability
void GreensFunction3DRadAbs::createPsurvTable(DoubleVector& table) const
{
   const DoubleVector& alphaTable_0 = getAlphaTable(0);
   table.clear();
   table.reserve(alphaTable_0.size());
   std::transform(alphaTable_0.begin(), alphaTable_0.end(), std::back_inserter(table), std::bind(&GreensFunction3DRadAbs::p_survival_i, this, std::placeholders::_1));
}

// --------------------------------------------------------------------------------------------------------------------------------

uint GreensFunction3DRadAbs::guess_maxi(double t) const
{
   const uint safety = 2;
   if (!std::isfinite(t)) return safety;

   double alpha0 = getAlpha0(0);
   double Dt = D_ * t;
   double threshold = std::exp(-Dt * alpha0 * alpha0) * GfCfg.TOLERANCE * 0.1;
   if (threshold <= 0.0) return MAX_ALPHA_SEQ;

   double max_alpha = std::sqrt(alpha0 * alpha0 - std::log(threshold) / Dt);
   return std::min(safety + static_cast<uint>(max_alpha * (a_ - sigma_) / M_PI), MAX_ALPHA_SEQ);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::p_survival_table(double t, DoubleVector& psurvTable) const
{
   double max_dist = H * std::sqrt(6.0 * D_ * t);
   double a_dist = a_ - r0_;
   double s_dist = r0_ - sigma_;

   if (a_dist > max_dist)
   {
      if (s_dist > max_dist) // far from anything; it'll survive.
         return 1.0;

      // close only to s, ignore a
      return p_survival_irr(t, r0_, kf_, D_, sigma_);
   }

   if (s_dist > max_dist)  // close only to a.
      return p_survival_nocollision(t, r0_, D_, a_);

   // close to both boundaries.  do the normal calculation.
   const uint i_max = guess_maxi(t);
   if (psurvTable.size() < i_max)
   {
      getAlpha0(i_max);  // this updates the table
      createPsurvTable(psurvTable);
   }
   return funcSum_all(std::bind(&GreensFunction3DRadAbs::p_survival_i_exp_table, this, std::placeholders::_1, t, psurvTable), i_max);
}

// --------------------------------------------------------------------------------------------------------------------------------

// calculates the ith term with exponent and time for the survival probability
double GreensFunction3DRadAbs::p_survival_i_exp_table(uint i, double t, const DoubleVector& table) const
{
   double alpha = getAlpha0(i);
   return std::exp(-D_ * t * alpha * alpha) * table[i];
}

// --------------------------------------------------------------------------------------------------------------------------------

// calculates the flux leaving through the inner interface at a given moment
// FIXME: This is inaccurate for small t's!!
double GreensFunction3DRadAbs::leaves(double t) const
{
   return funcSum(std::bind(&GreensFunction3DRadAbs::leaves_i_exp, this, std::placeholders::_1, t), MAX_ALPHA_SEQ);
}

// adds the exponential with the time to the sum. Needed for the inner interface (reaction)
double GreensFunction3DRadAbs::leaves_i_exp(uint i, double t) const
{
   double alpha = getAlpha0(i);
   return std::exp(-D_ * t * alpha * alpha) * leaves_i(alpha);
}

// --------------------------------------------------------------------------------------------------------------------------------

// calculates the flux leaving through the outer interface at a given moment
double GreensFunction3DRadAbs::leavea(double t) const
{
   return funcSum(std::bind(&GreensFunction3DRadAbs::leavea_i_exp, this, std::placeholders::_1, t), MAX_ALPHA_SEQ);
}

// adds the exponential with the time to the sum. Needed for the calculation of the flux through the outer interface
double GreensFunction3DRadAbs::leavea_i_exp(uint i, double t) const
{
   double alpha = getAlpha0(i);
   return std::exp(-D_ * t * alpha * alpha) * leavea_i(alpha);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::p_int_r(double r, double t) const
{
   return funcSum(std::bind(&GreensFunction3DRadAbs::p_int_r_i_exp, this, std::placeholders::_1, t, r), MAX_ALPHA_SEQ);
}

double GreensFunction3DRadAbs::p_int_r_i_exp(uint i, double t, double r) const
{
   double alpha = getAlpha0(i);
   return std::exp(-D_ * t * alpha * alpha) * p_int_r_i(r, alpha, num_r0(alpha));
}

// --------------------------------------------------------------------------------------------------------------------------------

// Draws a first passage time, this could be an escape (through the outer boundary) or a reaction (through
// the inner boundary)
double GreensFunction3DRadAbs::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_ && r0_ <= a_);

   if (r0_ == a_ || a_ == sigma_) return 0.0;

   double dist = kf_ != 0 ? std::min(a_ - r0_, r0_ - sigma_) : a_ - r0_;
   double t_guess = 0.1 * dist * dist / (6.0 * D_);

   DoubleVector table;
   auto f = [rnd, &table, this](double t) { return rnd - p_survival_table(t, table); };
   gsl_lambda<decltype(f)> F(f);

   double low = t_guess;
   double high = t_guess;
   double minT = std::min(sigma_ * sigma_ / D_ * GfCfg.MIN_T_FACTOR, t_guess * 1e-6);
   double value = GSL_FN_EVAL(&F, t_guess);
   if (value < 0.0)
   {
      for (;;)
      {
         high *= 10;
         value = GSL_FN_EVAL(&F, high);

         if (value >= 0.0) break;

         if (std::fabs(high) >= 1e10)
         {
            std::stringstream msg;
            msg << type_name() << ": couldn't adjust high. F(" << high << ")=" << value << "; " << dump();
            throw std::runtime_error(msg.str());
         }
      }
   }
   else
   {
      for (double prev_value = value;; prev_value = value)
      {
         low *= .1;
         value = GSL_FN_EVAL(&F, low);

         if (value <= 0.0) break;

         if (std::fabs(low) <= minT || std::fabs(value - prev_value) < GfCfg.TOLERANCE)
         {
            _log.warn() << "drawTime: couldn't adjust low. F(" << low << ")=" << value << ", " << dump();
            return low;
         }
      }
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(0.0, GfCfg.TOLERANCE, "GreensFuntion3DRadAbs::drawTime");
}

// --------------------------------------------------------------------------------------------------------------------------------

// This determines based on the flux at a certain time, if the 'escape' was a reaction or a proper escape
GreensFunction::EventKind GreensFunction3DRadAbs::drawEventType(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_ && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t > 0.0);

   if (kf_ == 0) return EventKind::IV_ESCAPE;

   // First, check if r0 is close only either to a or sigma relative to Dt.  In such cases, the event type is always IV_ESCAPE or 
   // IV_REACTION, respectively. This avoids numerical instability in calculating leavea() and/or leaves().

   double max_dist = H * std::sqrt(6.0 * D_ * t);
   double a_dist = a_ - r0_;
   double s_dist = r0_ - sigma_;

   if (a_dist > max_dist)
   {
      if (s_dist < max_dist)
         return EventKind::IV_REACTION;
   }
   else
   {
      if (s_dist > max_dist)
         return EventKind::IV_ESCAPE;
   }

   double reaction = leaves(t) * 4.0 * M_PI * sigma_ * sigma_;
   double escape = leavea(t) * 4.0 * M_PI * a_ * a_;
   double value = reaction / (reaction + escape);
   return rnd <= value ? EventKind::IV_REACTION : EventKind::IV_ESCAPE;
}

// --------------------------------------------------------------------------------------------------------------------------------

// This draws a radius R at a given time, provided that the particle was at r0 at t=0
double GreensFunction3DRadAbs::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_ && r0_ < a_);

   if (t == 0.0) return r0_;

   DoubleVector table;
   double psurv = p_survival_table(t, table);

   const double target = psurv * rnd;
   auto f = [t, target, this](double r) { return p_int_r(r, t) - target; };
   auto F = gsl_lambda<decltype(f)>(f);

   double low = r0_;
   double high = r0_;
   double sqrt6Dt = std::sqrt(6.0 * D_ * t);
   if (GSL_FN_EVAL(&F, r0_) < 0.0)
   {
      for (uint hh = 3;; ++hh)
      {
         high = r0_ + hh * sqrt6Dt;
         if (high > a_)
         {
            if (GSL_FN_EVAL(&F, a_) < 0.0)
            {
               _log.info() << "drawR: p_int_r(a) < 0.0. returning a";
               return a_;
            }

            high = a_;
            break;
         }

         double value = GSL_FN_EVAL(&F, high);
         if (value > 0.0)
            break;
      }

   }
   else
   {
      for (uint hh = 3;; ++hh)
      {
         low = r0_ - hh * sqrt6Dt;
         if (low < sigma_)
         {
            if (GSL_FN_EVAL(&F, sigma_) > 0.0)
            {
               _log.info() << "drawR: p_int_r(sigma) > 0.0. returning sigma";
               return sigma_;
            }

            low = sigma_;
            break;
         }

         if (GSL_FN_EVAL(&F, low) < 0.0) break;
      }
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-15, GfCfg.TOLERANCE, "GreensFuntion3DRadAbs::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::p_n_alpha(uint i, uint n, double r, double t) const
{
   double alpha = getAlpha(n, i);
   double alphasq = alpha * alpha;

   const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());
   double js1 = s.j(n, sigma_ * alpha);
   double js2 = s.j(n + 1, sigma_ * alpha);
   double ja = s.j(n, a_ * alpha);
   double ya = s.y(n, a_ * alpha);
   double jr = s.j(n, r * alpha);
   double yr = s.y(n, r * alpha);
   double jr0 = s.j(n, r0_ * alpha);
   double yr0 = s.y(n, r0_ * alpha);
   double J = (h_ * sigma_ - n) * js1 + sigma_ * alpha * js2;
   double Jsq = J * J;
   double denominator = a_ * (n + n * n - sigma_ * (h_ + h_ * h_ * sigma_ + sigma_ * alphasq)) * ja * ja + sigma_ * Jsq;
   return alphasq * alphasq * std::exp(-D_ * t * alphasq) * (Jsq * (ja * yr - ya * jr) * (ja * yr0 - ya * jr0)) / denominator;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::dp_n_alpha_at_a(uint i, uint n, double t) const
{
   double alpha = getAlpha(n, i);
   double alphasq = alpha * alpha;
   const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());
   double js1 = s.j(n, sigma_ * alpha);
   double js2 = s.j(n + 1, sigma_ * alpha);
   double ja = s.j(n, a_ * alpha);
   double ya = s.y(n, a_ * alpha);
   double jr0 = s.j(n, r0_ * alpha);
   double yr0 = s.y(n, r0_ * alpha);
   double J = (h_ * sigma_ - n) * js1 + sigma_ * alpha * js2;
   double Jsq = J * J;
   double term1 = (exp(-D_ * t * alphasq) * (alphasq * alpha));
   double den1 = (a_ * (n + n * n - sigma_ * (h_ + h_ * h_ * sigma_ + sigma_ * alphasq)) * ja * ja);
   double den2 = (sigma_ * Jsq);
   double den = (den1 + den2);
   double JY(-jr0 * ya + ja * yr0);
   double num = (Jsq * JY);
   double result = (term1 * num / den);
   return result;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::funcSumMaxAlpha(int n, double max_alpha, const std::function<double(uint, uint)>& func) const
{
   const uint min_i = 2;
   double p = 0.0;
   uint i;
   for (i = 0; i < MAX_ALPHA_SEQ; ++i)
   {
      p += func(i, n);;
      if (getAlpha(n, i) >= max_alpha && i >= min_i) break;
   }
   return p;
}

void GreensFunction3DRadAbs::make_pn_table(DoubleVector& p_nTable, double t, double factor, const std::function<double(uint, uint)>& func, bool max_Dt) const
{
   p_nTable.clear();

   double alpha00 = getAlpha(0, 0);
   double max_alpha = std::sqrt((max_Dt ? D_ * t : 1.0) * alpha00 * alpha00 - std::log(GfCfg.THETA_TOLERANCE * 0.1) / (D_ * t));           // factor D_*t, was not included in make_p_n_table, but was in make_dp_n_at_a_table
   double p_0 = funcSumMaxAlpha(0, max_alpha, func) * factor;
   p_nTable.emplace_back(p_0);

   if (p_0 == 0) return;

   double threshold = std::fabs(GfCfg.THETA_TOLERANCE * p_0);
   double p_n_prev_abs = std::fabs(p_0);
   for (uint n = 1; n <= GfCfg.MAX_ORDER(); ++n)
   {
      if (getAlpha(n, 0) >= max_alpha) break;

      double p_n = funcSumMaxAlpha(n, max_alpha, func) * factor;
      p_nTable.emplace_back(p_n);

      // truncate when converged enough.
      double p_n_abs = std::fabs(p_n);
      if (p_n_abs < threshold && p_n_prev_abs < threshold && p_n_abs <= p_n_prev_abs)
         break;

      p_n_prev_abs = p_n_abs;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

static double ip_theta_n(uint n, const DoubleVector& p_nTable, const DoubleVector& lgndTable1)
{
   // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.
   // the term (1 + 2 n) is canceled out.
   return p_nTable[n] * (lgndTable1[n] - lgndTable1[n + 2]);
}

static double ip_theta_table(double theta, const DoubleVector& p_nTable)
{
   uint tableSize = static_cast<uint>(p_nTable.size());

   // LgndTable is offset by 1 to incorporate the n=-1 case.
   // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

   DoubleVector lgndTable(tableSize + 2);
   lgndTable[0] = 1.0;  // n = -1
   gsl_sf_legendre_Pl_array(tableSize, std::cos(theta), &lgndTable[1]);

   return funcSum_all(std::bind(&ip_theta_n, std::placeholders::_1, p_nTable, lgndTable), tableSize);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadAbs::ip_theta(double theta, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, theta >= 0.0 && theta <= M_PI);
   THROW_UNLESS(std::invalid_argument, r >= sigma_ && r < a_);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_ && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0 || theta == 0.0) return 0.0;
   DoubleVector p_nTable;
   // makep_nTable
   make_pn_table(p_nTable, t, a_ * sigma_ / (M_PI * 2), std::bind(&GreensFunction3DRadAbs::p_n_alpha, this, std::placeholders::_1, std::placeholders::_2, r, t));
   return ip_theta_table(theta, p_nTable);
}

// --------------------------------------------------------------------------------------------------------------------------------

// This method draws a theta given a certain r and time (and initial condition of course)
double GreensFunction3DRadAbs::drawTheta(double rnd, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_ && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, r >= sigma_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0) return 0.0;

   DoubleVector p_nTable;
   if (r >= a_)
      //makedp_n_at_aTable(p_nTable, t);
      make_pn_table(p_nTable, t, D_ * sigma_ / (a_ * M_PI * 2), std::bind(&GreensFunction3DRadAbs::dp_n_alpha_at_a, this, std::placeholders::_1, std::placeholders::_2, t), true);
   else
      //makep_nTable(p_nTable, r, t);
      make_pn_table(p_nTable, t, a_ * sigma_ / (M_PI * 2), std::bind(&GreensFunction3DRadAbs::p_n_alpha, this, std::placeholders::_1, std::placeholders::_2, r, t));

   double ip_theta_pi_rnd = rnd * ip_theta_table(M_PI, p_nTable);
   auto f = [p_nTable, ip_theta_pi_rnd, this](double theta) { return ip_theta_table(theta, p_nTable) - ip_theta_pi_rnd; };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, 0.0, M_PI);
   return solver.findRoot(1e-11, GfCfg.THETA_TOLERANCE, "GreensFunction3DRadAbs::drawTheta");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction3DRadAbs::dump() const
{
   std::ostringstream ss;
   ss << std::scientific << std::setprecision(16) << "D=" << D_ << ", kf=" << kf_ << ", r0=" << r0_ << ", sigma=" << sigma_ << ", a=" << a_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
