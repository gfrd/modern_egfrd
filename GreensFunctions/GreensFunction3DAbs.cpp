#include <cmath>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "helperFunctionsGf.hpp"
#include "freeFunctions.hpp"
#include "SphericalBesselGenerator.hpp"
#include "GreensFunction3DAbs.hpp"
#include "Logger.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

const double GreensFunction3DAbs::MIN_T = 1e-18;
const uint GreensFunction3DAbs::MAX_ALPHA_SEQ = 1005;

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction3DAbs");

// --------------------------------------------------------------------------------------------------------------------------------

GreensFunction3DAbs::GreensFunction3DAbs(double D, double r0, double a) : PairGreensFunction(D, 0., r0, 0.), a_(a)
{
   THROW_UNLESS(std::invalid_argument, a >= 0.0);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::p_survival(double t) const
{
   return p_survival_nocollision(t, r0_, D_, a_);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::dp_survival(double t) const
{
   return dp_survival_nocollision(t, r0_, D_, a_);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::p_int_r(double r, double t) const
{
   double Dt = D_ * t;
   double asq = a_ * a_;
   double a_r = 1.0 / a_;
   double asq_r = a_r * a_r;
   double PIr0 = M_PI * r0_;
   double PIr = M_PI * r;
   double r0_angle_factor = PIr0 * a_r;
   double r_angle_factor = PIr * a_r;
   double exp_factor = -Dt * M_PI * M_PI * asq_r;

   uint i_max = std::max(static_cast<uint>(ceil(sqrt(1.0 - asq / M_PI / M_PI * log(GfCfg.TOLERANCE) / Dt))), 2u);

   double p = 0.0;
   for (uint i = 1;; ++i)
   {
      double isq = i * i;
      double term1 = std::exp(exp_factor * isq) * sin(r0_angle_factor * i);
      double term2 = a_ * sin(r_angle_factor * i) - PIr * i * cos(r_angle_factor * i);
      p += term1 * term2 / isq;
      if (i >= i_max) break;
   }
   return p * M_2_PI / PIr0;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::p_n_alpha(uint i, uint n, double r, double t) const
{
   double mDt = -D_ * t;
   double aalpha = gsl_sf_bessel_zero_Jnu(n + 0.5, i + 1);
   double alpha = aalpha / a_;
   double term1 = std::exp(mDt * alpha * alpha);

   auto& s = SphericalBesselGenerator::instance();
   double jr = s.j(n, r * alpha);
   double jr0 = s.j(n, r0_ * alpha);
   double ja2 = s.j(n + 1, aalpha);

   double num = jr * jr0;
   double den = ja2 * ja2;
   return term1 * num / den;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::p_n(int n, double r, double t) const
{
   return funcSum(std::bind(&GreensFunction3DAbs::p_n_alpha, this, std::placeholders::_1, n, r, t), MAX_ALPHA_SEQ);
}

// --------------------------------------------------------------------------------------------------------------------------------

void GreensFunction3DAbs::makep_nTable(DoubleVector& p_nTable, double r, double t) const
{
   p_nTable.clear();

   double factor = 1.0 / (2.0 * M_PI * gsl_pow_3(a_));
   double p_0 = p_n(0, r, t) * factor;
   p_nTable.emplace_back(p_0);

   if (p_0 == 0) return;
   double threshold = std::fabs(p_0 * GfCfg.THETA_TOLERANCE * 1e-1);

   double p_n_prev_abs = std::fabs(p_0);
   for (uint n = 1;;)
   {
      double pn = p_n(n, r, t) * factor;
      if (!std::isfinite(pn)) { _log.error() << "makep_nTable: invalid value: " << std::setprecision(16) << pn << " (n=" << n << ")"; break; }
      p_nTable.emplace_back(pn);

      double p_n_abs = std::fabs(pn);
      // truncate when converged enough.
      if (p_n_abs <= threshold && p_n_prev_abs <= threshold && p_n_abs <= p_n_prev_abs)
         break;

      n++;
      if (n >= GfCfg.MAX_ORDER()) break;
      p_n_prev_abs = p_n_abs;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

static double p_theta_i(uint n, const DoubleVector& p_nTable, const DoubleVector& lgndTable)
{
   return p_nTable[n] * lgndTable[n] * (2 * n + 1);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::p_theta_table(double theta, const DoubleVector& p_nTable) const
{
   uint tableSize = static_cast<uint>(p_nTable.size());
   DoubleVector lgndTable(tableSize);
   gsl_sf_legendre_Pl_array(tableSize - 1, cos(theta), &lgndTable[0]);
   return funcSum_all(std::bind(&p_theta_i, std::placeholders::_1, p_nTable, lgndTable), tableSize) * sin(theta);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::p_theta(double theta, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, theta >= 0.0 && theta <= M_PI);
   THROW_UNLESS(std::invalid_argument, r >= 0 && r < a_);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0 && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0) return 0.0;

   DoubleVector p_nTable;
   makep_nTable(p_nTable, r, t);
   return p_theta_table(theta, p_nTable);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::ip_theta(double theta, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, theta >= 0.0 && theta <= M_PI);
   THROW_UNLESS(std::invalid_argument, r >= 0 && r < a_);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0 && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0 || theta == 0.0) return 0.0;

   DoubleVector p_nTable;
   makep_nTable(p_nTable, r, t);
   return ip_theta_table(theta, p_nTable);
}

// --------------------------------------------------------------------------------------------------------------------------------

static double ip_theta_i(uint n, const DoubleVector& p_nTable, const DoubleVector& lgndTable1)
{
   // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.
   double lgnd_n_m1 = lgndTable1[n];   // n-1
   double lgnd_n_p1 = lgndTable1[n + 2]; // n+1
   return p_nTable[n] * (lgnd_n_m1 - lgnd_n_p1);// / (1.0 + 2.0 * n);
}

double GreensFunction3DAbs::ip_theta_table(double theta, const DoubleVector& p_nTable) const
{
   const uint tableSize = static_cast<uint>(p_nTable.size());
   DoubleVector pTable;
   pTable.reserve(tableSize);
   const double cos_theta = cos(theta);

   // LgndTable is offset by 1 to incorporate the n=-1 case.
   // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

   DoubleVector lgndTable1(tableSize + 2);
   lgndTable1[0] = 1.0;  // n = -1
   gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable1[1]);
   return funcSum_all(std::bind(&ip_theta_i, std::placeholders::_1, p_nTable, lgndTable1), tableSize);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::dp_n_alpha(uint i, uint n, double t) const
{
   double mDt = -D_ * t;
   double aalpha = gsl_sf_bessel_zero_Jnu(n + 0.5, i + 1);
   double alpha = aalpha / a_;
   double term1 = std::exp(mDt * alpha * alpha) * alpha;

   const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());
   double jr0 = s.j(n, r0_ * alpha);
   double ja2 = s.j(n + 1, aalpha);
   return term1 * jr0 / ja2;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::dp_n(int n, double t) const
{
   return funcSum(std::bind(&GreensFunction3DAbs::dp_n_alpha, this, std::placeholders::_1, n, t), MAX_ALPHA_SEQ);
}

// --------------------------------------------------------------------------------------------------------------------------------

void GreensFunction3DAbs::makedp_nTable(DoubleVector& p_nTable, double t) const
{
   p_nTable.clear();

   double factor = -D_ / (2.0 * M_PI * gsl_pow_3(a_));

   double p_0 = dp_n(0, t) * factor;
   p_nTable.emplace_back(p_0);

   if (p_0 == 0) return;

   double threshold = std::fabs(GfCfg.THETA_TOLERANCE * p_0 * 1e-1);

   double p_n_prev_abs = std::fabs(p_0);
   for (uint n = 1;;)
   {
      double p_n = dp_n(n, t) * factor;
      if (!std::isfinite(p_n)) { _log.error() << "makedp_nTable: invalid value: " << std::setprecision(16) << p_n << " (n=" << n << ")"; break; }
      p_nTable.emplace_back(p_n);

      double p_n_abs = std::fabs(p_n);
      // truncate when converged enough.
      if (p_n_abs <= threshold && p_n_prev_abs <= threshold && p_n_abs <= p_n_prev_abs)
         break;

      n++;
      if (n >= GfCfg.MAX_ORDER()) break;
      p_n_prev_abs = p_n_abs;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::dp_theta(double theta, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, theta >= 0.0 && theta <= M_PI);
   THROW_UNLESS(std::invalid_argument, r >= 0 && r <= a_);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0 && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0) return 0.0;

   DoubleVector p_nTable;
   makedp_nTable(p_nTable, t);
   return p_theta_table(theta, p_nTable);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::idp_theta(double theta, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, theta >= 0.0 && theta <= M_PI);
   THROW_UNLESS(std::invalid_argument, r >= 0 && r <= a_);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0 && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0 || theta == 0.0) return 0.0;

   DoubleVector p_nTable;
   makedp_nTable(p_nTable, t);
   return ip_theta_table(theta, p_nTable);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0.0 && r0_ <= a_);

   if (r0_ == a_ || a_ == 0.0) return 0.0;

   double low = 1e-6;
   double high = 1.0;

   auto f = [rnd, this](double t) { return rnd - p_survival(t); };
   gsl_lambda<decltype(f)> F(f);

   // adjust high and low to make sure that f(low) and f(high) straddle.
   while (GSL_FN_EVAL(&F, high) < 0.0)
   {
      high *= 10;
      _log.info() << "drawTime: adjusting high: " << std::setprecision(16) << high;
      if (std::fabs(high) >= 1e10)
      {
         std::stringstream msg;
         msg << type_name() << ": couldn't adjust high. F(" << high << ")=" << GSL_FN_EVAL(&F, high) << "; " << dump();
         throw std::runtime_error(msg.str());
      }
   }

   double low_value = GSL_FN_EVAL(&F, low);
   while (low_value > 0.0)
   {
      low *= .1;
      double low_value_new = GSL_FN_EVAL(&F, low);
      _log.info() << "drawTime: adjusting low: " << std::setprecision(16) << low << ", F=" << low_value_new;
      if (std::fabs(low) <= MIN_T || std::fabs(low_value - low_value_new) < GfCfg.TOLERANCE)
      {
         _log.warn() << "drawTime: couldn't adjust low. F(" << low << ")=" << low_value_new << ", " << dump();
         return MIN_T;
      }

      low_value = low_value_new;
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(MIN_T, GfCfg.TOLERANCE, "GreensFunction3DAbs: drawTime");
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0.0 && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0) return r0_;

   double target = rnd * p_survival(t);
   auto f = [target, t, this](double r) { return p_int_r(r, t) - target; };
   gsl_lambda<decltype(f)> F(f);

   double low = 0.0;
   double high = a_;

   // No initial range guess, except the negative value check below
   // as evaluation of p_int_r in this GF seems pretty robust.

   double highvalue = GSL_FN_EVAL(&F, high);
   if (highvalue < 0.0)
   {
      _log.info() << "drawR: high value < 0.0 (" << std::setprecision(16) << highvalue << "). returning a (" << a_ << ")";
      return a_;
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-15, GfCfg.TOLERANCE, "GreensFunction3DAbs: drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbs::drawTheta(double rnd, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0.0 && r0_ < a_);
   THROW_UNLESS(std::invalid_argument, r >= 0.0 && r <= a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   // t == 0 means no move.
   if (t == 0.0) return 0.0;

   DoubleVector p_nTable;
   if (r == a_ || r < 0.0)
   {
      makedp_nTable(p_nTable, t);
   }
   else
   {
      makep_nTable(p_nTable, r, t);
   }

   double target = rnd * ip_theta_table(M_PI, p_nTable);
   auto f = [target, &p_nTable, this](double theta) { return ip_theta_table(theta, p_nTable) - target; };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, 0.0, M_PI);
   return solver.findRoot(1e-11, GfCfg.THETA_TOLERANCE, "GreensFunction3DAbs: drawTheta");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction3DAbs::dump() const
{
   std::ostringstream ss;
   ss << std::setprecision(16) << "D=" << D_ << ", r0=" << r0_ << ", a=" << a_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
