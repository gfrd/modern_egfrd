#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "freeFunctions.hpp"
#include "helperFunctionsGf.hpp"
#include "SphericalBesselGenerator.hpp"
#include "GreensFunction3DRadInf.hpp"
#include "Logger.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction3DRadInf");

// --------------------------------------------------------------------------------------------------------------------------------

GreensFunction3DRadInf::GreensFunction3DRadInf(double D, double kf, double r0, double sigma)
   : PairGreensFunction(D, kf, r0, sigma), kD_(4.0 * M_PI * sigma * D), alpha_((1.0 + (kf / kD_)) * (std::sqrt(D) / sigma))
{
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadInf::p_corr_R(double alpha, uint n, double r, double t) const
{
   double ks = kf_ * sigma_;
   double ks_m_n = ks - n;
   double alphasq = alpha * alpha;
   double term1 = std::exp(-D_ * t * alphasq);
   double sAlpha = sigma_ * alpha;
   double rAlpha = r * alpha;
   double r0Alpha = r0_ * alpha;

   const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());
   double js = s.j(n, sAlpha);
   double ys = s.y(n, sAlpha);
   double js1 = s.j(n + 1, sAlpha);
   double ys1 = (s.y(n + 1, sAlpha));
   double jr = s.j(n, rAlpha);
   double yr = s.y(n, rAlpha);
   double jr0 = s.j(n, r0Alpha);
   double yr0 = s.y(n, r0Alpha);

   double R1 = ks_m_n * js + sAlpha * js1;
   double R2 = ks_m_n * ys + sAlpha * ys1;
   double F1R1 = R1 * jr * jr0 - R1 * yr * yr0;
   double F2 = jr0 * yr + jr * yr0;

   double num = 2.0 * std::sqrt(r * r0_) * alphasq * R1 * (F1R1 + F2 * R2);
   double den = M_PI * (R1 * R1 + R2 * R2);
   double result = term1 * num / den;
   ASSERT(std::isfinite(result));
   return result;
}

// --------------------------------------------------------------------------------------------------------------------------------

//double GreensFunction3DRadInf::p_corr(double theta, double r, double t) const
//{
//   DoubleVector RnTable;
//   makeRnTable(RnTable, r, t);
//   return p_corr_table(theta, r, RnTable);
//}
//
//double GreensFunction3DRadInf::ip_corr(double theta, double r, double t) const
//{
//   DoubleVector RnTable;
//   makeRnTable(RnTable, r, t);
//   return ip_corr_table(theta, r, RnTable);
//}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadInf::p_int_r(double r, double t) const
{
   double Dt = D_ * t;
   double kf_kD = kf_ + kD_;
   double sqrtDt4 = std::sqrt(4.0 * Dt);
   double alphasqrtt = alpha_ * std::sqrt(t);
   double r_r0_2s_sqrtDt4 = (r - 2.0 * sigma_ + r0_) / sqrtDt4;
   double r_r0_sqrtDt4 = (r - r0_) / sqrtDt4;
   double r0_s_sqrtDt4 = (r0_ - sigma_) / sqrtDt4;
   double term1 = (std::expm1(-gsl_pow_2(r_r0_2s_sqrtDt4)) - std::expm1(-gsl_pow_2(r_r0_sqrtDt4))) * std::sqrt(Dt / M_PI);
   double erf_r_r0_2s_sqrtDt4 = std::erf(r_r0_2s_sqrtDt4);
   double term2 = kf_kD * r0_ * std::erf(r_r0_sqrtDt4) + kf_kD * r0_ * erf_r_r0_2s_sqrtDt4 + 2.0 * kf_ * sigma_ * (std::erf(r0_s_sqrtDt4) - erf_r_r0_2s_sqrtDt4);
   double term3 = kf_ * sigma_ * W(r0_s_sqrtDt4, alphasqrtt) - (kf_ * r + kD_ * (r - sigma_)) * W(r_r0_2s_sqrtDt4, alphasqrtt);
   return 1 / r0_ * (term1 + 1 / kf_kD * (0.5 * term2 + term3));
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadInf::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_);

   double low = 1e-100;
   double high = 100;

   {
      double maxp = __p_reaction_irr(INFINITY, r0_, kf_, D_, sigma_, alpha_, kD_);
      if (rnd >= maxp) return INFINITY;
   }

   auto f = [rnd, this](double t) { return __p_reaction_irr(t, r0_, kf_, D_, sigma_, alpha_, kD_) - rnd; };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-18, 1e-12, "GreensFunction3DRadInf::drawTime");

}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadInf::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0) return r0_;

   double target = rnd * (1.0 - __p_reaction_irr(t, r0_, kf_, D_, sigma_, alpha_, kD_));
   auto f = [t, target, this](double r) { return p_int_r(r, t) - target; };
   gsl_lambda<decltype(f)> F(f);

   double low = r0_;
   double high = r0_;
   double sqrt6Dt = std::sqrt(6.0 * D_ * t);
   if (GSL_FN_EVAL(&F, r0_) < 0.0)
   {
      for (uint H = 3;; ++H)
      {
         if (H > 20) throw std::runtime_error("GreensFunction3DRadInf: drawR: H > 20 while adjusting upper bound of r");

         high = r0_ + H * sqrt6Dt;
         double value = GSL_FN_EVAL(&F, high);
         if (value > 0.0) break;
      }
   }
   else
   {
      for (uint H = 3;; ++H)
      {
         low = r0_ - H * sqrt6Dt;
         if (low < sigma_)
         {
            if (GSL_FN_EVAL(&F, sigma_) > 0.0)
            {
               _log.info() << "drawR: p_int_r(sigma) > 0.0. " "returning sigma.";
               return sigma_;
            }

            low = sigma_;
            break;
         }

         double value(GSL_FN_EVAL(&F, low));
         if (value < 0.0) break;
      }
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-15, GfCfg.TOLERANCE, "GreensFunction3DRadInf::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadInf::Rn(uint n, double r, double t, gsl_integration_workspace* workspace, double tol) const
{
   double integral, error;
   auto f = [n, r, t, this](double alpha) { return p_corr_R(alpha, n, r, t); };
   gsl_lambda<decltype(f)> F(f);

   double umax = std::sqrt(40.0 / (D_ * t));
   gsl_integration_qag(&F, 0.0, umax, tol, GfCfg.THETA_TOLERANCE, 2000, GSL_INTEG_GAUSS61, workspace, &integral, &error);
   return integral;
}

// --------------------------------------------------------------------------------------------------------------------------------

//double GreensFunction3DRadInf::p_corr_table(double theta, double r, const DoubleVector& RnTable) const
//{
//   const uint tableSize = static_cast<uint>(RnTable.size());
//   if (tableSize == 0) return 0.0;
//
//   DoubleVector lgndTable(tableSize);
//   gsl_sf_legendre_Pl_array(tableSize - 1, std::cos(theta), &lgndTable[0]);
//
//   auto f = [&RnTable, &lgndTable](uint n) {return RnTable[n] * lgndTable[n] * (2.0 * n + 1.0); };
//   double p = funcSum_all(f, tableSize);
//   double result = -p * std::sin(theta);
//   result /= 4.0 * M_PI * std::sqrt(r * r0_);
//   return result;
//}

double GreensFunction3DRadInf::ip_corr_table(double theta, double r, const DoubleVector& RnTable) const
{
   const uint tableSize = static_cast<uint>(RnTable.size());
   if (tableSize == 0) return 0.0;

   DoubleVector lgndTable(tableSize + 2);
   lgndTable[0] = 1.0; // n = -1
   gsl_sf_legendre_Pl_array(tableSize, std::cos(theta), &lgndTable[1]);

   auto f = [&RnTable, &lgndTable](uint n) {return RnTable[n] * (lgndTable[n] - lgndTable[n + 2]); };
   double p = funcSum_all(f, tableSize);
   return -p / (4.0 * M_PI * std::sqrt(r * r0_));
}

// --------------------------------------------------------------------------------------------------------------------------------

//double GreensFunction3DRadInf::p_theta(double theta, double r, double t) const
//{
//   DoubleVector RnTable;
//   makeRnTable(RnTable, r, t);
//   return p_theta_free(theta, r, r0_, t, D_) + p_corr_table(theta, r, RnTable);
//}
//
//double GreensFunction3DRadInf::ip_theta(double theta, double r, double t) const
//{
//   DoubleVector RnTable;
//   makeRnTable(RnTable, r, t);
//   return ip_theta_table(theta, r, t, RnTable);
//}


double GreensFunction3DRadInf::ip_theta_table(double theta, double r, double t, const DoubleVector& RnTable) const
{
   return ip_theta_free(theta, r, r0_, t, D_) + ip_corr_table(theta, r, RnTable);
}

static double p_free_max(double r, double r0, double t, double D)
{
   double Dt4 = 4.0 * D * t;
   double Dt4Pi = Dt4 * M_PI;
   double term1 = std::exp(-gsl_pow_2(r - r0) / Dt4);
   double term2 = 1.0 / std::sqrt(Dt4Pi * Dt4Pi * Dt4Pi);
   return term1 * term2;
}

// --------------------------------------------------------------------------------------------------------------------------------

DoubleVector GreensFunction3DRadInf::makeRnTable(double r, double t) const
{
   DoubleVector RnTable;

   {
      // First, estimate the size of p_corr, and if it's small enough
      // we don't need to calculate it in the first place.
      double pirr = p_irr(r, t, r0_, kf_, D_, sigma_);
      double ipfree_max = ip_theta_free(M_PI, r, r0_, t, D_) * 2 * M_PI * r * r;

      if (std::fabs((pirr - ipfree_max) / ipfree_max) < 1e-8)
         return RnTable;
   }

   double pfreemax = p_free_max(r, r0_, t, D_);
   double Rn_prev = 0.0;
   double RnFactor = 1.0 / (4.0 * M_PI * std::sqrt(r * r0_));
   double integrationTolerance = pfreemax / RnFactor * GfCfg.THETA_TOLERANCE;
   double truncationTolerance = pfreemax * GfCfg.THETA_TOLERANCE * 1e-1;
   gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(2000);

   uint MAX_ORDER = 70;    // differs from GfCfg.MAX_ORDER
   for (uint n = 0;; ++n)
   {
      double rn = Rn(n, r, t, workspace, integrationTolerance);
      RnTable.push_back(rn);

      // truncate when converged enough.
      double absRn = std::fabs(rn);
      if (absRn * RnFactor < truncationTolerance && absRn < Rn_prev) break;

      if (n >= MAX_ORDER)
      {
         _log.info() << "Rn didn't converge, t=" << std::setprecision(16) << t << ", r=" << r;
         break;
      }

      Rn_prev = absRn;
   }

   gsl_integration_workspace_free(workspace);
   return RnTable;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DRadInf::drawTheta(double rnd, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, r >= sigma_);
   THROW_UNLESS(std::invalid_argument, r0_ >= sigma_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   // t == 0 means no move.
   if (t == 0.0) return 0.0;

   DoubleVector RnTable = makeRnTable(r, t);
   double target = rnd * ip_theta_table(M_PI, r, t, RnTable);
   auto f = [r, t, &RnTable, target, this](double theta) { return ip_theta_table(theta, r, t, RnTable) - target; };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, 0.0, M_PI);
   return solver.findRoot(1e-15, GfCfg.THETA_TOLERANCE, "GreensFunction3DRadInf::drawTheta");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction3DRadInf::dump() const
{
   std::ostringstream ss;
   ss << std::setprecision(16) << "D=" << D_ << ", kf=" << kf_ << ", r0=" << r0_ << ", sigma=" << sigma_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
