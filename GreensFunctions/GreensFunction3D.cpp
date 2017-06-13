#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "freeFunctions.hpp"
#include "helperFunctionsGf.hpp"
#include "GreensFunction3D.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3D::p_r(double r, double t) const
{
   double Dt = D_ * t;
   double Dt4 = 4.0 * Dt;
   double rr04 = 4.0 * r * r0_;

   double mrr0sq_over_4Dt = -gsl_pow_2(r + r0_) / Dt4;
   double num1 = std::expm1(mrr0sq_over_4Dt);
   double num2 = std::expm1(mrr0sq_over_4Dt + rr04 / Dt4);
   double den = rr04 * std::sqrt(M_PI * M_PI * M_PI * Dt);
   double jacobian = 2.0 * r * r * M_PI;
   return jacobian * (-num1 + num2) / den;
}

double GreensFunction3D::ip_r(double r, double t) const
{
   double Dt4 = 4.0 * D_ * t;
   double Dt4r = 1.0 / Dt4;
   double sqrtDt4 = std::sqrt(Dt4);
   double sqrtDt4r = 1.0 / sqrtDt4;
   double num1a = std::exp(-gsl_pow_2(r - r0_) * Dt4r);
   double num1b = std::exp(-gsl_pow_2(r + r0_) * Dt4r);
   double den1 = r0_ * std::sqrt(M_PI);
   double term1 = sqrtDt4 * (-num1a + num1b) / den1;
   double term2 = std::erf((r - r0_) * sqrtDt4r);
   double term3 = std::erf((r + r0_) * sqrtDt4r);
   return 0.5 * (term1 + term2 + term3);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3D::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0) return r0_;
  
   auto f = [rnd, t, this](double r) { return ip_r(r, t) - rnd; };
   gsl_lambda<decltype(f)> F(f);

   const uint H = 7;
   double r_range = H * std::sqrt(6.0 * D_ * t);
   double low_r = std::max(r0_ - r_range, 0.0);
   double max_r = r0_ + r_range;

   if (GSL_FN_EVAL(&F, low_r) >= 0.0) return low_r;
   if (GSL_FN_EVAL(&F, max_r) <= 0.0) return max_r;

   root_fsolver_wrapper solver;
   solver.set(&F, low_r, max_r);
   return solver.findRoot(1e-15, GfCfg.TOLERANCE, "GreensFunction3D::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3D::drawTheta(double rnd, double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, r >= 0.0);
   THROW_UNLESS(std::invalid_argument, r0_ >= 0.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0) return 0.0;

   double target = rnd * ip_theta_free(M_PI, r, r0_, t, D_);
   auto f = [r, t, target, this](double theta) { return ip_theta_free(theta, r, r0_, t, D_) - target; };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, 0.0, M_PI + std::numeric_limits<double>::epsilon());
   return solver.findRoot(1e-15, GfCfg.THETA_TOLERANCE, "GreensFunction3D::drawTheta");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction3D::dump() const
{
   std::ostringstream ss;
   ss << std::setprecision(16) << "D=" << D_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
