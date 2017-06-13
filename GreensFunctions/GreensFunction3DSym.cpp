#include <sstream>
#include <cmath>
#include "helperFunctionsGf.hpp"
#include "GreensFunction3DSym.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DSym::p_r(double r, double t) const
{
   double Dt4 = 4.0 * D_ * t;
   double Dt4Pi = Dt4 * M_PI;
   double term1 = 1.0 / std::sqrt(gsl_pow_3(Dt4Pi));
   double term2 = std::exp(-r * r / Dt4);
   double jacobian = 4.0 * r * r * M_PI;
   return jacobian * term1 * term2;
}

double GreensFunction3DSym::ip_r(double r, double t) const
{
   double Dt = D_ * t;
   double sqrtDt_r = 1.0 / std::sqrt(D_ * t);
   double sqrtPi_r = 1.0 / std::sqrt(M_PI);
   double term1 = std::exp(-r * r / (4.0 * Dt)) * r * sqrtDt_r * sqrtPi_r;
   double term2 = std::erf(r * 0.5 * sqrtDt_r);
   return term2 - term1;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DSym::drawR(double rnd, double t) const
{
   // input parameter range checks.
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0 || D_ == 0.0) return 0.0;

   auto f = [rnd, t, this](double r) { return ip_r(r, t) - rnd; };
   gsl_lambda<decltype(f)> F(f);

   double max_r = 4.0 * std::sqrt(6.0 * D_ * t);
   while (GSL_FN_EVAL(&F, max_r) < 0.0)
      max_r *= 10;

   root_fsolver_wrapper solver;
   solver.set(&F, 0.0, max_r);
   return solver.findRoot(1e-15, GfCfg.TOLERANCE, "GreensFunction3DSym::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction3DSym::dump() const
{
   std::ostringstream ss;
   ss << "D=" << D_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
