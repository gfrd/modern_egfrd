#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "helperFunctionsGf.hpp"
#include "GreensFunction3DAbsSym.hpp"
#include "Logger.hpp"
#include <gsl/gsl_sf.h>

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction3DAbsSym");

// --------------------------------------------------------------------------------------------------------------------------------

// Efficiently calculate EllipticTheta[4,0,q] for q < 1.0.
double GreensFunction3DAbsSym::ellipticTheta4Zero(double q) const
{
   THROW_UNLESS(std::invalid_argument, std::fabs(q) <= 1.0);

   double value = 1.0;
   double q_n = q;
   double q_2n = 1.0;

   uint i = 1000; // max iterations
   while (--i)
   {
      double term2 = 1.0 - q_2n * q;  // q^(2n-1) = (q^(n-1))^2 * q

      q_2n = q_n * q_n;

      double term1 = 1.0 - q_2n;      // q^2n

      double term = term1 * term2 * term2;
      double value_prev = value;
      value *= term;

      if (std::fabs(value - value_prev) < 1e-18)
         return value;

      q_n *= q;  // q_(++n)
   }

   _log.warn() << "ellipticTheta4Zero: didn't converge";
   return value;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbsSym::p_survival(double t) const //eq. 11 implementation notes.
{
   double q = -D_ * M_PI * M_PI * t / (a_ * a_);
   return 1.0 - ellipticTheta4Zero(std::exp(q)); // = S(t)
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbsSym::p_int_r_free(double r, double t) const //eq.4 implementation notes.
{
   double sqrtDt = std::sqrt(D_ * t);
   return std::erf(r / (sqrtDt + sqrtDt)) - r * std::exp(-r * r / (4.0 * D_ * t)) / (std::sqrt(M_PI) * sqrtDt);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbsSym::p_int_r(double r, double t) const //eqS. 13 implementation notes.
{
   double p_free = p_int_r_free(r, t);

   // p_int_r is always smaller than p_free.
   if (std::fabs(p_free) < GfCfg.CUTOFF) return 0.0;

   double DtPIsq_asq = D_ * t * M_PI * M_PI / (a_ * a_);
   double maxn = (a_ / M_PI) * std::sqrt(std::log(std::exp(DtPIsq_asq) / GfCfg.CUTOFF) / (D_ * t));

   const int N_MAX = 10000;
   int N = std::min(static_cast<int>(std::ceil(maxn) + 1), N_MAX);
   THROW_UNLESS_MSG(std::invalid_argument, N < N_MAX, "p_int_r: didn't converge");

   double angle = M_PI * r / a_;
   double value = 0.0;
   for (int n = 1; n <= N; ++n)
      value += std::exp(-n * n * DtPIsq_asq) * (a_ * std::sin(n * angle) - n * M_PI * r * std::cos(n * angle)) / n;

   return value * 2.0 / (a_ * M_PI);
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbsSym::p_r_fourier(double r, double t) const
{
   double value = 0.0;
   double asq = a_ * a_;
   double PIsq = M_PI * M_PI;

   const uint maxIter = 100;
   for (uint n = 1;; ++n)
   {
      double term1 = std::exp(-(PIsq * r * r + asq * n*n) / (4.0 * D_ * PIsq * t));
      double term2 = M_PI * r * std::exp(gsl_sf_lncosh(a_ * r * n / (2.0 * D_ * M_PI * t)));
      double term3 = a_ * n * std::exp(gsl_sf_lnsinh(a_ * r * n / (2.0 * D_ * M_PI * t)));

      double term = term1 * r * (term2 - term3);
      value += term;

      if (std::fabs(value) * 1e-8 > std::fabs(term))
         break;

      if (n > maxIter)
      {
         _log.warn() << "p_r_fourier didn't converge; n=" << n << ", value=" << value;
         break;
      }
   }

   return value / (M_SQRT2 * PIsq * std::pow(D_ * t, 1.5));
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbsSym::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);

   if (D_ == 0.0 || a_ == INFINITY) return INFINITY;
   if (a_ == 0.0) return 0.0;

   auto f = [rnd, this](double t) { return rnd - p_survival(t); };
   gsl_lambda<decltype(f)> F(f);

   double t_guess = a_ * a_ / (6.0 * D_);
   double low = t_guess;
   double high = t_guess;
   double value = GSL_FN_EVAL(&F, t_guess);

   if (value < 0.0)
   {
      high *= 10;
      for (;;)
      {
         double high_value = GSL_FN_EVAL(&F, high);
         if (high_value >= 0.0) break;

         if (std::fabs(high) >= t_guess * 1e6)
         {
            std::stringstream msg;
            msg << type_name() << ": couldn't adjust high. F(" << high << ")=" << high_value << "; " << dump();
            throw std::runtime_error(msg.str());
         }
         high *= 10;
      }
   }
   else
   {
      double low_value_prev = value;
      low *= .1;
      for (;;)
      {
         double low_value = GSL_FN_EVAL(&F, low);
         if (low_value <= 0.0) break;

         if (std::fabs(low) <= t_guess * 1e-6 || std::fabs(low_value - low_value_prev) < GfCfg.CUTOFF)
         {
            _log.warn() << "drawTime: couldn't adjust low. F(" << low << ")=" << GSL_FN_EVAL(&F, low) << ", " << dump();
            return low;
         }

         low_value_prev = low_value;
         low *= .1;
      }
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-18, 1e-12, "GreensFunction3DAbsSym::drawTime");
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction3DAbsSym::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (a_ == 0.0 || t == 0.0 || D_ == 0.0) return 0.0;

   root_fsolver_wrapper solver;

   double thresholdDistance = GfCfg.CUTOFF_H * std::sqrt(6.0 * D_ * t);
   if (a_ <= thresholdDistance)
   {
      //psurv = p_survival(t);  // this causes a problem when p_survival is very small.
      double psurv = p_int_r(a_, t);
      if (psurv == 0.0) return a_;
      ASSERT(psurv >= 0.0);

      const double target = psurv * rnd;
      auto f = [t, target, this](double r) { return p_int_r(r, t) - target; };
      auto F = gsl_lambda<decltype(f)>(f);
      solver.set(&F, 0, a_);
      return solver.findRoot(1e-18, 1e-12, "GreensFunction3DAbsSym::drawR");
   }

   // p_int_r < p_int_r_free
   if (p_int_r_free(a_, t) < rnd)
   {
      _log.info() << "p_int_r_free(a_, t) < rnd, returning a";
      return a_;
   }

   auto f = [t, rnd, this](double r) { return p_int_r_free(r, t) - rnd; };
   auto F = gsl_lambda<decltype(f)>(f);
   solver.set(&F, 0, a_);
   return solver.findRoot(1e-18, 1e-12, "GreensFunction3DAbsSym::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction3DAbsSym::dump() const
{
   std::ostringstream ss;
   ss << std::scientific << std::setprecision(16) << "D=" << D_ << ", a=" << a_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
