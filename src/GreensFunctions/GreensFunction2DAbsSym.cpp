#include <sstream>
#include <gsl/gsl_sf_bessel.h>
#include "helperFunctionsGf.hpp"
#include "GreensFunction2DAbsSym.hpp"
#include "Logger.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction2DAbsSym");

// --------------------------------------------------------------------------------------------------------------------------------

// an alternative form, which is not very convergent.
double GreensFunction2DAbsSym::p_survival(double t) const
{
   double Dt(-D_ * t);
   int N(100);	// number of terms to use
   double sum(0.);
   double threshold(GfCfg.CUTOFF);
   for (int n(1); n <= N; ++n)
   {
      double aAn = gsl_sf_bessel_zero_J0(n);		// gsl roots of J0(aAn)
      double An = aAn / a_;
      double J1_aAn = gsl_sf_bessel_J1(aAn);
      double term = (exp(An*An*Dt)) / (An*J1_aAn);
      sum += term;

      if (fabs(term / sum) < threshold) break;
   }
   return (2.0 / a_) * sum;
}

double GreensFunction2DAbsSym::p_int_r_free(double r, double t) const
{
   double Dt(D_ * t);
   double sqrtDt(sqrt(Dt));
   double sqrtPI(sqrt(M_PI));
   return erf(r / (sqrtDt + sqrtDt)) - r * exp(-r * r / (4.0 * Dt)) / (sqrtPI * sqrtDt);
}

double GreensFunction2DAbsSym::p_int_r(double r, double t) const
{
   double Dt(-D_ * t);
   double term;
   double sum(0.0);
   int n(1);

   //    double maxn( ( a / M_PI ) * sqrt( log( exp( DtPIsq_asq ) / GfCfg.CUTOFF ) / 
   //                                          ( D * t ) ) );

   int N_MAX(10000);
   double threshold(GfCfg.CUTOFF);

   do
   {
      double aAn = gsl_sf_bessel_zero_J0(n);         // gsl roots of J0(aAn)
      double An = aAn / a_;
      double rAn = r*An;
      double J1_aAn = gsl_sf_bessel_J1(aAn);
      double J1_rAn = gsl_sf_bessel_J1(rAn);
      term = (exp(An*An*Dt) * r * J1_rAn) / (An*J1_aAn*J1_aAn);
      sum += term;
      n++;
   } while (fabs(term / sum) > threshold && n <= N_MAX);

   return (2.0 / (a_*a_)) * sum;
}

double GreensFunction2DAbsSym::p_survival_F(double t, const p_survival_params* params)
{
   const GreensFunction2DAbsSym* const gf(params->gf);
   return 1 - gf->p_survival(t) - params->rnd;
}

double GreensFunction2DAbsSym::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   if (D_ == 0.0 || a_ == INFINITY) return INFINITY;
   if (a_ == 0.0) return 0.0;

   p_survival_params params = { this, rnd };
   gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&p_survival_F), &params };


   // Find a good interval to determine the first passage time in
   double t_guess(a_ * a_ / (4. * D_));   // construct a guess: msd = sqrt (2*d*D*t)
   double value(GSL_FN_EVAL(&F, t_guess));
   double low(t_guess);
   double high(t_guess);

   // scale the interval around the guess such that the function straddles
   if (value < 0.0)               // if the guess was too low
   {
      do
      {
         high *= 10;     // keep increasing the upper boundary until the function straddles
         value = GSL_FN_EVAL(&F, high);

         if (fabs(high) >= t_guess * 1e6)
         {
            _log.warn() << "Couldn't adjust high. F(" << high << ")=" << value;
               throw std::exception();
         }
      } while (value <= 0.0);
   }
   else                            // if the guess was too high
   {
      double value_prev(value);
      do
      {
         low *= .1;      // keep decreasing the lower boundary until the function straddles
         value = GSL_FN_EVAL(&F, low);     // get the accompanying value

         if (fabs(low) <= t_guess * 1e-6 || fabs(value - value_prev) < GfCfg.CUTOFF)
         {
            _log.warn() << "Couldn't adjust low. F(" << low << ")=" << value;
            return low;
         }
         value_prev = value;
      } while (value >= 0.0);
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-18, 1e-12, "GreensFunction2DAbsSym::drawTime");
}

double GreensFunction2DAbsSym::p_r_free_F(double r, const p_r_params* params)
{
   const GreensFunction2DAbsSym* const gf(params->gf);
   return gf->p_int_r_free(r, params->t) - params->target;
}

double GreensFunction2DAbsSym::p_r_F(double r, const p_r_params* params)
{
   const GreensFunction2DAbsSym* const gf(params->gf);
   return gf->p_int_r(r, params->t) - params->target;
}

double GreensFunction2DAbsSym::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (a_ == 0.0 || t == 0.0 || D_ == 0.0) return 0.0;

   //double thresholdDistance( GfCfg.CUTOFF_H * sqrt( 4.0 * D * t ) );

   //  if( a_ <= thresholdDistance )	// if the domain is not so big, the boundaries are felt
   //  {
   double psurv = p_survival(t);
   //psurv = p_int_r( a_, t );
   //printf("dr %g %g\n",psurv, p_survival( t ));
   //assert( fabs(psurv - p_int_r( a_, t )) < psurv * 1e-8 );

   ASSERT(psurv > 0.0);

   gsl_function F;
   F.function = reinterpret_cast<double(*)(double, void*)>(&p_r_F);
   /*  }
       else				// if the domain is very big, just use the free solution
       {
       // p_int_r < p_int_r_free
       if( p_int_r_free( a_, t ) < rnd )	// if the particle is outside the domain?
       {
       std::cerr << "p_int_r_free( a_, t ) < rnd, returning a_."
       << std::endl;
       return a_;
       }

       psurv = 1.0;
       F.function = reinterpret_cast<double(*)(double, void*)>( &p_r_free_F );
       }

       */
   double target(psurv * rnd);
   p_r_params params = { this, t, target };
   F.params = &params;

   double low(0.0);
   double high(a_);
   //double high( std::min( thresholdDistance, a_ ) );

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-18, 1e-12, "GreensFunction2DAbsSym::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction2DAbsSym::dump() const
{
   std::ostringstream ss;
   ss << "D=" << D_ << ", a=" << a_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
