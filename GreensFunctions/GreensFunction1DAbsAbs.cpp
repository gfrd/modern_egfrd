#include <cmath>
#include <algorithm>
#include <sstream>
#include "freeFunctions.hpp"
#include "helperFunctionsGf.hpp"
#include "GreensFunction1DAbsAbs.hpp"
#include "Logger.hpp"
#include <iomanip>

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction1DAbsAbs");

// --------------------------------------------------------------------------------------------------------------------------------

/* returns a guess for the number of terms needed for
   the greensfunction to converge at time t */
uint GreensFunction1DAbsAbs::guess_maxi(double t) const
{
   uint safety = 2;
   if (!std::isfinite(t)) return safety;

   double L = std::fabs(a_ - sigma_);
   double root0 = M_PI / L;
   double Dt = D_ * t;
   double thr = std::exp(-Dt * root0 * root0) * GfCfg.EPSILON * 1e-1;
   if (thr <= 0.0) return GfCfg.MAX_TERMS();

   double max_root = std::sqrt(root0 * root0 - log(thr) / Dt);
   uint maxi = std::max(safety + static_cast<uint>(max_root * L / M_PI), GfCfg.MIN_TERMS());
   return std::min(maxi, GfCfg.MAX_TERMS());
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction1DAbsAbs::p_survival(double t) const
{
   DoubleVector table;
   return p_survival_table(t, table);
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Calculates survival probability using a table.
Switchbox for which greensfunction to use. */
double GreensFunction1DAbsAbs::p_survival_table(double t, DoubleVector& psurvTable) const
{
   THROW_UNLESS(std::invalid_argument, t >= 0.0);
   double L = a_ - sigma_;

   if (std::fabs(r0_ - sigma_) < L*GfCfg.EPSILON || std::fabs(a_ - r0_) < L*GfCfg.EPSILON || L < 0.0)
      return 0.0;

   if (t == 0.0 || (D_ == 0.0 && v_ == 0.0)) return 1.0;     //particle can't escape.

  /* First check if we need full solution.
     Else we use approximation. */
   double distToa = a_ - r0_;
   double distTos = r0_ - sigma_;
   double maxDist = GfCfg.CUTOFF_H * (std::sqrt(2.0 * D_ * t) + std::fabs(v_*t));

   if (distToa > maxDist) //Absorbing boundary 'not in sight'.
   {
      if (distTos > maxDist)//And radiation boundary 'not in sight'.
         return 1.0; //No prob. outflux.
      return XS10(t, distTos, D_, v_); //Only absorbing BCn of s.
   }
   if (distTos > maxDist)
      return XS10(t, distToa, D_, -v_); //Only absorbing BCn of a.

   uint maxi = guess_maxi(t);
   if (maxi >= GfCfg.MAX_TERMS())
      _log.warn() << "drawT: maxi was cut to MAX_TERMS for t=" << std::setprecision(16) << t;

   if (psurvTable.size() < maxi)
      createPsurvTable(maxi, psurvTable);

   double p = funcSum_all(std::bind(&GreensFunction1DAbsAbs::p_survival_i, this, std::placeholders::_1, t, psurvTable), maxi);
   if (v_ == 0.0)
      p *= 2.0;
   else
   {
      double vexpo = -v_*v_*t / 4.0 / D_ - v_*r0_ / 2.0 / D_;
      p *= 2.0 * std::exp(vexpo);
   }

   return p;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Calculates the i'th term of the p_survival sum */
double GreensFunction1DAbsAbs::p_survival_i(uint i, double t, const DoubleVector& table) const
{
   double L = a_ - sigma_;
   return std::exp(-D_ * t * gsl_pow_2((i + 1) * M_PI / L)) * table[i];
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Calculates the part of the i'th term of p_surv not dependent on t, with drift */
double GreensFunction1DAbsAbs::p_survival_table_i_v(uint i) const
{
   double nPI = (i + 1) * M_PI;
   double L = a_ - sigma_;
   double r0s_L = (r0_ - sigma_) / L;
   double sigmav2D = sigma_*v_ / 2.0 / D_;
   double av2D = a_*v_ / 2.0 / D_;
   double Lv2D = L*v_ / 2.0 / D_;

   return (std::exp(sigmav2D) - std::cos(nPI)*std::exp(av2D)) * nPI / (Lv2D * Lv2D + nPI * nPI) * std::sin(nPI*r0s_L);
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Calculates the part of the i'th term of p_surv not dependent on t, without drift */
double GreensFunction1DAbsAbs::p_survival_table_i_nov(uint i) const
{
   double nPI = (i + 1) * M_PI;
   double L = a_ - sigma_;
   double r0s_L = (r0_ - sigma_) / L;
   return std::sin(nPI * r0s_L) * (1.0 - std::cos(nPI)) / nPI;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Fills table with terms in the p_survival sum which don't depend on t */
void GreensFunction1DAbsAbs::createPsurvTable(uint maxi, DoubleVector& table) const
{
   uint i = static_cast<uint>(table.size());
   table.reserve(maxi);
   if (v_ == 0.0)
   {
      while (i < maxi)
         table.emplace_back(p_survival_table_i_nov(i++));
   }
   else
   {
      while (i < maxi)
         table.emplace_back(p_survival_table_i_v(i++));
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Calculates the probability density of finding the particle at location r at time t.
double GreensFunction1DAbsAbs::prob_r(double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= (r - sigma_) && r <= a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   double L(a_ - sigma_);
   // if there was no time change or no diffusivity => no movement
   if (t == 0 || D_ == 0)
      return r == r0_ ? INFINITY : 0.0;        // the probability density function is a delta function
   if (std::fabs(r - sigma_) < L*GfCfg.EPSILON || std::fabs(a_ - r) < L*GfCfg.EPSILON || L < 0.0)
      return 0.0;

   // Set values that are constant in this calculation
   double expo = -D_*t / (L*L);
   double rs_L = (r - sigma_) / L;
   double r0s_L = (r0_ - sigma_) / L;
   double vexpo = -v_*v_*t / 4.0 / D_ + v_*(r - r0_) / 2.0 / D_;	// exponent of the drift-prefactor

   // Initialize summation
   double sum = 0, term = 0, prev_term;
   uint n = 0;
   do
   {
      if (n >= GfCfg.MAX_TERMS())
      {
         _log.warn("Too many terms for prob_r. N: %6u", n);
         break;
      }

      prev_term = term;
      double nPI = (n + 1) * M_PI;
      term = std::exp(nPI*nPI*expo) * std::sin(nPI * r0s_L) * std::sin(nPI * rs_L);
      sum += term;
      n++;
   } while (std::fabs(term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || std::fabs(prev_term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || n < GfCfg.MIN_TERMS());
   return 2.0 / L * std::exp(vexpo) * sum;
} */


// --------------------------------------------------------------------------------------------------------------------------------

/* Calculates the amount of flux leaving the left boundary at time t */
double GreensFunction1DAbsAbs::leaves(double t) const
{
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   double L = a_ - sigma_;
   if (std::fabs(r0_ - sigma_) < L*GfCfg.EPSILON || std::fabs(a_ - r0_) < L*GfCfg.EPSILON || L < 0.0)
   {
      // The flux of a zero domain is INFINITY. Also if the particle 
      // started on the left boundary (leaking out immediately).
      return INFINITY;
   }
   if (t < GfCfg.EPSILON*tscale_)
      return 0.0; // if t=0.0 the flux must be zero

   double sum = 0, term = 0, prev_term;
   double D_L_sq = D_ / (L*L);
   double expo = -D_L_sq*t;
   double r0s_L = (r0_ - sigma_) / L;
   double vexpo = -v_*v_*t / 4.0 / D_ - v_*(r0_ - sigma_) / 2.0 / D_;
   uint n = 0;
   do
   {
      if (n >= GfCfg.MAX_TERMS())
      {
         _log.warn() << "Too many terms for leaves. N=" << n;
         break;
      }

      double nPI = (n + 1) * M_PI;
      prev_term = term;
      term = nPI * std::exp(nPI * nPI * expo) * std::sin(nPI * r0s_L);
      sum += term;
      n++;
   } while (std::fabs(term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || std::fabs(prev_term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || n < GfCfg.MIN_TERMS());
   return 2.0 * D_L_sq * std::exp(vexpo) * sum;
}

// --------------------------------------------------------------------------------------------------------------------------------

// Calculates the amount of flux leaving the right boundary at time t
double GreensFunction1DAbsAbs::leavea(double t) const
{
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   double L = a_ - sigma_;

   if (std::fabs(r0_ - sigma_) < L*GfCfg.EPSILON || std::fabs(a_ - r0_) < L*GfCfg.EPSILON || L < 0.0)
   {
      // The flux of a zero domain is INFINITY. Also if the particle 
      // started on the right boundary (leaking out immediately).
      return INFINITY;
   }
   if (t < GfCfg.EPSILON*tscale_)
      return 0.0;     // if t=0.0 the flux must be zero

   double sum = 0, term = 0, prev_term;
   double D_L_sq = D_ / (L*L);
   double expo = -D_L_sq*t;		// exponent -D n^2 PI^2 t / l^2
   double r0s_L = (r0_ - sigma_) / L;
   double vexpo = -v_*v_*t / 4.0 / D_ + v_*(a_ - r0_) / 2.0 / D_;
   uint n = 0;
   do
   {
      if (n >= GfCfg.MAX_TERMS())
      {
         _log.warn() << "Too many terms for leavea. N=" << n;
         break;
      }

      double nPI = (n + 1) * M_PI;
      prev_term = term;
      term = nPI * std::exp(nPI * nPI * expo) * std::cos(nPI) * std::sin(nPI * r0s_L);
      sum += term;
      n++;
   } while (std::fabs(term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || std::fabs(prev_term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || n < GfCfg.MIN_TERMS());

   return -2.0 * D_L_sq * std::exp(vexpo) * sum;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* This draws an eventtype of time t based on the flux through the left (z=sigma)
   and right (z=a) boundary. Although not completely accurate, it returns an
   IV_ESCAPE for an escape through the right boundary and a IV_REACTION for an
   escape through the left boundary. */
GreensFunction::EventKind GreensFunction1DAbsAbs::drawEventType(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, t > 0.0);

   // if t=0 nothing has happened => no event

   double L = a_ - sigma_;

   // For particles at the boundaries
   if (std::fabs(a_ - r0_) < GfCfg.EPSILON*L)
      return EventKind::IV_ESCAPE;        // if the particle started on the right boundary
   if (std::fabs(r0_ - sigma_) < GfCfg.EPSILON*L)
      return EventKind::IV_REACTION;         // if the particle started on the left boundary

   double leaves_s = leaves(t);
   double leaves_a = leavea(t);
   double flux_total = leaves_s + leaves_a;
   double fluxratio = leaves_s / flux_total;
   return rnd > fluxratio ? EventKind::IV_ESCAPE : EventKind::IV_REACTION;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction1DAbsAbs::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd < 1.0);
   double L = a_ - sigma_;

   if (D_ == 0.0) return INFINITY;
   if (L < 0.0 || std::fabs(a_ - r0_) < GfCfg.EPSILON*L || std::fabs(r0_ - sigma_) > (1.0 - GfCfg.EPSILON)*L) return 0.0;

   DoubleVector psurvTable;
   auto f = [rnd, &psurvTable, this](double t) { return rnd - p_survival_table(t, psurvTable); };
   gsl_lambda<decltype(f)> F(f);

   /* Find a good interval to determine the first passage time in */
   double dist = std::min(r0_ - sigma_, a_ - r0_);
   double t_guess = 0;
   if (v_ == 0.0)
   {
      t_guess = dist * dist / (2.0 * D_);
   }
   else
   {
      // When drifting towards the closest boundary...
      if ((r0_ - sigma_ >= L / 2.0 && v_ > 0.0) || (r0_ - sigma_ <= L / 2.0 && v_ < 0.0))
         t_guess = std::sqrt(D_*D_ / (v_*v_*v_*v_) + dist*dist / (v_*v_)) - D_ / (v_*v_);

      // When drifting away from the closest boundary...
      if ((r0_ - sigma_  < L / 2.0 && v_ > 0.0) || (r0_ - sigma_ > L / 2.0 && v_ < 0.0))
         t_guess = D_ / (v_*v_) - std::sqrt(D_*D_ / (v_*v_*v_*v_) - dist*dist / (v_*v_));
   }
   t_guess *= .1;

   double value = GSL_FN_EVAL(&F, t_guess);
   double low = t_guess;
   double high = t_guess;
   if (value < 0.0)
   {
      // scale the interval around the guess such that the function 
      // straddles if the guess was too low
      do
      {
         if (std::fabs(high) >= t_guess * 1e6)
         {
            std::stringstream msg;
            msg << type_name() << ": couldn't adjust high. F(" << high << ")=" << value << "; " << dump();
            throw std::runtime_error(msg.str());
         }
         // keep increasing the upper boundary until the 
         // function straddles
         high *= 10.0;
         value = GSL_FN_EVAL(&F, high);
      } while (value <= 0.0);
   }
   else
   {
      // if the guess was too high initialize with 2 so the test 
      // below survives the first iteration
      double value_prev = 2.0;
      do
      {
         if (std::fabs(low) <= t_guess * 1.0e-6 || std::fabs(value - value_prev) < GfCfg.EPSILON*tscale_)
         {
            _log.warn() << type_name() << ": couldn't adjust low. F(" << std::setprecision(16) << low << ")=" << value << ", " << dump();
            return low;
         }

         value_prev = value;
         // keep decreasing the lower boundary until the 
         // function straddles
         low *= 0.1;
         // get the accompanying value
         value = GSL_FN_EVAL(&F, low);

      } while (value >= 0.0);
   }

   // find the intersection on the y-axis between the random number and the function
   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(GfCfg.EPSILON*tscale_, GfCfg.EPSILON, "GreensFunction1DAbsAbs::drawTime");
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction1DAbsAbs::p_int_r_table(double r, double t, DoubleVector& table) const
{
   double distToa = a_ - r0_;
   double distTos = r0_ - sigma_;
   double maxDist = GfCfg.CUTOFF_H * (std::sqrt(2.0 * D_ * t) + std::fabs(v_ * t));

   if (distToa > maxDist) //Absorbing boundary a 'not in sight'.
   {
      if (distTos > maxDist) //Absorbing boundary sigma 'not in sight'.
         return XI00(r, t, r0_, D_, v_); //free particle.
     //Only absorbing BCn at sigma.
      return XI10(r - sigma_, t, distTos, D_, v_);
   }
   if (distTos > maxDist)
      //Only absorbing BCn at a.
      return XI10(a_ - r, t, distToa, D_, -v_);

   double vexpo = -v_*v_*t / 4.0 / D_ - v_*r0_ / 2.0 / D_;
   double prefac = 2.0 * std::exp(vexpo);

   uint maxi = guess_maxi(t);
   if (table.size() < maxi)
      create_p_int_r_Table(t, maxi, table);

   if (maxi >= GfCfg.MAX_TERMS())
      _log.warn() << "p_int_r_table: maxi was cut to MAX_TERMS for t=" << std::setprecision(16) << t;

   double p = funcSum(std::bind(&GreensFunction1DAbsAbs::p_int_r_i, this, std::placeholders::_1, r, t, table), GfCfg.MAX_TERMS());
   return prefac * p;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction1DAbsAbs::p_int_r_i(uint i, double r, double t, DoubleVector& table) const
{
   double L = a_ - sigma_;
   double v2D = v_ / (2 * D_);
   double n_L = (i + 1.0) * M_PI / L;

   double term;
   if (v2D == 0.0)
      term = 1.0 - std::cos(n_L*(r - sigma_));
   else
      term = std::exp(v2D*sigma_) + std::exp(v2D*r) * (v2D / n_L*std::sin(n_L*(r - sigma_)) - std::cos(n_L*(r - sigma_)));

   return term * get_p_int_r_Table_i(i, t, table);
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Fills table for p_int_r of factors independent of r. */
void GreensFunction1DAbsAbs::create_p_int_r_Table(double t, uint maxi, DoubleVector& table) const
{
   uint n = static_cast<uint>(table.size());
   double L = a_ - sigma_;
   double expo = -D_*t / (L*L);
   double r0s_L = (r0_ - sigma_) / L;
   double Lv2D = L*v_ / 2.0 / D_;
   table.reserve(maxi);

   double nPI, term;
   while (n < maxi)
   {
      nPI = (n + 1)*M_PI;
      if (v_ == 0.0)
         term = std::exp(nPI*nPI*expo) * std::sin(nPI*r0s_L) / nPI;
      else
         term = std::exp(nPI*nPI*expo) * std::sin(nPI*r0s_L) * nPI / (nPI*nPI + Lv2D*Lv2D);
      table.emplace_back(term);
      n++;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Draws the position of the particle at a given time from p(r,t), assuming
   that the particle is still in the domain */
double GreensFunction1DAbsAbs::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   double L = a_ - sigma_;

   // the trivial case: if there was no movement or the domain was zero
   if ((D_ == 0.0 && v_ == 0.0) || L < 0.0 || t == 0.0) return r0_;

   // if the initial condition is at the boundary, raise an error
   // The particle can only be at the boundary in the ABOVE cases
   THROW_UNLESS(std::invalid_argument, (r0_ - sigma_) >= L*GfCfg.EPSILON && (r0_ - sigma_) <= L*(1.0 - GfCfg.EPSILON));

   DoubleVector pintTable;
   double rndpsurf = rnd * p_survival(t);
   auto f = [t, &pintTable, rndpsurf, this](double r) { return p_int_r_table(r, t, pintTable) - rndpsurf; };
   gsl_lambda<decltype(f)> F(f);

   // find the intersection on the y-axis between the random number and the function
   root_fsolver_wrapper solver;
   solver.set(&F, sigma_, a_);
   return solver.findRoot(L*GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction1DAbsAbs::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction1DAbsAbs::dump() const
{
   std::ostringstream ss;
   ss << "D=" << D_ << ", sigma=" << sigma_ << ", a=" << a_ << ", v=" << v_ << ", r0=" << r0_ << std::endl;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
