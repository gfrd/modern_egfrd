#include <cmath>
#include <algorithm>
#include <sstream>
#include "freeFunctions.hpp"
#include "helperFunctionsGf.hpp"
#include "GreensFunction1DRadAbs.hpp"
#include "Logger.hpp"
#include <iomanip>

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction1DRadAbs");

// --------------------------------------------------------------------------------------------------------------------------------

// This is the appropriate definition of the function defining
// the roots of our Green's functions in GSL.
// Later needed by the rootfinder.
//
// It expects a reaction rate h=k/D already divided by D_.
double GreensFunction1DRadAbs::tan_f(double x, void *p)
{
   struct tan_f_params *params = static_cast<struct tan_f_params *>(p);
   const double h_a(params->h * params->a);

   // h = k/D
   if (fabs(h_a) < 1) return 1 / tan(x) + (h_a) / x;
   return tan(x) + x / (h_a);
}

/* return the n + 1'th root */
double GreensFunction1DRadAbs::get_root(uint n) const
{
   if (n >= rootList_size())
      calculate_n_roots(n + 1);
   return rootList_[n];
}

/* Fills the rootList_ with all the roots of tan(x*a)=-x/h up to n */
void GreensFunction1DRadAbs::calculate_n_roots(uint n) const
{
   uint i(rootList_size());
   if (n <= i) return;

   // lift the const veil and fill the table
   auto ucthis = const_cast<GreensFunction1DRadAbs*>(this);
   ucthis->fill_table_to_n(i, n);
}

void GreensFunction1DRadAbs::fill_table_to_n(uint i, uint n)
{
   const double L(a_ - sigma_);
   const double h((k_ + v_ / 2.0) / D_);
   // the drift v also comes into this constant, h=(k+v/2)/D
   double upper, lower, root_i;

   // reserve space
   rootList_.reserve(n);

   //No drift, and k = 0, use reflective solution.
   if (k_ < GfCfg.EPSILON && fabs(v_) < GfCfg.EPSILON)
   {
      while (i++ < n)
         rootList_.emplace_back(M_PI * (i - 0.5) / L);
      return;
   }

   gsl_function F;
   struct tan_f_params params = { L, h };

   F.function = &GreensFunction1DRadAbs::tan_f;
   F.params = &params;

   root_fsolver_wrapper solver;

   /* Find all the roots up to the nth */
   if (h*L < 1)
   {
      lower = i*M_PI + 1E-10;
      upper = (i + 1) * M_PI - 1E-10;
   }
   else
   {
      lower = i * M_PI + M_PI_2 + 1E-10;
      upper = (i + 1) * M_PI + M_PI_2 - 1E-10;
   }

   while (i++ < n)
   {
      solver.set(&F, lower, upper);
      root_i = solver.findRoot(GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction1DRadAbs::root_tan");
      //root_i = findRoot(F, solver, lower, upper, 1.0*GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction1DRadAbs::root_tan");
      rootList_.emplace_back(root_i / L);
      lower += M_PI;
      upper += M_PI;
   }
}

/* returns a guess for the number of terms needed for
   the greensfunction to converge at time t */
uint GreensFunction1DRadAbs::guess_maxi(double t) const
{
   const uint safety(2);

   if (!std::isfinite(t)) return safety;

   const double L(fabs(a_ - sigma_));
   const double root0(get_root(0));
   const double Dt(D_ * t);

   const double thr(exp(-Dt * root0 * root0) * GfCfg.EPSILON * 1e-1);

   if (thr <= 0.0) return GfCfg.MAX_TERMS();
   const double max_root(sqrt(root0 * root0 - log(thr) / Dt));
   const uint maxi(std::max(safety + static_cast<uint> (max_root * L / M_PI), GfCfg.MIN_TERMS()));
   return std::min(maxi, GfCfg.MAX_TERMS());
}

// This is the non-exponential factor in the Green's function sum, not
// including the factor containing the explicit r-dependency (The latter
// is given by the Bn's, see below).
//
// r0 is here still in the interval from 0 to a (and supposed to be the
// starting point of the particle at t0).
//
// The root a_n also must be the specific one for that interval, thus
// the one rescaled by a (see comments in function a_n(n) ).
//
// The factor calculated here is identical for the cases w. or w/o drift
// only h changes.
double GreensFunction1DRadAbs::An(double root_n) const
{
   const double h((k_ + v_ / 2.0) / D_);

   const double L(a_ - sigma_);
   const double rootn_r0_s = root_n*(r0_ - sigma_);

   return (root_n*cos(rootn_r0_s) + h*sin(rootn_r0_s)) / (h + (root_n*root_n + h*h)*L);
}

// This factor appears in the survival prob.
double GreensFunction1DRadAbs::Bn(double root_n) const
{
   const double h((k_ + v_ / 2.0) / D_);
   const double L(a_ - sigma_);

   const double rootnL(root_n*L);
   const double rootn2(root_n*root_n);
   const double h2(h*h);
   const double v2D(v_ / 2.0 / D_);

   if (v_ == 0.0) return (h2 - (rootn2 + h2)*cos(rootnL)) / (h*root_n);
   return (exp(v2D*sigma_)*h*k_ / D_ - exp(v2D*a_)*(rootn2 + h2)*cos(rootnL)) / (h / root_n*(rootn2 + v2D*v2D));
}

// This is the exponential factor in the Green's function sum, also
// appearing in the survival prob. and prop. function.
//
// Also here the root is the one referring to the interval of length L.
double GreensFunction1DRadAbs::Cn(double root_n, double t) const
{
   return exp(-D_*root_n*root_n*t);
}

double GreensFunction1DRadAbs::p_survival(double t) const
{
   DoubleVector table;
   return p_survival_table(t, table);
}

/* Calculates survival probability using a table.
   Switchbox for which greensfunction to use. */
double GreensFunction1DRadAbs::p_survival_table(double t, DoubleVector& psurvTable) const
{
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   const double L(a_ - sigma_);

   if (fabs(a_ - r0_) < L*GfCfg.EPSILON || L < 0.0)
      return 0.0; // The survival probability of a zero domain is zero

   if (t == 0.0 || (D_ == 0.0 && v_ == 0.0))
      return 1.0;     //particle can't escape.

  /* First check if we need full solution.
     Else we use approximation. */
   const double distToa(a_ - r0_);
   const double distTos(r0_ - sigma_);
   const double maxDist(GfCfg.CUTOFF_H * (sqrt(2.0 * D_ * t) + fabs(v_ * t)));

   if (distToa > maxDist) //Absorbing boundary 'not in sight'.
   {
      if (distTos > maxDist) //Radiation boundary 'not in sight'.
         return 1.0; //No prob. outflux.
      return XS30(t, distTos, k_, D_, v_); //Only radiation BCn.
   }

   if (distTos > maxDist)
      return XS10(t, distToa, D_, -v_); //Only absorbing BCn.

   const uint maxi(guess_maxi(t));

   if (maxi >= GfCfg.MAX_TERMS())
      _log.warn() << "drawT: maxi was cut to MAX_TERMS for t=" << std::scientific << std::setprecision(16) << t;

   if (psurvTable.size() < maxi)
   {
      calculate_n_roots(maxi);
      createPsurvTable(psurvTable);
   }

   double p = funcSum_all(std::bind(&GreensFunction1DRadAbs::p_survival_i, this, std::placeholders::_1, t, psurvTable), maxi);
   if (v_ == 0.0)
   {
      p *= 2.0;
   }
   else
   {
      const double vexpo(-v_*v_*t / 4.0 / D_ - v_*r0_ / 2.0 / D_);
      p *= 2.0 * exp(vexpo);
   }
   return p;
}

/* Calculates the i'th term of the p_survival sum */
double GreensFunction1DRadAbs::p_survival_i(uint i, double t, const DoubleVector& table) const
{
   return exp(-D_ * t * gsl_pow_2(get_root(i))) * table[i];
}

/* Calculates the part of the i'th term of p_surv not dependent on t, with drift */
double GreensFunction1DRadAbs::p_survival_table_i_v(uint i) const
{
   const double L(a_ - sigma_);
   const double h((k_ + v_ / 2.0) / D_);
   const double v2D(v_ / 2.0 / D_);
   const double exp_av2D(exp(a_*v2D));
   const double exp_sigmav2D(exp(sigma_*v2D));
   const double root_n(get_root(i));
   const double root_n2 = root_n * root_n;
   const double root_n_r0_s = root_n * (r0_ - sigma_);
   const double root_n_L = root_n * L;
   const double h_root_n = h / root_n;
   return (h * sin(root_n_r0_s) + root_n * cos(root_n_r0_s)) / (L * (root_n2 + h * h) + h) * (exp_sigmav2D * h * k_ / D_ - exp_av2D * (root_n2 + h * h) * cos(root_n_L)) / (h_root_n * (root_n2 + v2D * v2D));
}

/* Calculates the part of the i'th term of p_surv not dependent on t, without drift */
double GreensFunction1DRadAbs::p_survival_table_i_nov(uint i) const
{
   const double L(a_ - sigma_);
   const double h(k_ / D_);
   const double root_n(get_root(i));
   const double root_n2(root_n * root_n);
   const double root_n_r0_s(root_n * (r0_ - sigma_));
   const double root_n_L(root_n * L);
   const double h_root_n(h / root_n);
   return (h*sin(root_n_r0_s) + root_n*cos(root_n_r0_s)) / (L*(root_n2 + h*h) + h) * (h_root_n + sin(root_n_L) - h_root_n*cos(root_n_L));
}

/* Fills table with terms in the p_survival sum which don't depend on t */
void GreensFunction1DRadAbs::createPsurvTable(DoubleVector& table) const
{
   const uint root_nbr(rootList_size());
   uint i(static_cast<uint>(table.size()));

   if (v_ == 0.0)
   {
      while (i < root_nbr)
         table.emplace_back(p_survival_table_i_nov(i++));
   }
   else
   {
      while (i < root_nbr)
         table.emplace_back(p_survival_table_i_v(i++));
   }
}

/* Calculates the probability density of finding the particle at location r
   at time t. */
double GreensFunction1DRadAbs::prob_r(double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, t >= 0.0);
   THROW_UNLESS(std::invalid_argument, (r - sigma_) >= 0.0 && r <= a_ && (r0_ - sigma_) >= 0.0 && r0_ <= a_);

   const double L(a_ - sigma_);
   const double h((k_ + v_ / 2.0) / D_);

   const double vexpo(-v_*v_*t / D_ / 4.0 + v_*(r - r0_) / D_ / 2.0);

   // if there was no time change or zero diffusivity => no movement
   if (t == 0 || D_ == 0) return r == r0_ ? INFINITY : 0.0;    // the probability density function is a delta function

   // if r is at the absorbing boundary
   if (fabs(a_ - r) < GfCfg.EPSILON*L) return 0.0;

   double root_n, root_n_r_s;
   double sum = 0, term = 0, prev_term;

   const uint maxi(guess_maxi(t));
   calculate_n_roots(maxi);

   uint n = 0;
   do
   {
      if (n >= GfCfg.MAX_TERMS())
      {
         _log.warn() << "Too many terms needed for prob_r. N=" << n;
         break;
      }

      root_n = get_root(n);
      root_n_r_s = root_n*(r - sigma_);

      prev_term = term;
      term = Cn(root_n, t) * An(root_n) * (h*sin(root_n_r_s) + root_n*cos(root_n_r_s));
      sum += term;

      n++;
   } while (fabs(term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || fabs(prev_term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || n < GfCfg.MIN_TERMS());
   return 2.0*exp(vexpo)*sum;
}

/* Calculates the probability density of finding the particle at location z at
   timepoint t, given that the particle is still in the domain. */
double GreensFunction1DRadAbs::calcpcum(double r, double t) const
{
   return prob_r(r, t) / p_survival(t);
}

/* Calculates the total probability flux leaving the domain at time t
   This is simply the negative of the time derivative of the survival prob.
   at time t [-dS(t')/dt' for t'=t]. */
double GreensFunction1DRadAbs::flux_tot(double t) const
{
   const double vexpo(-v_*v_*t / 4.0 / D_ - v_*r0_ / 2.0 / D_);
   const double D2 = D_*D_;
   const double v2Dv2D = v_*v_ / 4.0 / D2;
   double sum = 0, term = 0, prev_term;

   const uint maxi(guess_maxi(t));
   calculate_n_roots(maxi);

   uint n = 0;
   do
   {
      if (n >= GfCfg.MAX_TERMS())
      {
         _log.warn() << "Too many terms needed for flux_tot. N=" << n;
         break;
      }

      double root_n = get_root(n);
      prev_term = term;
      term = (root_n * root_n + v2Dv2D) * Cn(root_n, t) * An(root_n) * Bn(root_n);
      n++;
      sum += term;
   } while (fabs(term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || fabs(prev_term / sum) > GfCfg.EPSILON*GfCfg.PDENS_TYPICAL || n < GfCfg.MIN_TERMS());

   return 2.0*D_*exp(vexpo)*sum;
}

/* Calculates the probability flux leaving the domain through the radiative
   boundary at time t */
double GreensFunction1DRadAbs::flux_rad(double t) const
{
   return k_ * prob_r(sigma_, t);
}

/* Calculates the flux leaving the domain through the radiative boundary as a
   fraction of the total flux. This is the probability that the particle left
   the domain through the radiative boundary instead of the absorbing
   boundary. */
double GreensFunction1DRadAbs::fluxRatioRadTot(double t) const
{
   return flux_rad(t) / flux_tot(t);
}

/* Determine which event has occurred, an escape or a reaction. Based on the
   fluxes through the boundaries at the given time. Beware: if t is not a
   first passage time you still get an answer! */
GreensFunction::EventKind GreensFunction1DRadAbs::drawEventType(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, t > 0.0);
   // if t=0 nothing has happened => no event

   const double L(a_ - sigma_);

   // if the radiative boundary is impermeable (k==0) or
   // the particle is at the absorbing boundary (at a) => IV_ESCAPE event
   if (k_ == 0 || fabs(a_ - r0_) < GfCfg.EPSILON*L) return EventKind::IV_ESCAPE;

   /* First check if we need to compare flux ratio's.
      If only one boundary is 'visible' to the particle, use this boundary as escape.*/
   const double distToa(a_ - r0_);
   const double distTos(r0_ - sigma_);
   const double maxDist(GfCfg.CUTOFF_H * (sqrt(2.0 * D_ * t) + fabs(v_ * t)));

   if (distToa > maxDist) //Absorbing boundary 'not in sight'.
   {
      if (distTos < maxDist) //Only radiation boundary 'in sight'.
         return EventKind::IV_REACTION;
   }
   else
   {
      if (distTos > maxDist) //Only absorbing boundary 'in sight'.
         return EventKind::IV_ESCAPE;
   }

   // Else the event is sampled from the flux ratio
   const double fluxratio(fluxRatioRadTot(t));

   return rnd > fluxratio ? EventKind::IV_ESCAPE : EventKind::IV_REACTION;
}

/* This function is needed to cast the math. form of the function
   into the form needed by the GSL root solver. */
double GreensFunction1DRadAbs::drawT_f(double t, void *p)
{
   struct drawT_params *params = static_cast<struct drawT_params *>(p);
   return params->rnd - params->gf->p_survival_table(t, params->psurvTable);
}

/* Draws the first passage time from the survival probability
   using an assistance function drawT_f that casts the math. function
   into the form needed by the GSL root solver. */
double GreensFunction1DRadAbs::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd < 1.0);

   const double L(a_ - sigma_);

   if (D_ == 0.0 || L == INFINITY) return INFINITY;

   if (rnd > 1 - GfCfg.EPSILON || L < 0.0 || fabs(a_ - r0_) < GfCfg.EPSILON*L) return 0.0;

   /* Find a good interval to determine the first passage time. */
   double t_guess;
   if (k_ != 0.0)
   {
      const double t_Abs(gsl_pow_2(a_ - r0_) / D_);
      const double t_Rad(D_ / (k_ * k_) + gsl_pow_2(r0_ - sigma_) / D_);

      t_guess = std::min(t_Abs, t_Rad);
   }
   else
   {
      t_guess = gsl_pow_2(a_ - r0_) / D_;
   }
   t_guess *= .1;

   // A different guess has to be made in case of nonzero drift to account for the displacement due to it
   // TODO: This does not work properly in this case yet, but we don't know why...
   // When drifting towards the closest boundary
   //if( (r0 >= a/2.0 && v > 0.0) || (r0 <= a/2.0 && v < 0.0) )	t_guess = sqrt(D*D_/(v*v*v*v)+dist*dist/(v*v)) - D_/(v*v);
   // When drifting away from the closest boundary
   //if( ( r0 < a/2.0 && v > 0.0) || ( r0 > a/2.0 && v < 0.0) )	t_guess = D/(v*v) - sqrt(D_*D_/(v*v*v*v)-dist*dist/(v*v));

   /* Set params structure. */
   DoubleVector psurvTable;
   struct drawT_params parameters = { this, psurvTable, rnd };

   /* Define the function for the rootfinder. */
   gsl_function F;
   F.function = &GreensFunction1DRadAbs::drawT_f;
   F.params = &parameters;

   double value(GSL_FN_EVAL(&F, t_guess));
   double low(t_guess);
   double high(t_guess);

   // scale the interval around the guess such that the function straddles
   if (value < 0.0)
   {
      // if the guess was too low
      do
      {
         if (fabs(high) >= t_guess * 1e10)
         {
            std::stringstream msg;
            msg << type_name() << ": couldn't adjust high. F(" << high << ")=" << value << "; " << dump();
            throw std::runtime_error(msg.str());
         }
         // keep increasing the upper boundary until the
         // function straddles
         high *= 10;
         value = GSL_FN_EVAL(&F, high);
      } while (value <= 0.0);
   }
   else
   {
      // if the guess was too high
      // initialize with 2 so the test below survives the first
      // iteration
      double value_prev(2);
      do
      {
         if (fabs(low) <= t_guess * 1e-10 || fabs(value - value_prev) < GfCfg.EPSILON)
         {
            _log.warn() << type_name() << ": couldn't adjust low. F(" << low << ")=" << value << ", " << dump();
            return low;
         }
         value_prev = value;
         // keep decreasing the lower boundary until the function straddles
         low *= 0.1;
         // get the accompanying value
         value = GSL_FN_EVAL(&F, low);
      } while (value >= 0.0);
   }

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(tscale_*GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction1DRadAbs::drawTime");
}

/* Returns c.d.f. for drawR */
double GreensFunction1DRadAbs::p_int_r_table(double r, double t, DoubleVector& table) const
{
   /* If not all boundaries are 'visible' to the particle
      use approximation. */
   const double distToa(a_ - r0_);
   const double distTos(r0_ - sigma_);
   const double maxDist(GfCfg.CUTOFF_H * (sqrt(2.0 * D_ * t) + fabs(v_ * t)));

   //TODO: include XI30 (c.d.f) with drift.
   if (distToa > maxDist) //Absorbing boundary 'not in sight'.
   {
      if (distTos > maxDist) //Radiation boundary 'not in sight'.
         return XI00(r, t, r0_, D_, v_); //free particle.
      else
      {
         if (k_ != 0.0 && v_ == 0.0)
            //Only radiation BCn.
            return XI30(r - sigma_, t, distTos, k_, D_, 0.0);
         else if (k_ == 0.0 && v_ == 0.0)
            //Only reflecting BCn.
            return XI20(r - sigma_, t, distTos, D_, 0.0);
      }
   }
   else
   {
      if (distTos > maxDist)
         //Only absorbing BCn.
         return XI10(a_ - r, t, distToa, D_, -v_);
   }

   const uint maxi(guess_maxi(t));
   const double vexpo(-v_*v_*t / 4.0 / D_ - v_*r0_ / 2.0 / D_);
   const double prefac(2.0*exp(vexpo));

   if (maxi >= GfCfg.MAX_TERMS())
      _log.warn() << "p_int_r_table: maxi was cut to MAX_TERMS for t=" << std::scientific << std::setprecision(16) << t;

   if (table.size() < maxi)
   {
      calculate_n_roots(maxi);
      create_p_int_r_Table(t, table);
   }

   double p = funcSum(std::bind(&GreensFunction1DRadAbs::p_int_r_i, this, std::placeholders::_1, r, t, table), GfCfg.MAX_TERMS());
   return prefac * p;
}

double GreensFunction1DRadAbs::p_int_r_i(uint i, double r, double t, DoubleVector& table) const
{
   const double h((k_ + v_ / 2.0) / D_);
   const double v2D(v_ / (2 * D_));
   const double costerm(k_ / D_);
   const double sinterm(h * v2D);
   const double expsigma(exp(sigma_*v2D));
   const double zs(r - sigma_);
   double root_n(get_root(i));
   double term((expsigma*costerm - exp(v2D*r) * (costerm*cos(root_n*zs) - (root_n + sinterm / root_n)*sin(root_n*zs))));
   return get_p_int_r_Table_i(i, t, table) * term;
}

/* Fills table for p_int_r of factors independent of r. */
void GreensFunction1DRadAbs::create_p_int_r_Table(double t, DoubleVector& table) const
{
   const uint root_nmbr(rootList_size());
   uint n(static_cast<uint>(table.size()));
   const double L(a_ - sigma_);
   const double h((k_ + v_ / 2.0) / D_);
   const double v2D(v_ / 2.0 / D_);
   const double v2Dv2D(v2D*v2D);

   double root_n2, root_n_r0_s, root_n;
   while (n < root_nmbr)
   {
      root_n = get_root(n);
      root_n2 = root_n * root_n;
      root_n_r0_s = root_n * (r0_ - sigma_);

      double term = exp(-D_*root_n2*t) * (root_n*cos(root_n_r0_s) + h*sin(root_n_r0_s)) / (L*(root_n2 + h*h) + h) * root_n / (root_n2 + v2Dv2D);
      table.emplace_back(term);
      n++;
   }
}

/* Function for GSL rootfinder of drawR. */
double GreensFunction1DRadAbs::drawR_f(double r, void *p)
{
   struct drawR_params *params = static_cast<struct drawR_params *>(p);
   return params->gf->p_int_r_table(r, params->t, params->table) - params->rnd;
}

/* Return new position */
double GreensFunction1DRadAbs::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd < 1.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   const double L(a_ - sigma_);

   if (t == 0.0 || (D_ == 0.0 && v_ == 0.0)) return r0_;
   if (a_ < 0.0) return 0.0;

   // the structure to store the numbers to calculate the numbers for 1-S
   DoubleVector pintTable;
   struct drawR_params parameters = { this, t, pintTable, rnd * p_survival(t) };

   // define gsl function for rootfinder
   gsl_function F;
   F.function = &GreensFunction1DRadAbs::drawR_f;
   F.params = &parameters;

   root_fsolver_wrapper solver;
   solver.set(&F, sigma_, a_);
   return solver.findRoot(L*GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction1DRadAbs::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction1DRadAbs::dump() const
{
   std::ostringstream ss;
   ss << "D=" << D_ << ", sigma=" << sigma_ << ", a=" << a_ << ", r0=" << r0_ << ", k=" << k_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
