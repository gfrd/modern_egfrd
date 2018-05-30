#include <cmath>
#include <algorithm>
#include <sstream>
#include "freeFunctions.hpp"
#include "helperFunctionsGf.hpp"
#include "GreensFunction1DAbsSinkAbs.hpp"
#include "Logger.hpp"
#include <iomanip>

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("GreensFunction1DAbsSinkAbs");

// --------------------------------------------------------------------------------------------------------------------------------

/* This is the appropriate definition of the function defining
   the roots of our Green's functions in GSL.
   Later needed by the rootfinder. */
double GreensFunction1DAbsSinkAbs::root_f(double x, double Lm_L, double h) const
{
   // L   = Lr + Ll
   // h    = L * k / (2 * D)
   // L_Lm = Lr + Ll / Lr - Ll
   // x    = q * L
   return x * std::sin(x) + h * (std::cos(x * Lm_L) - std::cos(x));
}

/* return the n + 1'th root */
double GreensFunction1DAbsSinkAbs::get_root(uint n) const
{
   if (n >= rootList_.size())
      calculate_n_roots(n + 1);
   return rootList_[n];
}

/* Calculates the first n roots of root_f */
void GreensFunction1DAbsSinkAbs::calculate_n_roots(uint n) const
{
   uint i = rootList_size();
   if (n <= i) return;

   // lift the const veil and fill the table
   auto ucthis = const_cast<GreensFunction1DAbsSinkAbs*>(this);
   ucthis->fill_table_to_n(i, n);
}

void GreensFunction1DAbsSinkAbs::fill_table_to_n(uint i, uint n)
{
   double L = Lr_ + Ll_;
   double Lm = Lr_ - Ll_;
   double Lm_L = Lm / L;
   double h = k_ * L / (2 * D_);

   // reserve space
   rootList_.reserve(n);

   /* Define the root function. */
   auto f = [Lm_L, h, this](double x) { return root_f(x, Lm_L, h); };
   gsl_lambda<decltype(f)> F(f);


   /* define and ad a new solver type brent */
   root_fsolver_wrapper solver;

   /* If this is the first run, set parameters.*/
   if (i == 0)
   {
      luparams.h = h;
      luparams.Lm_L = Lm_L;
      luparams.long_period = std::max(L / Lr_ * M_PI, L / Ll_ * M_PI);
      luparams.short_period = std::min(L / Lr_ * M_PI, L / Ll_ * M_PI);
      luparams.last_long_root = 0.0;
      luparams.last_short_root = 0.0;
   }

   /* Find all the roots up to the nth */
   while (i++ < n)
   {
      DoublePair lower_upper_pair = get_lower_and_upper();
      solver.set(&F, lower_upper_pair.first, lower_upper_pair.second);
      double root_i = solver.findRoot(GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction1DAbsSinkAbs::root_f");

      ASSERT(root_i > std::max(luparams.last_long_root, luparams.last_short_root) - GfCfg.EPSILON);
      rootList_.emplace_back(root_i / L);
      if (luparams.last_was_long)
         luparams.last_long_root = root_i;
      else
         luparams.last_short_root = root_i;
   }
}

/* returns two points on the x-axis which straddle the next root. */
DoublePair GreensFunction1DAbsSinkAbs::get_lower_and_upper()
{
   double root_n = std::max(luparams.last_long_root, luparams.last_short_root);
   double safety = 0.75;

   double lower, upper, next_root_est, left_offset, right_offset;
   double last_root = root_n == 0.0 ? M_PI : root_n;

   if (luparams.h / last_root < 1)
   {
      right_offset = M_PI;
      next_root_est = root_n + M_PI;
   }
   else
   {
      double next_root_long = luparams.last_long_root + luparams.long_period;
      double next_root_short = luparams.last_short_root + luparams.short_period;

      if (next_root_long < next_root_short)
      {
         next_root_est = next_root_long;
         right_offset = std::min(next_root_short - next_root_est, luparams.long_period);
         luparams.last_was_long = true;
      }
      else
      {
         next_root_est = next_root_short;
         right_offset = std::min(next_root_long - next_root_est, luparams.short_period);
         luparams.last_was_long = false;
      }
   }

   left_offset = next_root_est - root_n - 1000 * GfCfg.EPSILON;
   lower = next_root_est - left_offset;
   upper = next_root_est + safety * right_offset;

   double f_lower = root_f(lower, luparams.Lm_L, luparams.h);
   double f_upper = root_f(upper, luparams.Lm_L, luparams.h);

   /* set the parity operator for the next root:
      +1 for an even root.
      -1 for an odd root. */
   int parity_op = 2 * (rootList_size() % 2) - 1;

   /* f_lower must have correct sign. */
   if (f_lower * parity_op > 0)
   {
      _log.warn() << "f(lower) has wrong sign at root# " << (rootList_size() + 1) << ", for h=" << std::scientific << std::setprecision(16) << luparams.h << ", Lm/L=" << luparams.Lm_L;
      _log.warn() << "f_low( " << lower << " )=" << f_lower << " , f_high(" << upper << ")=" << f_upper;
   }

   /* If the parity is incorrect -> correct it */
   if (f_upper * parity_op < 0)
   {
      int cntr(0);

      double delta = .1 * safety * right_offset;
      double save_upper = upper;

      /* Assuming the upper point has overshoot the straddle region
         subtract from the upper limit, until endpoints do straddle.
         */
      while (f_upper * parity_op < 0 && cntr++ < 10)
      {
         //Correct for overshoot.
         upper -= delta;
         f_upper = root_f(upper, luparams.Lm_L, luparams.h);
      }

      /* If upper point still doesn't straddle the root, increase the old estimate
         of upper, until it does straddle*/
      if (cntr >= 10)
      {
         cntr = 0;
         upper = save_upper;
         f_upper = root_f(upper, luparams.Lm_L, luparams.h);

         while (f_upper * parity_op < 0 && cntr++ < 10)
         {
            //Correct for undershoot.
            upper += delta;
            f_upper = root_f(upper, luparams.Lm_L, luparams.h);
         }

         // Still no straddle? => Crash is imminent.
         if (cntr >= 10)
         {
            auto io = _log.warn();
            io << "Failed to straddle root# " << rootList_size() + 1;
            io << " next root est.=" << next_root_est << " with stepsize=" << delta << ", ll=" << luparams.last_long_root << ", ls=" << luparams.last_short_root;
            io << "f_low(" << lower << ")= " << f_lower << ", f_high(" << upper << ")=" << f_upper;
         }
      }
   }

   return DoublePair(lower, upper);
}

/* returns a guess for the number of terms needed for
   the greensfunction to converge at time t */
uint GreensFunction1DAbsSinkAbs::guess_maxi(double t) const
{
   const uint safety(2);
   if (!std::isfinite(t)) return safety;

   double root0(get_root(0));
   double Dt(D_ * t);
   double thr(exp(-Dt * root0 * root0) * GfCfg.EPSILON * 1e-1);
   if (thr <= 0.0) return GfCfg.MAX_TERMS();

   double max_root(sqrt(root0 * root0 - log(thr) / Dt));
   const uint maxi(safety + static_cast<uint>(max_root * (Lr_ + Ll_) / M_PI));
   return std::min(maxi, GfCfg.MAX_TERMS());
}

/* Standard form of the greensfunction without numerator */
inline double GreensFunction1DAbsSinkAbs::p_exp_den_i(double t, double root_i, double root_i2) const
{
   return exp(-D_ * root_i2 * t) / p_denominator_i(root_i);
}

/* Denominator of the greensfunction. */
inline double GreensFunction1DAbsSinkAbs::p_denominator_i(double root_n) const
{
   double Lm(Lr_ - Ll_);
   double L(Lr_ + Ll_);

   double term1(root_n * L * std::cos(root_n * L) + std::sin(root_n * L));
   double term2(L * std::sin(root_n * L) - Lm * std::sin(root_n * Lm));

   return D_ * term1 + k_ / 2. * term2;
}

double GreensFunction1DAbsSinkAbs::p_survival(double t) const
{
   DoubleVector table;
   return p_survival_table(t, table);
}

/* Calculates survival probability using a table.
   Switchbox for which greensfunction to use. */
double GreensFunction1DAbsSinkAbs::p_survival_table(double t, DoubleVector& psurvTable) const
{
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   if (t == 0.0 || D_ == 0.0) return 1.0; //particle can't escape.

   /* First check if we need full solution.
      Else we use approximation. */
   double maxDist(GfCfg.CUTOFF_H * sqrt(2.0 * D_ * t));
   // dist to nearest absorbing boundary.
   double distToAbs(std::min(a_ - r0_, r0_ - sigma_));

   if (L0_ > maxDist) //Sink not in sight
   {
      if (distToAbs > maxDist)
         return 1.0;
      else
         return XS10(t, distToAbs, D_, 0.0);
   }
   else
   {
      if (distToAbs > maxDist) //Only sink in sight.
      {
         return XS030(t, L0_, k_, D_);
      }
   }

   const uint maxi(guess_maxi(t));

   if (maxi == GfCfg.MAX_TERMS())
      _log.info() << "drawT: maxi was cut to MAX_TERMS for t=" << std::scientific << std::setprecision(16) << t;

   if (psurvTable.size() < maxi)
   {
      calculate_n_roots(maxi);  // this updates the table
      createPsurvTable(psurvTable);
   }

   double p = funcSum_all(std::bind(&GreensFunction1DAbsSinkAbs::p_survival_i, this, std::placeholders::_1, t, psurvTable), maxi);
   return p;
}

/* Calculates the i'th term of the p_survival sum. */
double GreensFunction1DAbsSinkAbs::p_survival_i(uint i, double t, const DoubleVector& table) const
{
   double root_i(get_root(i));
   return exp(-D_ * t * gsl_pow_2(root_i)) * table[i];
}

/* Calculates the part of the i'th term of p_surv not dependent on t */
double GreensFunction1DAbsSinkAbs::p_survival_table_i(double root_i) const
{
   double L(Lr_ + Ll_);
   double LrmL0(Lr_ - L0_);
   double term1(std::sin(root_i * L) - std::sin(root_i * LrmL0) - std::sin(root_i * (Ll_ + L0_)));
   double term2(std::sin(root_i * Lr_) - std::sin(root_i * L0_) - std::sin(root_i * LrmL0));

   double numerator(D_ * term1 + k_ * std::sin(root_i * Ll_) * term2 / root_i);
   numerator *= 2.0;
   return numerator / p_denominator_i(root_i);
}

/* Fills table with terms in the p_survival sum which don't depend on t. */
void GreensFunction1DAbsSinkAbs::createPsurvTable(DoubleVector& table) const
{
   uint const root_nbr(rootList_size());
   uint i(static_cast<uint>(table.size()));
   table.reserve(root_nbr);

   while (i < root_nbr)
      table.emplace_back(p_survival_table_i(get_root(i++)));
}

/* Returns i'th term of prob_r (domain containing r0) function */
double GreensFunction1DAbsSinkAbs::prob_r_r0_i(uint i, double rr, double t) const
{
   double root_i = get_root(i);
   double L02 = getL0();       // copy!
   double rr2 = rr;

   /* Check if r is left or right of the starting position r0.
      If so, interchange rr with L0. */
   if (rr < L02)
   {
      rr2 = L02;
      L02 = rr;
   }

   double LlpL0 = Ll_ + L02;
   double Lrmrr = Lr_ - rr2;

   double numerator = D_ * root_i * std::sin(root_i * LlpL0) + k_ * std::sin(root_i * Ll_) * std::sin(root_i * L02);
   numerator *= std::sin(root_i * Lrmrr);
   return -2.0 * p_exp_den_i(t, root_i, gsl_pow_2(root_i)) * numerator;
}

/* Returns i'th term of prob_r (domain not containing r0) function */
double GreensFunction1DAbsSinkAbs::prob_r_nor0_i(uint i, double rr, double t) const
{
   double root_i(get_root(i));
   double LrmL0(Lr_ - L0_);
   double Llprr(Ll_ + rr);
   double numerator(D_ * root_i * std::sin(root_i * Llprr) * std::sin(root_i * LrmL0));
   return -2.0 * p_exp_den_i(t, root_i, gsl_pow_2(root_i)) * numerator;
}

/* Calculates the probability density of finding the particle at location r at time t. */
double GreensFunction1DAbsSinkAbs::prob_r(double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, t >= 0.0);
   THROW_UNLESS(std::invalid_argument, (r - sigma_) >= 0.0 && r <= a_ && (r0_ - sigma_) >= 0.0 && r0_ <= a_);

   double L(Lr_ + Ll_);

   // if there was no time change or zero diffusivity => no movement
   if (t == 0 || D_ == 0) return r == r0_ ? INFINITY : 0.0;            // the probability density function is a delta function

   // if r is at one of the the absorbing boundaries
   if (fabs(a_ - r) < GfCfg.EPSILON * L || fabs(r - sigma_) < GfCfg.EPSILON * L) return 0.0;

   double rr(r0_ - rsink_ >= 0 ? r - rsink_ : rsink_ - r);

   double p;

   /* Determine wether rr lies in the same sub-domain as r0. A different function is calculated when this is the case. */
   if (rr >= 0)
   {
      p = funcSum(std::bind(&GreensFunction1DAbsSinkAbs::prob_r_r0_i, this, std::placeholders::_1, rr, t), GfCfg.MAX_TERMS());
   }
   else
   {
      p = funcSum(std::bind(&GreensFunction1DAbsSinkAbs::prob_r_nor0_i, this, std::placeholders::_1, rr, t), GfCfg.MAX_TERMS());
   }

   return p;
}

/* Calculates the probability density of finding the particle at location r at
   time t, given that the particle is still in the domain. */
double GreensFunction1DAbsSinkAbs::calcpcum(double r, double t) const
{
   return prob_r(r, t) / p_survival(t);
}

/* Function returns flux at absorbing boundary sigma. */
double GreensFunction1DAbsSinkAbs::flux_leaves(double t) const
{
   if (t == 0 || D_ == 0) return 0.0;
   //const uint maxi(guess_maxi(t));

   if (r0_ >= rsink_)
      return flux_abs_Ll(t);
   return -flux_abs_Lr(t);
}

/* Function returns flux at absorbing boundary a. */
double GreensFunction1DAbsSinkAbs::flux_leavea(double t) const
{
   if (t == 0 || D_ == 0) return 0.0;

   //const uint maxi(guess_maxi(t));
   if (r0_ < rsink_)
      return -flux_abs_Ll(t);
   return flux_abs_Lr(t);
}

/* Calculates the total probability flux leaving the domain at time t
   This is simply the negative of the time derivative of the survival prob.
   at time t [-dS(t')/dt' for t'=t]. */
double GreensFunction1DAbsSinkAbs::flux_tot(double t) const
{
   if (t == 0 || (D_ == 0 && r0_ != rsink_)) return 0.0;

   double p;
   p = funcSum(std::bind(&GreensFunction1DAbsSinkAbs::flux_tot_i, this, std::placeholders::_1, t), GfCfg.MAX_TERMS());
   return D_ * p;
}

/* Calculates i'th term of total flux leaving at time t. */
double GreensFunction1DAbsSinkAbs::flux_tot_i(uint i, double t) const
{
   double root_i(get_root(i));
   double root_i2(gsl_pow_2(root_i));
   return root_i2 * exp(-D_ * t * root_i2) * p_survival_table_i(root_i);
}

/* Flux leaving through absorbing boundary from sub-domain containing r0. */
double GreensFunction1DAbsSinkAbs::flux_abs_Lr(double t) const
{
   double p;
   p = funcSum(std::bind(&GreensFunction1DAbsSinkAbs::flux_abs_Lr_i, this, std::placeholders::_1, t), GfCfg.MAX_TERMS());
   return -D_ * 2 * p;
}

/* Calculates the i'th term of the flux at Lr. */
double GreensFunction1DAbsSinkAbs::flux_abs_Lr_i(uint i, double t) const
{
   double root_i(get_root(i));
   double LlpL0(Ll_ + L0_);

   double numerator(k_ * std::sin(root_i * Ll_) * std::sin(root_i * L0_) + D_ * root_i * std::sin(root_i * LlpL0));
   numerator *= root_i;
   return p_exp_den_i(t, root_i, gsl_pow_2(root_i)) * numerator;
}

/* Flux leaving through absorbing boundary from sub-domain not containing r0. */
double GreensFunction1DAbsSinkAbs::flux_abs_Ll(double t) const
{
   double D2(gsl_pow_2(D_));
   double p;
   p = funcSum(std::bind(&GreensFunction1DAbsSinkAbs::flux_abs_Ll_i, this, std::placeholders::_1, t), GfCfg.MAX_TERMS());
   return 2 * D2 * p;
}

/* Calculates the i'th term of the flux at Ll. */
double GreensFunction1DAbsSinkAbs::flux_abs_Ll_i(uint i, double t) const
{
   double root_i(get_root(i));
   double root_i2(gsl_pow_2(root_i));
   double LrmL0(Lr_ - L0_);

   double numerator(root_i2 * std::sin(root_i * LrmL0));
   return p_exp_den_i(t, root_i, root_i2) * numerator;
}

/* Calculates the probability flux leaving the domain through the sink
   at time t */
double GreensFunction1DAbsSinkAbs::flux_sink(double t) const
{
   if (t == 0 || (D_ == 0 && r0_ != rsink_)) return 0.0;
   return k_ * prob_r(rsink_, t);
}

/* Determine which event has occurred, an escape or a reaction. Based on the
   fluxes through the boundaries and the sink at the given time. */
GreensFunction::EventKind GreensFunction1DAbsSinkAbs::drawEventType(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, rnd >= 0.0 && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, t > 0.0);

   double L(a_ - sigma_);

   /* If the sink is impermeable (k==0) or
      the particle is at one the absorbing boundaries (sigma or a) => IV_ESCAPE event */
   if (k_ == 0 || fabs(a_ - r0_) < GfCfg.EPSILON * L || fabs(sigma_ - r0_) < GfCfg.EPSILON * L) return EventKind::IV_ESCAPE;

   /* The event is sampled from the flux ratios.
      Two possibilities:
      (1) Leave through left or right boundary - IV_ESCAPE
      (2) Leave through sink - IV_REACTION
      */

      /* First check if we need to compare flux ratio's.
         If particle is near only one boundary or sink, this is the escape event. */
   double maxDist(GfCfg.CUTOFF_H * sqrt(2.0 * D_ * t));
   // dist to nearest absorbing boundary.
   double distToAbs(std::min(a_ - r0_, r0_ - sigma_));

   if (L0_ > maxDist) //Sink not in sight
   {
      if (distToAbs < maxDist) return EventKind::IV_ESCAPE;
   }
   else
   {
      if (distToAbs > maxDist) return EventKind::IV_REACTION; //Only sink in sight.
   }

   /* Already fill rootList with needed roots. */
   const uint maxi(guess_maxi(t));

   if (maxi == GfCfg.MAX_TERMS()) _log.info() << "drawEventType: maxi was cut to MAX_TERMS for t=" << std::scientific << std::setprecision(16) << t;

   calculate_n_roots(maxi);
   rnd *= flux_tot(t);
   double p_sink(flux_sink(t));
   return rnd < p_sink ? EventKind::IV_REACTION : EventKind::IV_ESCAPE;
}

double GreensFunction1DAbsSinkAbs::drawTime(double rnd) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd < 1.0);

   double L(Lr_ + Ll_);

   if (D_ == 0.0 || L == INFINITY) return INFINITY;

   if (rnd > (1 - GfCfg.EPSILON) || L < 0.0 || fabs(a_ - r0_) < GfCfg.EPSILON * L) return 0.0;

   /* Find a good interval to determine the first passage time.
      First we get the distance to one of the absorbing boundaries or the sink. */
   double dist(std::min(Lr_ - L0_, Ll_ + L0_));
   double t_guess;

   /* For a particle at contact, the time for the particle to be absorbed by the
      radiation boundary is 4 D / (k*k). If the particle is at a distance x0 from the
      radiation boundary, we approximate the guess for the next event-time by:
      t_guess = D / (k*k) + x0 * x0 / D_.
      */
   double t_Abs(gsl_pow_2(dist) / D_);
   double t_Rad(4 * D_ / (k_ * k_) + gsl_pow_2(L0_) / D_);

   t_guess = std::min(t_Abs, t_Rad);
   t_guess *= .1;

   /* the structure to store the numbers to calculate the numbers for 1-S */
   DoubleVector psurvTable;
   auto f = [rnd, &psurvTable, this](double t) { return rnd - p_survival_table(t, psurvTable); };
   gsl_lambda<decltype(f)> F(f);

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
            _log.warn() << "drawTime Couldn't adjust low. F(" << low << ")=" << value;
            /*
              std::cerr << "GF1DSink::drawTime Couldn't adjust low. F(" << low << ")="
              << value << " t_guess: " << t_guess << " diff: "
              << (value - value_prev) << " value: " << value
              << " value_prev: " << value_prev << " rnd: "
              << rnd << std::endl;
              */
            return low;
         }
         value_prev = value;
         // keep decreasing the lower boundary until the function straddles
         low *= 0.1;
         // get the accompanying value
         value = GSL_FN_EVAL(&F, low);
      } while (value >= 0.0);
   }

   /* find the intersection on the y-axis between the random number and the function */
   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(tscale_*GfCfg.EPSILON, GfCfg.EPSILON, "GreensFunction1DAbsSinkAbs::drawTime");
}

/* Returns the c.d.f. with respect to the position at time t.
   - Also a switchbox for which GF to choose. */
double GreensFunction1DAbsSinkAbs::p_int_r_table(double r, double t, DoubleVector& table) const
{
   double p;

   /* when r0 lies left of the sink, mirror the domain around rsink
      : rr -> -rr. */
   double rr(r0_ - rsink_ >= 0 ? r - rsink_ : rsink_ - r);

   /* First check if we need full solution.
      Else we use approximation. */
   double maxDist(GfCfg.CUTOFF_H * sqrt(2.0 * D_ * t));
   double distToa(a_ - r0_);
   double distTos(r0_ - sigma_);

   if (L0_ > maxDist) //Sink not in sight
   {
      if (distTos > maxDist)
      {
         if (distToa > maxDist)
            return XI00(r, t, r0_, D_, 0.0);  //Nothing in sight.
         else
            return XI10(a_ - r, t, distToa, D_, 0.0); //Only a boundary in sight.
      }
      else
      {
         if (distToa > maxDist)
            return XI10(r - sigma_, t, distTos, D_, 0.0); //Only s boundary in sight.
      }
   }
   else
   {
      if (distToa > maxDist && distTos > maxDist) //Only sink in sight.
      {
         return XI030(rr, t, L0_, k_, D_);
      }
   }

   const uint maxi(guess_maxi(t));

   if (maxi == GfCfg.MAX_TERMS()) _log.warn() << "p_int_r_table: maxi was cut to MAX_TERMS for t=" << t;

   if (table.size() < maxi)
   {
      calculate_n_roots(maxi);  // this updates the table
      create_p_int_r_Table(t, table);
   }

   /* Determine in which part of the domain rr lies, and
      thus which function to use. */
   double(GreensFunction1DAbsSinkAbs::*p_int_r_i) (uint, double, double, DoubleVector&) const;

   if (rr <= 0)
      p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_leftdomain;
   else if (rr < L0_)
      p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA;
   else
      p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB;

   p = funcSum(std::bind(p_int_r_i, this, std::placeholders::_1, rr, t, table), GfCfg.MAX_TERMS());

   return 2.0 * p;
}

double GreensFunction1DAbsSinkAbs::p_int_r(double r, double t) const
{
   THROW_UNLESS(std::invalid_argument, r >= sigma_ && r <= a_);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   DoubleVector table;
   return p_int_r_table(r, t, table);
}

void GreensFunction1DAbsSinkAbs::create_p_int_r_Table(double t, DoubleVector& table) const
{
   const uint root_nbr(rootList_size());
   uint i(static_cast<uint>(table.size()));

   while (i < root_nbr)
   {
      double root_i(get_root(i));
      table.emplace_back(p_exp_den_i(t, root_i, gsl_pow_2(root_i)));
      i++;
   }
}

//Integrated Greens function for rr part of [-Ll, 0]
double GreensFunction1DAbsSinkAbs::p_int_r_leftdomain(uint i, double rr, double t, DoubleVector& table) const
{
   double root_i(get_root(i));
   double LrmL0(Lr_ - L0_);
   double Llprr(Ll_ + rr);
   double temp(D_ * std::sin(root_i * LrmL0) * (std::cos(root_i * Llprr) - 1.0));
   return get_p_int_r_Table_i(i, t, table) * temp;
}

//Integrated Greens function for rr part of (0, L0]
double GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA(uint i, double rr, double t, DoubleVector& table) const
{
   double root_i(get_root(i));
   double LrmL0(Lr_ - L0_);
   double Llprr(Ll_ + rr);
   double root_i_rr(root_i * rr);
   double temp(D_ * (std::cos(root_i * Llprr) - 1.0) + k_ / root_i * (std::cos(root_i_rr) - 1.0) * std::sin(root_i * Ll_));
   return get_p_int_r_Table_i(i, t, table) * std::sin(root_i * LrmL0) * temp;
}

//Integrated Greens function for rr part of (L0, Lr]
double GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB(uint i, double rr, double t, DoubleVector& table) const
{
   double root_i(get_root(i));
   double L(Lr_ + Ll_);
   double LrmL0(Lr_ - L0_);
   double Lrmrr(Lr_ - rr);
   double LlpL0(Ll_ + L0_);

   double term1(std::sin(root_i * L) - std::sin(root_i * LrmL0) - std::sin(root_i * LlpL0) * std::cos(root_i * Lrmrr));
   double term2(std::sin(root_i * Lr_) - std::sin(root_i * LrmL0) - std::sin(root_i * L0_) * std::cos(root_i * Lrmrr));
   double temp(D_ * term1 + k_ * std::sin(root_i * Ll_) * term2 / root_i);

   return get_p_int_r_Table_i(i, t, table) * temp;
}

// --------------------------------------------------------------------------------------------------------------------------------

double GreensFunction1DAbsSinkAbs::drawR(double rnd, double t) const
{
   THROW_UNLESS(std::invalid_argument, 0.0 <= rnd && rnd <= 1.0);
   THROW_UNLESS(std::invalid_argument, t >= 0.0);

   double L(Lr_ + Ll_);

   if (t == 0.0 || D_ == 0.0) return r0_;// the trivial case
   if (L < 0.0) return 0.0;// if the domain had zero size
   if (rnd <= GfCfg.EPSILON) return sigma_;
   if (rnd >= (1 - GfCfg.EPSILON)) return a_;

   // the structure to store the numbers to calculate r.
   DoubleVector pintTable;
   double rndpsurf = rnd * p_survival(t);
   auto f = [t, &pintTable, rndpsurf, this](double r) { return p_int_r_table(r, t, pintTable) - rndpsurf; };
   gsl_lambda<decltype(f)> F(f);

   root_fsolver_wrapper solver;
   solver.set(&F, sigma_, a_);
   return solver.findRoot(GfCfg.EPSILON*L, GfCfg.EPSILON, "GreensFunction1AbsSinkAbs::drawR");
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string GreensFunction1DAbsSinkAbs::dump() const
{
   std::ostringstream ss;
   ss << "D=" << D_ << ", sigma=" << sigma_ << ", a=" << a_ << ", r0=" << r0_ << ", rsink=" << rsink_ << ", k=" << k_;
   return ss.str();
}

// --------------------------------------------------------------------------------------------------------------------------------
