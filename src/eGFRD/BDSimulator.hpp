#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "ParticleSimulator.hpp"
#include "BDPropagator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class BDSimulator : public ParticleSimulator
{
public:
   using base_type = ParticleSimulator;
   using double_pair = std::pair<const double, const double>;

   // --------------------------------------------------------------------------------------------------------------------------------

   BDSimulator() = delete;

   explicit BDSimulator(World& world, ReactionRuleCollection& reaction_rules, RandomNumberGenerator& rng) noexcept
      : base_type(world, reaction_rules, rng), dt_factor_(0.1), num_retries_(1), reaction_length_factor_(GfrdCfg.MULTI_SHELL_FACTOR)
   {
      dt_ = determine_dt();
      reaction_length_ = determine_reaction_length();
   }


   // --------------------------------------------------------------------------------------------------------------------------------

   double dt_factor() const { return dt_factor_; }

   double get_reaction_length() const { return reaction_length_; }

   double get_num_retries() const { return num_retries_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Returns the upper limit of dt calculated using the  */
   double determine_dt() const
   {
      auto maxD_minr = Dmax_minsigma();
      double k_max = get_max_rate();
      double r_min = maxD_minr.second;

      double Pacc_max = 0.01; // maximum value of the acceptance probability.
      double tau_D = 2.0 * r_min * r_min / maxD_minr.first; // ~step over particle - time.
      if (k_max <= 0.0) return dt_factor_ * tau_D;

      double dt_temp = 2.0 * Pacc_max * reaction_length_factor_ * r_min / k_max;
      return std::min(dt_temp, dt_factor_ * tau_D); // dt_factor * tau_D is upper limit of dt.
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Returns the reaction length calculated from the smallest particle radius in the multi. */
   double determine_reaction_length() const
   {
      auto dmax_rmin = Dmax_minsigma();
      return reaction_length_factor_ * dmax_rmin.second;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Sets the reaction factor and reaction length and calculates the time step size and the reaction length. */
   void set_reaction_length_factor(double rl, double dt)
   {
      dt_factor_ = dt;
      dt_ = determine_dt();

      reaction_length_factor_ = rl;
      reaction_length_ = determine_reaction_length();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Returns largest diffusion constant and smallest particle radius in the multi. */
   double_pair Dmax_minsigma() const
   {
      double D_max = 0.0;
      double radius_min = std::numeric_limits<double>::max();

      for (auto&& st : world_.get_species())
      {
         SpeciesType sid = st.second;
         if (D_max < sid.D()) D_max = sid.D();
         if (radius_min > sid.radius()) radius_min = sid.radius();
      }

      return double_pair(D_max, radius_min);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Returns the largest _1D_ intrinsic reaction rate in the world. */
   double get_max_rate() const
   {
      double k_max = 0.0;
      int i = 0;
      int j = 0;

      // surface rates not divided by 2 to compensate for double reaction attempts.
      k_max *= 2.0;
      for (auto s0 : world_.get_species())
      {
         j = 0;
         for (auto s1 : world_.get_species())
         {
            if (j++ < i) continue;

            auto rules = reaction_rules_.query_reaction_rules(s0.second.id(), s1.second.id());
            if (rules.size() == 0) continue;

            i++;
            for (auto rr : rules)
            {
               double k = 0.001;
               if (k_max < k) k_max = k;
            }
         }
      }
      return k_max;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void step()
   {
      _step(dt_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // step upto 
   bool step(double upto)
   {
      auto lt = upto - time_;
      if (lt <= 0.0) return false;

      if (lt > dt_) { _step(lt); }
      else { _step(dt_); }

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Runs the BD propagator with stepsize dt untill the queue is empty. */
   void _step(double dt)
   {
      std::vector<ParticleID> particles;
      class VolumeClearer vc;

      auto propagator = BDPropagator(world_, particles, reaction_rules_, rng_, dt, 0, vc, 0.0, nullptr, 0);
      propagator.step_all();
      ++num_steps_;
      time_ += dt;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:
   double dt_factor_;
   int num_retries_;
   double reaction_length_;
   double reaction_length_factor_;
};

#endif
