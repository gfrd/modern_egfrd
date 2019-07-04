#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "Shell.hpp"
#include "Domain.hpp"
#include "Particle.hpp"
#include "DefsEgfrd.hpp"
#include "ReactionRuleCollection.hpp"
#include "BDPropagator.hpp"
#include "ReactionRecorder.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Multi : public Domain
{

public:
   using shell_matrix_type = MatrixSpace<Shell, ShellID, World::particle_matrix_type::SizeX, World::particle_matrix_type::SizeY, World::particle_matrix_type::SizeZ>;

   using shell_id_pair = Shell::shell_id_pair;
   using particle_id_pair = Particle::particle_id_pair;
   using particle_id_pair_generator = agi::iteratorRange<particle_id_pair>;
   using particle_id_list = const std::vector<ParticleID>&;
   using species_map = std::unordered_set<SpeciesTypeID>;

   explicit Multi(const DomainID id, World& world, shell_matrix_type& shell_matrix, ReactionRuleCollection& reactions, RandomNumberGenerator& rng, reaction_recorder* rrec, double starttime) noexcept :
      Domain(id), world_(world), reactions_(reactions), rng_(rng), rrec_(rrec), shell_map_(), shell_container_(shell_matrix), reaction_length_(), start_time_(starttime) { }

   size_t num_shells() const override { return particles_.size(); }

   Multiplicity multiplicity() const override { return Multiplicity::MULTI; }

   shell_id_pair_ref get_shell() const override
   {
      if (num_shells() > 1) throw unsupported(make_string() << "For multi's use get_shell_list.");
      auto sid = shell_map_.begin();
      return std::make_pair(sid->second, std::ref<const Shell>(shell_container_.find(sid->second)->second));
   }

   std::vector<shell_id_pair_ref> get_shell_list() const override
   {
      auto shells = std::vector<shell_id_pair_ref>();
      shells.reserve(shell_map_.size());
      for (auto& sid : shell_map_)
         shells.emplace_back(std::make_pair(sid.second, std::ref<const Shell>(shell_container_.find(sid.second)->second)));
      return shells;
   }

   std::string as_string() const override { return std::string("empty"); }

   const char* type_name() const override { return "Multi"; }

   double start_time() const { return start_time_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void initialize(double t)
   {
      last_time_ = t;
      set_dt_and_reaction_length(GfrdCfg.DEFAULT_STEPSIZE_FACTOR, GfrdCfg.BD_DT_HARDCORE_MIN);
      eventType_ = EventType::MULTI_DIFFUSION;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_shell(shell_id_pair sid_pair, const ParticleID pid)
   {
      shell_map_[pid] = sid_pair.first;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void queue_particle(particle_id_pair pid_pair)
   {
      particles_.emplace_back(pid_pair.first);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_dt_and_reaction_length(double step_size_factor, double dt_min = -1.0)
   {
      const auto& species = get_species();

      double k_max = get_max_rate(species);
      double r_min = get_min_radius(species);
      reaction_length_ = step_size_factor * r_min;

      // The following gives us the typical timescale to travel the length step_size_factor * r_min,
      // taking into account all species in the Multi and automatically comparing whether their movement
      // is dominated by convection or diffusion, respectively. tau_Dv is the upper limit of time step dt.
      auto tau_Dv = get_min_tau_Dv(step_size_factor * r_min, species);

      auto Pacc_max = 0.1;  // Maximum allowed value of the acceptance probability. // TESTING was 0.01
                            // This should be kept very low (max. 0.01), otherwise the approximation of
                            // treating the reaction as two sequential attempts (first move, then react)
                            // might break down!

      // Here it is assumed that the reaction length is linear in any dimension,
      // which requires it to be very small!
      auto dt_motion = tau_Dv; // timescale of motion
      auto dt_reaction = k_max > 0 ? 2.0 * Pacc_max * reaction_length_ / k_max : 0; // timescale of reacting

      if(dt_reaction == dt_motion)
      { Log("GFRD").debug() << "BD time step set by timescale of motion";  }
      else
      { Log("GFRD").debug() << "BD time step set by largest reaction rate, k_max = " << k_max << ", r_min = " << r_min;}

      dt_ = dt_reaction > 0 ? std::min(dt_reaction, dt_motion) : dt_motion;
      if (dt_ < dt_min)
      {
         dt_ = dt_min;
         Log("EGFRD").warn() << "Setting Multi time step to hard-coded minimum, dt = " << dt_ << ", to allow continuing the simulation. This breaks detailed balance!";
      }

   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void step()
   {
      int num_retries = 1; //TODO dissociation_retry_moves ???

      volume_clearer vc(*this);
      auto propagator = BDPropagator(world_, particles_, reactions_, rng_, dt_, num_retries, vc, reaction_length_, rrec_, last_time_);
      propagator.step_all();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   //TODO
   void check()
   {

   }

   // --------------------------------------------------------------------------------------------------------------------------------

   const std::vector<ParticleID>& get_particles() const { return particles_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   const Shell& get_shell_for_particle(const ParticleID pid) const
   {
      auto sid = shell_map_.find(pid);
      THROW_UNLESS_MSG(not_found, sid != shell_map_.end(), "Particle has no shell!");
      return (*shell_container_.find(sid->second)).second;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   double get_max_rate(const species_map& species) const
   {
      double k_max = 0.0;

      // TODO Rates for particle-surface interactions
      // since surface rates are not divided by 2 to compensate for double reaction attempts.
      k_max *= 2.0;


      // mono-molecular reactions: (NO CHECKED in original)
      //for (auto sid : species)
      //{
      //   const auto& rules = reactions_.query_reaction_rules(sid);
      //   k_max = std::accumulate(rules.begin(), rules.end(), k_max, [](double s, const ReactionRule& rule) { return std::max(s, rule.getK()); });
      //}

      // bi-molecular reactions:
      for (auto sid_a : species) {
          for (auto sid_b : species) {
              auto scale = 1.0;

              const auto &s0 = world_.get_species(sid_a);
              const auto &s1 = world_.get_species(sid_b);
              double r01 = s0.radius() + s1.radius();
              StructureTypeID sdef = world_.get_def_structure_type_id();

//            THROW_UNLESS_MSG(not_implemented, s0.structure_type_id() == sdef && s1.structure_type_id() == sdef, "Structures not yet supported!");

              if(s0.structure_type_id() == sdef || s1.structure_type_id() == sdef)
              {
                  // for Cuboidal Surface (world)
                  scale = 1.0 / (4 * M_PI * r01 * r01); // 4 pi r^2 = Surface of a sphere
              }
              else
              {
                  // We assume PlanarSurface here
                  Log("GFRD").info() << "Assuming structure that a BD-simulated particle is attached to is a PlanarSurface. Other types of structures not yet supported.";
                  scale = 1.0 / (2 * M_PI * r01); // 2 pi r = circumference of a circle
              }

              const auto &rules = reactions_.query_reaction_rules(sid_a, sid_b);
              k_max = std::accumulate(rules.begin(), rules.end(), k_max, [scale](double s, const ReactionRule &rule) {
                  return std::max(s, scale * rule.getK());
              });
          }
      }

      return k_max;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   species_map get_species() const
   {
      species_map species;
      for (auto pid : particles_)
      {
         const auto& pid_pair = world_.get_particle(pid);
         species.insert(pid_pair.second.sid());
      }
      return species;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double get_min_radius(const species_map& species) const
   {
      double radius_min(std::numeric_limits<double>::max());
      for (auto sid : species)
      {
         const auto& sp = world_.get_species(sid);
         auto radius = sp.radius();
         if (radius < radius_min) radius_min = radius;
      }
      return radius_min;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Determines the largest diffusion constant and largest drift coefficient inside the multi-particle container. Finds whether the
   // movement is dominated by diffusion or drift for a given length scale. Returns the typical time to travel that length scale based
   // on which type of movement is dominant.
   double get_min_tau_Dv(double r_typical, const species_map& species) const
   {
      ASSERT(r_typical > 0.0);

      double tau_dominant, tau_min = -1.0;
      double LARGE_PREF = std::numeric_limits<double>::infinity();

      for (auto sid : species)
      {
         const auto& sp = world_.get_species(sid);

         if (sp.D() == 0.0 && sp.v() == 0.0)
            tau_dominant = std::numeric_limits<double>::infinity();
         else
         {
            auto tau_D = sp.D() > 0.0 ? 2.0 * std::pow(r_typical,2) / sp.D() : LARGE_PREF * r_typical / sp.v();
            auto tau_v = sp.v() > 0.0 ? r_typical / sp.v() : LARGE_PREF * 2.0 * std::pow(r_typical,2) / sp.D();
            tau_dominant = tau_v < 0.1 * tau_D ? tau_v : tau_D;
         }

         if (tau_min < 0.0 || tau_dominant < tau_min)
            tau_min = tau_dominant;
      }

      assert(tau_min > 0.0);
      return tau_min;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   class volume_clearer final : public VolumeClearer
   {
   public:
      volume_clearer(Multi& multi) : multi_(multi) { }

      bool check_move(const Sphere& particle, const ParticleID pid) const override
      {
         const auto& shell = multi_.get_shell_for_particle(pid);

         THROW_UNLESS_MSG(illegal_state, shell.code() == Shell::Code::MULTI, "Expected only Multi shells in multi shell container!");
         THROW_UNLESS_MSG(illegal_state, shell.shape() == Shell::Shape::SPHERE, "Expected only Spherical shells in multi shell container!");

         double distance = shell.get_sphere().distance(particle.position());

         if (distance > -(particle.radius() + multi_.reaction_length_))
         {
            // stop multi defusing, change to escape
            if (multi_.eventType_ == EventType::MULTI_DIFFUSION)
               multi_.eventType_ = EventType::MULTI_ESCAPE;

            // TODO burst (need EGFD sim)
            // burst_volume(shape.position, shape.radius, ignore = [self.outer_.domain_id, ]);
         }
         return true;      // we allow all moves, particle overlap already checked (we cant, because of singular ignore, not suited for bimolecularreaction)
      }


      void remove_particle(const ParticleID pid) const
      {
         auto sid = multi_.shell_map_.find(pid);
         if (sid != multi_.shell_map_.end())
         {
            // Particles that are created in the same Multi.Step() do not have a shell (yet)!
            multi_.shell_container_.erase(sid->second);
            multi_.shell_map_.erase(pid);
         }

         auto iter = std::find(multi_.particles_.begin(), multi_.particles_.end(), pid);
         if (multi_.particles_.cend() != iter) { std::iter_swap(iter, multi_.particles_.rbegin()); multi_.particles_.pop_back(); }   // swap and pop, faster then erase
      }

      void unimolecular_reaction(const ParticleID pid, const std::vector<ParticleID>& products) const override
      {
         remove_particle(pid);

         if (multi_.eventType_ == EventType::MULTI_DIFFUSION)
            multi_.eventType_ = EventType::MULTI_UNIMOLECULAR_REACTION;

         // add new particles to the Multi
         std::copy(products.begin(), products.end(), std::back_inserter(multi_.particles_));
      }

      void bimolecular_reaction(const ParticleID pid1, const ParticleID pid2, const std::vector<ParticleID>& products) const override
      {
         remove_particle(pid1);
         remove_particle(pid2);

         if (multi_.eventType_ == EventType::MULTI_DIFFUSION)
            multi_.eventType_ = EventType::MULTI_BIMOLECULAR_REACTION;

         // add new particles to the Multi
         std::copy(products.begin(), products.end(), std::back_inserter(multi_.particles_));
      }

   private:
      Multi& multi_;
   };

   // --------------------------------------------------------------------------------------------------------------------------------

protected:
   friend class Persistence;

   World& world_;
   ReactionRuleCollection& reactions_;
   RandomNumberGenerator& rng_;
   reaction_recorder* rrec_;
   std::vector<ParticleID> particles_;
   std::unordered_map<ParticleID, ShellID> shell_map_;
   shell_matrix_type& shell_container_;
   double reaction_length_, start_time_;
};

// --------------------------------------------------------------------------------------------------------------------------------
