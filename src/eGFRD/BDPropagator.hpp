#pragma once

#include "randomNumberGenerator.hpp"
#include "ReactionRuleCollection.hpp"
#include "ReactionRecorder.hpp"
#include "World.hpp"
#include <numeric>
#include <random>

// --------------------------------------------------------------------------------------------------------------------------------

class VolumeClearer
{
public:
   virtual ~VolumeClearer() = default;

   virtual bool check_move(const Sphere&, const ParticleID) const { return true; }
   virtual void unimolecular_reaction(const ParticleID, const std::vector<ParticleID>&) const { }
   virtual void bimolecular_reaction(const ParticleID, const ParticleID, const std::vector<ParticleID>&) const { }
};

// --------------------------------------------------------------------------------------------------------------------------------


class BDPropagator
{
public:
   using particle_id_pair = Particle::particle_id_pair;
   using particle_id_pair_and_distance = std::pair<particle_id_pair, double>;
   using particle_id_pair_and_distance_list = std::vector<particle_id_pair_and_distance>;
   using pos_structid_pair = std::pair<Vector3, StructureID>;
   using pos_sid_pair_pair = std::pair<pos_structid_pair, pos_structid_pair>;

   // --------------------------------------------------------------------------------------------------------------------------------

   enum class ReactionResult { Ok, None, Error };

   // --------------------------------------------------------------------------------------------------------------------------------

   BDPropagator() = delete;

   explicit BDPropagator(World& world, const std::vector<ParticleID>& particles, ReactionRuleCollection& rules, RandomNumberGenerator& rng, double dt, int max_retry_count, const VolumeClearer& vc, double rl, reaction_recorder* rrec, double time) noexcept
      : world_(world), reaction_rules_(rules), rng_(rng), dt_(dt), max_retry_count_(max_retry_count), reaction_length_(rl), queue_(), vc_(vc), rejected_move_count_(0), rrec_(rrec), time_(time)
   {
      queue_.reserve(particles.size());
      std::copy(particles.begin(), particles.end(), std::back_inserter(queue_));

      // shuffle the queue
      std::mt19937 urng(rng_.uniform_int(0, INT32_MAX-1));
      std::shuffle(queue_.begin(), queue_.end(), urng);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   size_t get_rejected_move_count() const { return rejected_move_count_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void step_all()
   {
      while (!queue_.empty())
      {
         auto pid = queue_.back();
         queue_.pop_back();
         step_particle(pid);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   ReactionResult attempt_unimolecular_reaction(const particle_id_pair& pip) const
   {
      const auto& rules = reaction_rules_.query_reaction_rules(pip.second.sid());
      if (rules.size() == 0) return ReactionResult::None; // no reaction rules.

      double rnd = rng_() / dt_;
      double k_cumm = 0.0;

      for (const auto& rr : rules)
      {
         k_cumm += rr.getK();
         if (k_cumm > rnd)
         {
            const auto& p = rr.get_products();
            switch (p.size())
            {
            case 0:
            {
               world_.remove_particle(pip.first);
               vc_.unimolecular_reaction(pip.first, std::vector<ParticleID>());
               if (rrec_) rrec_->StoreDecayReaction(time_, rr.id(), pip.first);
               break;
            }
            case 1:
            {
               if (!attempt_unimolecular_one_prodcut_reaction(pip, p[0], rr.id())) return ReactionResult::Error;
               break;
            }
            case 2:
            {
               if (!attempt_unimolecular_two_product_reaction(pip, p[0], p[1], rr.id())) return ReactionResult::Error;
               break;
            }
            default: throw not_implemented("more than two products are currently not supported!");
            }
            return ReactionResult::Ok;
         }
      }
      return ReactionResult::None;     // no reactions
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool attempt_unimolecular_one_prodcut_reaction(const particle_id_pair& pip, const SpeciesTypeID sid1, const ReactionRuleID rrid) const
   {
      auto reactant_position = pip.second.position();
      auto reactant_sid = pip.second.structure_id();

      auto product_species = world_.get_species(sid1);
      auto product_shape = Sphere(reactant_position, product_species.radius());

      // check for overlapping particles
      if (world_.test_particle_overlap(product_shape, pip.first))
         return false; // no space due to other particle.

      if (!vc_.check_move(product_shape, pip.first))
         return false;

      // no checking of movement to other structure types!!!
      // no checking whether particle has overlapping surface!!!

      world_.remove_particle(pip.first);
      auto pip_new = world_.new_particle(Particle(product_species, reactant_sid, reactant_position));

      vc_.unimolecular_reaction(pip.first, std::vector<ParticleID>{pip_new.first});
      
      if (rrec_) rrec_->StoreUniMolecularReaction(time_, rrid, pip.first,pip_new.first);
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool attempt_unimolecular_two_product_reaction(const particle_id_pair& pip, const SpeciesTypeID sid1, const SpeciesTypeID sid2, const ReactionRuleID rrid) const
   {
      auto reactant_pos = pip.second.position();
      auto reactant_structure = world_.get_structure(pip.second.structure_id());

      auto species1 = world_.get_species(sid1);
      auto species2 = world_.get_species(sid2);

      int retry_count = max_retry_count_;
      while (retry_count--)
      {
         pos_sid_pair_pair posAB_pair = reactant_structure->get_pos_sid_pair_pair(*reactant_structure, reactant_pos, species1, species2, reaction_length_, rng_);

         auto p1pip = world_.apply_boundary(posAB_pair.first);
         auto p2pip = world_.apply_boundary(posAB_pair.second);

         // check for particle overlap with product 1.
         Sphere new_shape1(p1pip.first, species1.radius());
         if (world_.test_particle_overlap(new_shape1, pip.first)) continue;

         // check for particle overlap with product 2.
         Sphere new_shape2(p2pip.first, species2.radius());
         if (world_.test_particle_overlap(new_shape2, pip.first)) continue;

         posAB_pair = make_pair_move(species1, species2, p1pip, p2pip, pip.first);

         if (!vc_.check_move(new_shape1, pip.first)) continue;
         if (!vc_.check_move(new_shape2, pip.first)) continue;

         // no checking of movement to other structure types!!!
         // no checking whether particle has overlapping surface!!!

         world_.remove_particle(pip.first);
         auto pip1_new = world_.new_particle(Particle(species1, posAB_pair.first.second, posAB_pair.first.first));
         auto pip2_new = world_.new_particle(Particle(species2, posAB_pair.second.second, posAB_pair.second.first));

         vc_.unimolecular_reaction(pip.first, std::vector<ParticleID>{ pip1_new.first, pip2_new.first});

         if (rrec_) rrec_->StoreUnbindingReaction(time_, rrid, pip.first, pip1_new.first, pip2_new.first);
         return true;
      }

      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   ReactionResult attempt_bimolecular_reaction(const particle_id_pair& pip1, const particle_id_pair& pip2)
   {
      const auto& rules = reaction_rules_.query_reaction_rules(pip1.second.sid(), pip2.second.sid());
      if (rules.size() == 0) return ReactionResult::None;

      double k_tot = k_total(pip1.second.sid(), pip2.second.sid());
      double rnd = k_tot * rng_();
      double k_cumm = 0;

      for (const auto& rr : rules)
      {
         k_cumm += rr.getK();
         if (k_cumm >= rnd)
         {
            const auto& p = rr.get_products();
            switch (p.size())
            {
            case 0:
               {
                  world_.remove_particle(pip1.first); 
                  world_.remove_particle(pip2.first); 
                  vc_.bimolecular_reaction(pip1.first, pip2.first, std::vector<ParticleID>()); 
                  if (rrec_) rrec_->StoreAnnihilationReaction(time_, rr.id(), pip1.first, pip2.first);
                  break;
               }
            case 1:
               {
                  if (!attempt_bimolecular_one_product_reaction(pip1, pip2, p[0], rr.id())) return ReactionResult::Error; 
                  break;
               }
            default: throw not_implemented("bimolecular reactions that produce more than one product are not supported");
            }

            dequeue_particle(pip2.first);
            return ReactionResult::Ok;
         }
      }
      return ReactionResult::None;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool attempt_bimolecular_one_product_reaction(const particle_id_pair& pip1, const particle_id_pair& pip2, const SpeciesTypeID product, const ReactionRuleID rrid) const
   {
      auto reactant1_sid = pip1.second.structure_id();
      auto reactant1_structure = world_.get_structure(reactant1_sid);
      auto reactant1_pos = pip1.second.position();

      auto reactant2_sid = pip2.second.structure_id();
      auto reactant2_structure = world_.get_structure(reactant2_sid);
      //auto reactant2_pos = pip2.second.position();

      auto product_sid = reactant1_sid;
      auto product_structure = world_.get_structure(product_sid);
      auto product_pos = reactant1_pos;
      auto product_species = world_.get_species(product);

      const Sphere new_shape(product_pos, product_species.radius());
      if (world_.test_particle_overlap(new_shape, pip1.first, pip2.first))
         return false; // no space due to particle.

      if (!vc_.check_move(new_shape, pip1.first))
         return false;  // no space due to other particle.

      world_.remove_particle(pip1.first);
      auto pip_new = world_.new_particle(Particle(product_species, product_sid, product_pos));
      world_.remove_particle(pip2.first);

      vc_.bimolecular_reaction(pip1.first, pip2.first, std::vector<ParticleID>{ pip_new.first });

      if (rrec_) rrec_->StoreBindingReaction(time_, rrid, pip1.first, pip2.first, pip_new.first);
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Moves structure if not overlapping with particles or surfaces. */
   pos_structid_pair make_move(const SpeciesType& species, const pos_structid_pair& psp, const ParticleID ignore) const
   {
      // no movement when diffusion constant is zero.
      if (species.D() == 0) return psp;

      auto structure = world_.get_structure(StructureID(psp.second));
      auto displacement = structure->bd_displacement(species.v() * dt_, std::sqrt(2.0 * species.D() * dt_), rng_);

      auto pos_structid = world_.apply_boundary(std::make_pair(psp.first + displacement, psp.second));
      auto shape = Sphere(pos_structid.first, species.radius());

      // return old pos if overlapping with particles.
      if (world_.test_particle_overlap(shape, ignore)) return psp;

      // return old pos if overlapping with surfaces.
      if (world_.test_surface_overlap(shape, psp.first, pos_structid.second, species.radius())) return psp;

      return pos_structid;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   pos_sid_pair_pair make_pair_move(const SpeciesType& s0, const SpeciesType& s1, const pos_structid_pair& psp0, const pos_structid_pair& psp1, const ParticleID ignore) const
   {
      double min_distance = s0.radius() + s1.radius();

      int  move_count = 1 + rng_.uniform_int(0, 1);
      bool move_first = rng_.uniform_int(0, 1) == 1;

      pos_structid_pair new_psp0 = psp0;
      pos_structid_pair new_psp1 = psp1;

      while (move_count--)
      {
         if (move_first)
         {
            new_psp0 = make_move(s0, psp0, ignore);
            if (world_.distance(new_psp0.first, new_psp1.first) < min_distance)
               new_psp0 = psp0;
         }
         else
         {
            new_psp1 = make_move(s1, psp1, ignore);
            if (world_.distance(new_psp0.first, new_psp1.first) < min_distance)
               new_psp1 = psp1;
         }
         move_first = !move_first;
      }
      return std::make_pair(new_psp0, new_psp1);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void step_particle(ParticleID pid)
   {
      auto pid_pair = world_.get_particle(pid);
      auto result = attempt_unimolecular_reaction(pid_pair);
      if (result == ReactionResult::Ok) return;
      if (result == ReactionResult::Error) { ++rejected_move_count_; return; }

      auto species = world_.get_species(pid_pair.second.sid());
      double radius = species.radius();
      auto old_pos = pid_pair.second.position();
      auto old_struct_id = pid_pair.second.structure_id();
      auto structure = world_.get_structure(old_struct_id);

      auto new_pos = old_pos;
      auto new_structure_id = old_struct_id;
      if (species.D() != 0.0 || species.v() != 0.0)
      {
         auto displacement = structure->bd_displacement(species.v() * dt_, std::sqrt(2.0 * species.D() * dt_), rng_);
         auto pos_structid = world_.apply_boundary(std::make_pair(old_pos + displacement, old_struct_id));
         new_pos = pos_structid.first;
         new_structure_id = pos_structid.second;
      }

      auto overlap = world_.check_particle_overlap(Sphere(new_pos, radius + reaction_length_), pid);
      bool bounced = std::any_of(overlap.begin(), overlap.end(), [&](const particle_id_pair_and_distance& e) {return e.second < -reaction_length_; });
      bool bounced_struct = world_.test_surfaces_overlap(Sphere(new_pos, radius), old_pos, 0,
              std::vector<StructureID>({world_.get_def_structure_id(), new_structure_id}));

       auto structures = world_.get_structures();
       for(const auto& structure : structures) {
           if (structure.get()->id() == world_.get_def_structure_id() || structure.get()->id() == new_structure_id) {
               continue;
           }

           auto distance = structure.get()->distance(world_.cyclic_transpose(old_pos, structure.get()->position()));
           if(fabs(distance) < radius) {
               // TODO: this scenario should never happen, but currently does. A root cause fix is necessary later in time.
               Log("GFRD").warn() << "Particle " << pid << " overlaps structure, moving it forcefully to structure boundary";

               auto plane = dynamic_cast<PlanarSurface*>(structure.get());
               if(plane == nullptr) {
                   THROW_EXCEPTION(unsupported, "Particle overlaps a structure that is not a PlanarSurface. This is not yet supported.");
               }

               // Move particle back to the side of the
               auto projected = plane->project_point(old_pos);
               auto z_distance = projected.second.first / GfrdCfg.SAFETY;
               auto to_move = z_distance >= 0 ? radius - z_distance : -radius - z_distance;
               new_pos = world_.apply_boundary(old_pos + plane->shape().unit_z() * to_move);

               // Force particle move, which would not happen if we say a bounce occurred
               bounced = false;
               bounced_struct = false;

               break;
           }
       }

      if(bounced_struct)
      {
          Log("GFRD").info() << "Particle " << pid << " bounced off of structure";
      }

      if (bounced || bounced_struct)
      {
         // restore old position and structure_id
         new_pos = old_pos;
         new_structure_id = old_struct_id;

         // re-get the reaction partners (particles), now on old position.
         auto overlap_nomove = world_.check_particle_overlap(Sphere(old_pos, radius + reaction_length_), pid);
         overlap.swap(overlap_nomove);
      }

      double rnd = rng_();
      double accumulated_prob = 0.0;
      for (auto& overlap_particle : overlap)
      {
         auto s0 = species;
         auto s1 = world_.get_species(overlap_particle.first.second.sid());
         auto s1_struct = world_.get_structure(overlap_particle.first.second.structure_id());
         double prob_increase = k_total(s0.id(), s1.id()) * dt_ / (2. * s1_struct->particle_reaction_volume(s0.radius() + s1.radius(), reaction_length_));
         accumulated_prob += prob_increase;

         if (accumulated_prob > rnd)
         {
            result = attempt_bimolecular_reaction(pid_pair, overlap_particle.first);
            if (result == ReactionResult::Ok) return;
            if (result == ReactionResult::Error) { ++rejected_move_count_; accumulated_prob -= prob_increase; }
         }
      }

      if(!bounced)
      {
         auto update = std::make_pair(pid, Particle(species.id(), Sphere(new_pos, radius), new_structure_id, species.D(), species.v()));
         if (vc_.check_move(update.second.shape(), pid)) world_.update_particle(update);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double k_total(const SpeciesTypeID sid_a, const SpeciesTypeID sid_b) const
   {
      const auto& rules = reaction_rules_.query_reaction_rules(sid_a, sid_b);
      return std::accumulate(rules.begin(), rules.end(), 0.0, [](double s, const ReactionRule& rule) { return s + rule.getK(); });
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void dequeue_particle(const ParticleID pid)
   {
      auto iter = std::find(queue_.begin(), queue_.end(), pid);
      if (queue_.cend() != iter) { std::iter_swap(iter, queue_.rbegin()); queue_.pop_back(); }   // swap and pop, faster then erase
   }

   // --------------------------------------------------------------------------------------------------------------------------------


protected:
   World& world_;
   const ReactionRuleCollection& reaction_rules_;
   RandomNumberGenerator& rng_;
   const double dt_;
   const int max_retry_count_;
   const double reaction_length_;

   std::vector<ParticleID> queue_;
   const VolumeClearer& vc_;
   size_t rejected_move_count_;
   reaction_recorder* rrec_;
   const double time_;

   // --------------------------------------------------------------------------------------------------------------------------------

};
