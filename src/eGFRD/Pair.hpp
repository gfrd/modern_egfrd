#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <utility>
#include "Domain.hpp"
#include "Particle.hpp"
#include "ParticleID.hpp"
#include "Shell.hpp"
#include "ShellID.hpp"
#include "ReactionRuleCollection.hpp"
#include "MatrixSpace.hpp"
#include "World.hpp"
#include "ShellCreateUtils.hpp"
#include "randomNumberGenerator.hpp"
#include <GreensFunction.hpp>
#include "GreenFunctionHelpers.hpp"
#include "Matrix4.hpp"
#include <GreensFunction3DRadInf.hpp>
#include <GreensFunction3DAbs.hpp>
#include <GreensFunction3D.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

class Pair : public Domain
{
public:
   using particle_id_pair = Particle::particle_id_pair;
   using shell_id_pair = Shell::shell_id_pair;
   using ReactionRules = ReactionRuleCollection::reaction_rule_set;
   using shell_matrix_type = MatrixSpace<Shell, ShellID, World::particle_matrix_type::SizeX, World::particle_matrix_type::SizeY, World::particle_matrix_type::SizeZ>;
   using shell_distance_checker = ShellCreateUtils::shell_distance_check<shell_matrix_type>;
   using position_structid_pair = World::position_structid_pair;

   // --------------------------------------------------------------------------------------------------------------------------------

   static const double CUTOFF_FACTOR;        // CUTOFF_FACTOR is a threshold to choose between the real and approximate Green's functions.

   // --------------------------------------------------------------------------------------------------------------------------------

   explicit Pair(const DomainID did, const particle_id_pair& pid_pair1, const particle_id_pair& pid_pair2, const shell_id_pair& sid_pair, const ReactionRules& reactions) noexcept : Domain(did),
      pid_pair1_(pid_pair1.second.D() <= pid_pair2.second.D() ? pid_pair1 : pid_pair2),       // store by smallest defusing coefficient
      pid_pair2_(pid_pair1.second.D() <= pid_pair2.second.D() ? pid_pair2 : pid_pair1),
      sid_pair_(sid_pair), rrules_(reactions), iv_(pid_pair2_.second.position() - pid_pair1_.second.position()),
      a_R_(0), a_r_(0), single1_ktotal_(0), single2_ktotal_(0), pid_reactingsingle_(0), rrule_(nullptr), gf_com_(), gf_iv_() { }

   // --------------------------------------------------------------------------------------------------------------------------------

   const particle_id_pair& pip1() const { return pid_pair1_; }            // this returns a value copy of the Particle, and it should! (because Domain is about to be erased)
   ParticleID particle1_id() const { return pid_pair1_.first; }
   const Particle& particle1() const { return pid_pair1_.second; }

   const particle_id_pair& pip2() const { return pid_pair2_; }
   ParticleID particle2_id() const { return pid_pair2_.first; }
   const Particle& particle2() const { return pid_pair2_.second; }

   const shell_id_pair& shell_pair() const { return sid_pair_; }        // this returns a value copy of the Shell !
   ShellID shell_id() const { return sid_pair_.first; }
   const Shell& shell() const { return sid_pair_.second; }

   // --------------------------------------------------------------------------------------------------------------------------------

   size_t num_shells() const override { return 1; }

   Multiplicity multiplicity() const override { return Multiplicity::PAIR; }

   virtual shell_id_pair_ref get_shell() const override { return std::make_pair(sid_pair_.first, std::ref<const Shell>(sid_pair_.second)); }

   virtual std::vector<shell_id_pair_ref> get_shell_list() const override { return std::vector<shell_id_pair_ref >({ get_shell() }); }

   // --------------------------------------------------------------------------------------------------------------------------------

   const Vector3& iv() const { return iv_; }

   double r0() const { return iv_.length(); }

   double a_R() const { return a_R_; }

   double a_r() const { return a_r_; }

   double sigma() const { return pid_pair1_.second.radius() + pid_pair2_.second.radius(); }

   double D_tot() const { return pid_pair1_.second.D() + pid_pair2_.second.D(); }

   double D_geom() const { return std::sqrt(pid_pair1_.second.D() * pid_pair2_.second.D()); }

   double D_R() const { return pid_pair1_.second.D() * pid_pair2_.second.D() / D_tot(); }

   ParticleID get_reacting_single() const noexcept { return pid_reactingsingle_; }

   void set_k_totals(SpeciesTypeID sid1, double ktot1, double ktot2)
   {
      if (sid1 == pid_pair1_.second.sid())
      {
         single1_ktotal_ = ktot1;
         single2_ktotal_ = ktot2;
      }
      else
      {
         single1_ktotal_ = ktot2;
         single2_ktotal_ = ktot1;
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual bool create_updated_shell(const shell_matrix_type& smat, const World& world, ShellID sid1, ShellID sid2) = 0;

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string as_string() const override {
      return make_string() << ", " << pid_pair1_.first << ", " << pid_pair1_.second << ", " <<
         ", " << pid_pair2_.first << ", " << pid_pair2_.second << ", " << sid_pair_.first << ", " << sid_pair_.second <<
         ", iv=" << iv_ << /*", reactions=" << rrules_ << */ ", a_r=" << a_r_ << ", a_R" << a_R_;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double get_min_pair_size() const
   {
      ASSERT(iv_.length());

      double dist_from_com1 = r0() * (D_tot() > 0 ? pid_pair1_.second.D() / D_tot() : 0.5);
      double dist_from_com2 = r0() * (D_tot() > 0 ? pid_pair2_.second.D() / D_tot() : 0.5);

      // Calculate total radii including the margin for the burst volume for the particles
      return std::max(dist_from_com1 + pid_pair1_.second.radius() * GfrdCfg.SINGLE_SHELL_FACTOR, dist_from_com2 + pid_pair2_.second.radius() * GfrdCfg.SINGLE_SHELL_FACTOR);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double k_total() const        // inter-particle
   {
      // calculates the total rate for a list of reaction rules
      // The probability for the reaction to happen is proportional to the sum of the rates of all the possible reaction types.
      double k_tot = 0;
      for (auto& rr : rrules_)
         k_tot += rr.getK();
      return k_tot;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void determine_next_event(RandomNumberGenerator& rng)
   {
      double dt_com = GreenFunctionHelper::draw_time_wrapper(rng, *gf_com_.get());                 // EventType.COM_ESCAPE
      double dt_iv = GreenFunctionHelper::draw_time_wrapper(rng, *gf_iv_.get());                   // EventType.IV_EVENT
      double dt_reaction1 = GreenFunctionHelper::draw_reaction_time(single1_ktotal_, rng);         // EventType.SINGLE_REACTION
      double dt_reaction2 = GreenFunctionHelper::draw_reaction_time(single2_ktotal_, rng);         // EventType.SINGLE_REACTION

      // get minimal time step
      dt_ = std::min(std::min(dt_com, dt_iv), std::min(dt_reaction1, dt_reaction2));
      if (dt_ == dt_com)
         eventType_ = EventType::COM_ESCAPE;
      else if (dt_ == dt_iv)
         eventType_ = EventType::IV_EVENT;
      else if (dt_ == dt_reaction1)
      {
         eventType_ = EventType::SINGLE_REACTION;
         pid_reactingsingle_ = pid_pair1_.first;
      }
      else if (dt_ == dt_reaction2)
      {
         eventType_ = EventType::SINGLE_REACTION;
         pid_reactingsingle_ = pid_pair2_.first;
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual const PairGreensFunction& choose_pair_greens_function() const
   {
      return *gf_iv_.get();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   Vector3 draw_new_com(RandomNumberGenerator& rng) const
   {
      // draws a new coordinate for the CoM in world coordinates
      double r;
      if (eventType_ == EventType::IV_ESCAPE)
         r = a_R_;
      else
         r = GreenFunctionHelper::draw_r_wrapper(rng, *gf_com_.get(), dt_, a_R_);

      // Add displacement to old CoM. This assumes(correctly) that r0 = 0 for the CoM. 
      // Compare this to 1D singles, where r0 is not necessarily 0.

       // note that we need to make sure that the com and com_vector are in the structure to prevent
       // the particle leaving the structure due to numerical errors
      return com_ + create_com_vector(r, rng);
   }

   Vector3 draw_new_iv(RandomNumberGenerator& rng) const
   {
      const auto& gf = choose_pair_greens_function();

      double r;
      if (eventType_ == EventType::IV_ESCAPE)
         r = a_r_;
      else if (eventType_ == EventType::IV_REACTION)
         r = sigma();      // maybe this should be zero(product particle is at CoM)
      else
         r = GreenFunctionHelper::draw_r_wrapper(rng, gf, dt_, a_r_, sigma());

      // note that we need to make sure that the interparticle vector is in the structure to prevent
      // the particle leaving the structure due to numerical errors
      return create_interparticle_vector(gf, r, rng);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   const ReactionRule& draw_reaction_rule(RandomNumberGenerator&rng) const
   {
      // draws a reaction rules out of a list of reaction rules based on their relative rates
      if (rrule_ != nullptr) return *rrule_;

      double rnd = rng.uniform(0, k_total());
      ReactionRules::iterator i = rrules_.begin();
      double k_sum = 0.0;
      while (k_sum + i->getK() < rnd) { k_sum += i->getK();  ++i; }

      rrule_ = &(*i);
      return *rrule_;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   const ReactionRule& draw_single_reaction_rule(RandomNumberGenerator&rng, const ReactionRuleCollection& rules) const
   {
      // for a SINGLE_RECTION event, we need a ReactionRule, draw one for the correct particle/species_type
      ASSERT(eventType_ == EventType::SINGLE_REACTION);
      ASSERT(pid_reactingsingle_);

      SpeciesTypeID sid = (pid_reactingsingle_ == pid_pair1_.first) ? pid_pair1_.second.sid() : pid_pair2_.second.sid();
      const auto& rr = rules.query_reaction_rules(sid);
      double k_total = 0;
      for (auto r : rr) k_total += r.getK();
      ASSERT(k_total > 0);    // the particle cannot SINGLE_REACT, because it has no rules!

      double rnd = rng.uniform(0, k_total);
      ReactionRules::iterator i = rr.begin();
      double k_sum = 0.0;
      while (k_sum + i->getK() < rnd) { k_sum += i->getK();  ++i; }
      return *i;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual Vector3 create_com_vector(double r, RandomNumberGenerator& rng) const = 0;
   virtual Vector3 create_interparticle_vector(const PairGreensFunction& gf, double r, RandomNumberGenerator& rng) const = 0;
   virtual std::pair<position_structid_pair, position_structid_pair> do_back_transform(Vector3 com, Vector3 iv, double D1, double D2, double r1, double r2, StructureID s1, StructureID s2, Vector3 unit, const World& world) const = 0;
   virtual EventType draw_iv_event_type(RandomNumberGenerator& rng) = 0;

   // --------------------------------------------------------------------------------------------------------------------------------

   std::pair<position_structid_pair, position_structid_pair> draw_new_position(RandomNumberGenerator& rng, const World& world) const
   {
      Vector3 new_com = draw_new_com(rng);
      Vector3 new_iv = draw_new_iv(rng);
      Vector3 unit_z = sid_pair_.second.shape() == Shell::Shape::SPHERE ? Vector3::uz : sid_pair_.second.get_cylinder().unit_z();

      return do_back_transform(new_com, new_iv, pid_pair1_.second.D(), pid_pair2_.second.D(),
         pid_pair1_.second.radius(), pid_pair2_.second.radius(),
         pid_pair1_.second.structure_id(), pid_pair2_.second.structure_id(), unit_z, world);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:
   friend class Persistence;

   const particle_id_pair pid_pair1_;
   const particle_id_pair pid_pair2_;
   shell_id_pair sid_pair_;
   const ReactionRules rrules_;
   Vector3 iv_, com_;
   double a_R_, a_r_;
   double single1_ktotal_, single2_ktotal_;
   ParticleID pid_reactingsingle_;
   mutable const ReactionRule *rrule_;

   std::unique_ptr<GreensFunction> gf_com_;
   std::unique_ptr<PairGreensFunction> gf_iv_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class PairSpherical : public Pair
{
public:
   PairSpherical(const DomainID did, const particle_id_pair& pid_pair1, const particle_id_pair& pid_pair2, const shell_id_pair& sid_pair, const ReactionRules& reactions)
      : Pair(did, pid_pair1, pid_pair2, sid_pair, reactions)
   {
      THROW_UNLESS(illegal_state, sid_pair.second.shape() == Shell::Shape::SPHERE);
   }

   const char* type_name() const override { return "PairSpherical"; }

   GFRD_EXPORT void determine_radii(double r0, double shellsize);

   GFRD_EXPORT bool create_updated_shell(const shell_matrix_type& smat, const World& world, ShellID sid1, ShellID sid2) override;

   Vector3 create_com_vector(double r, RandomNumberGenerator& rng) const override
   {
      return r * Vector3::random(rng);
   }

   GFRD_EXPORT const PairGreensFunction& choose_pair_greens_function() const override;

   EventType draw_iv_event_type(RandomNumberGenerator& rng) override
   {
      auto event = GreenFunctionHelper::draw_eventtype_wrapper(rng, *gf_iv_, dt_);
      if (event == GreensFunction::EventKind::IV_REACTION)
         eventType_ = EventType::IV_REACTION;
      else if (event == GreensFunction::EventKind::IV_ESCAPE)
         eventType_ = EventType::IV_ESCAPE;
      else THROW_EXCEPTION(illegal_state, "Unknown event type.");
      return eventType_;
   }

   Vector3 create_interparticle_vector(const PairGreensFunction& gf, double r, RandomNumberGenerator& rng) const override
   {
      Vector3 new_iv;
      double theta = GreenFunctionHelper::draw_theta_wrapper(rng, gf, r, dt_);
      if (theta == 0.0)
         new_iv = iv_;                 // no rotation is necessary->new_iv = new_iv
      else if (std::fmod(theta, M_PI) != 0.0)
      {
         // alternative calculation, rotate the old_iv to the new theta
         Vector3 rotation_axis = Vector3::cross(iv_, Vector3::uz);
         new_iv = Vector3::transformVector(iv_, Matrix4::createRotationA(theta, rotation_axis.normal()));
         // rotate the new_iv around the old_iv with angle phi
         double phi = rng.uniform(0, 1) * 2 * M_PI;
         new_iv = Vector3::transformVector(new_iv, Matrix4::createRotationA(phi, iv_.normal()));
      }
      else
         new_iv = -iv_;                //  theta == pi->just mirror the old_iv

      return (r / r0()) * new_iv;      // adjust length of the vector, note that r0 = length (old_iv)
   }


   virtual std::pair<position_structid_pair, position_structid_pair> do_back_transform(Vector3 com, Vector3 iv, double D1, double D2, double, double, StructureID s1, StructureID s2, Vector3, const World&) const override
   {
      // Here we assume that the com and iv are really in the structure and no adjustments have to be made

      // Since this is meant to be a general class, we do not check explicitly whether the structures
      // are really the same, but drop a warning if they are not (A possible situation where this matters
      // is when one particle is on a substructure of the other particle's structure; then the structure IDs
      // are different, but the calculations then still work perfectly fine).
      THROW_UNLESS(illegal_state, s1 == s2);      // log.warn('Particles live on different structures in StandardPair, structure1=%s, structure2=%s')

      Vector3 pos1, pos2;
      double D_tot = D1 + D2;
      if (D_tot != 0)
      {
         pos1 = com - iv * (D1 / D_tot);
         pos2 = com + iv * (D2 / D_tot);
      }
      else
      {
         pos1 = com - iv * 0.5;
         pos2 = com + iv * 0.5;
      }
      return std::make_pair<position_structid_pair, position_structid_pair>(position_structid_pair(pos1, s1), position_structid_pair(pos2, s2));
   }

private:
   friend class Persistence;

   mutable std::unique_ptr<PairGreensFunction> gf_tmp_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class PairCylindrical : public Pair
{
public:

   PairCylindrical(const DomainID did, const particle_id_pair& pid_pair1, const particle_id_pair& pid_pair2, const shell_id_pair& sid_pair, const ReactionRules& reactions)
      : Pair(did, pid_pair1, pid_pair2, sid_pair, reactions)
   {
      THROW_UNLESS(illegal_state, sid_pair.second.shape() == Shell::Shape::CYLINDER);
   }

   const char* type_name() const override { return "PairCylindrical"; }
};

// --------------------------------------------------------------------------------------------------------------------------------
