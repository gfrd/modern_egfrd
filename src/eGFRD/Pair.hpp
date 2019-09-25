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
#include "scaling.hpp"
#include <GreensFunction3DRadInf.hpp>
#include <GreensFunction3DAbs.hpp>
#include <GreensFunction3D.hpp>
#include <GreensFunction3DRadAbs.hpp>
#include <GreensFunction2DAbsSym.hpp>
#include <GreensFunction2DRadAbs.hpp>

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

    virtual double sigma() const { return pid_pair1_.second.radius() + pid_pair2_.second.radius(); }

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
         r = GreenFunctionHelper::draw_r_wrapper(rng, *gf_com_, dt_, a_R_);

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
      // for a SINGLE_REACTION event, we need a ReactionRule, draw one for the correct particle/species_type
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

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual EventType draw_iv_event_type(RandomNumberGenerator& rng)
   {
       auto event = GreenFunctionHelper::draw_eventtype_wrapper(rng, *gf_iv_, dt_);
       if (event == GreensFunction::EventKind::IV_REACTION)
           eventType_ = EventType::IV_REACTION;
       else if (event == GreensFunction::EventKind::IV_ESCAPE)
           eventType_ = EventType::IV_ESCAPE;
       else THROW_EXCEPTION(illegal_state, "Unknown event type.");
       return eventType_;
   };

   virtual void do_transform(const World& world)
   {
       Vector3 pos1 = pid_pair1_.second.position();
       Vector3 pos2 = pid_pair2_.second.position();
       Vector3 pos2c = world.cyclic_transpose(pos2, pos1);

       Vector3 com = D_tot() > 0 ? (pid_pair2_.second.D() * pos1 + pid_pair1_.second.D() * pos2c) / D_tot() : 0.5 * (pos1 + pos2c);
       com_ = world.apply_boundary(com);

       iv_ = pos2c - pos1;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual std::pair<position_structid_pair, position_structid_pair> do_back_transform(Vector3 com, Vector3 iv, double D1, double D2, double r1, double r2, StructureID s1, StructureID s2, Vector3 unit, const World& world) const
   {
       // Here we assume that the com and iv are really in the structure and no adjustments have to be made

       // Since this is meant to be a general class, we do not check explicitly whether the structures
       // are really the same, but drop a warning if they are not (A possible situation where this matters
       // is when one particle is on a substructure of the other particle's structure; then the structure IDs
       // are different, but the calculations then still work perfectly fine).
       if (s1 != s2)
       {
           Logger::get_logger("EGFRD").warn() << "Particles live on different structures in StandardPair";
       }

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

private:
   friend class Persistence;

   mutable std::unique_ptr<PairGreensFunction> gf_tmp_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class PairPlanar : public Pair
{
public:

   PairPlanar(const DomainID did, const particle_id_pair& pid_pair1, const particle_id_pair& pid_pair2, const std::shared_ptr<Structure>& structure, const shell_id_pair& sid_pair, const ReactionRules& reactions)
      : Pair(did, pid_pair1, pid_pair2, sid_pair, reactions)
   {
//      THROW_UNLESS(illegal_state, sid_pair.second.shape() == Shell::Shape::CYLINDER);

       structure_ = std::dynamic_pointer_cast<PlanarSurface>(structure);
       THROW_UNLESS_MSG(illegal_state, structure != nullptr, "PairPlanar created with structure that is not a PlanarSurface");
   }

   const char* type_name() const override { return "PairPlanar"; }

   void determine_radii(const double shell_size)
   {
    // Determine a_r_ and a_R_ from the size of the protective domain.

    double radius1 = pid_pair1_.second.radius();
    double radius2 = pid_pair2_.second.radius();

    double D1 = pid_pair1_.second.D();
    double D2 = pid_pair2_.second.D();

    ASSERT(r0() >= sigma());        // '%s;  r0 %g < sigma %g' % (self, r0, self.sigma)

    double Da = D1;
    double Db = D2;
    double radiusa = radius1;
    double radiusb = radius2;

    // equalize expected mean t_r and t_R.
    if ((D_geom() - D2) * r0() / D_tot() + shell_size + std::sqrt(D2 / D1) * (radius1 - shell_size) - radius2 < 0)
    {
        Da = D2;
        Db = D1;
        radiusa = radius2;
        radiusb = radius1;
    }

    a_R_ = (D_geom() * (Db * (shell_size - radiusa) + Da * (shell_size - r0() - radiusa))) / (Da * Da + Da * Db + D_geom() * D_tot());
    a_r_ = (D_geom() * r0() + D_tot() * (shell_size - radiusa)) / (Da + D_geom());

    ASSERT(a_R_ + a_r_ * Da / D_tot() + radiusa >= a_R_ + a_r_ * Db / D_tot() + radiusb);
    ASSERT(std::abs(a_R_ + a_r_ * Da / D_tot() + radiusa - shell_size) < 1e-12 * shell_size);                        // here the shell_size is the relevant scale

    ASSERT(a_r_ > 0);
    ASSERT(a_r_ > r0());
    ASSERT(a_R_ > 0 || (feq(a_R_, 0) && (D1 == 0 || D2 == 0)));
   }

    const PairGreensFunction& choose_pair_greens_function() const override
    {
        // Selects between the full solution or an approximation where one of
        // the boundaries is ignored
        auto D_r = D_tot();
        double distance_from_sigma = r0() - sigma();
        double distance_from_shell = a_r_ - r0();
        double threshold_distance = CUTOFF_FACTOR * std::sqrt(4.0 * D_r * dt_);

        // If sigma reachable
        if (distance_from_sigma < threshold_distance)
        {
            // If shell reachable
            if (distance_from_shell < threshold_distance)
            {
                // Near both a and sigma;
                return *gf_iv_;
            }
            else
            {
                return *gf_iv_;
                // TODO near sigma; use GreensFunction2DRadInf
                // gf_tmp_ = std::make_unique<GreensFunction2DRadInf>();
                // return *gf_tmp_;
            }
        }
        else
        {
            // Sigma unreachable
            if (distance_from_shell < threshold_distance)
            {
                return *gf_iv_;
                // TODO near a; use GreensFunction2DAbs
                // gf_tmp_ = std::make_unique<GreensFunction2DAbs>();
            }
            else
            {
                return *gf_iv_;
                // TODO distant from both a and sigma;
                //gf_tmp_ = std::make_unique<GreensFunction2D>();
            }
        }

        // return *gf_tmp_;
    }

    bool create_updated_shell(const shell_matrix_type &smat, const World &world, ShellID sid1, ShellID sid2) override
    {
        // Create IV and CoM vectors from particles
        do_transform(world);

        THROW_UNLESS(no_space, r0() >= sigma());        // distance_from_sigma (pair gap) between %s and %s = %s < 0' % \(self.single1, self.single2, (self.r0 - self.sigma)))


        auto pos = com_ ;
        auto max_part_radius = gsl_max(particle1().radius(), particle2().radius());
        double min_radius = max_part_radius * GfrdCfg.MULTI_SHELL_FACTOR;
        double max_radius = std::min(smat.cell_size() / std::sqrt(8.0),                       // any angle cylinder must fit into cell-matrix! 2*sqrt(2)
                                     scaling::dist_to_plane_edge(pos, structure_.get()->id(), world));  // and not exceed its plane edges
        auto height = 2 * min_radius;
        auto D_r = D_tot();

        // Prevent domain creation if neither particle can move, since no events will be possible
        if (D_r == 0.0)
        {
            return false;
        }

        auto plane = structure_;
        THROW_UNLESS(not_found, plane != nullptr);
        auto unit_z = plane->shape().unit_z();

        std::vector<ShellID> ignored_shells = {sid1, sid2};
        std::vector<StructureID> ignored_structures = {structure_.get()->id(), world.get_def_structure_id()};

        auto start1 = pos - unit_z * height / 2;
        auto end1 = pos + unit_z * height / 2;
        auto radius = scaling::find_maximal_cylinder_radius(start1, end1, smat, world, ignored_structures, ignored_shells);

        // Ensure we don't exceed the matrix cell dimensions
        radius = std::min(radius, max_radius);
        radius /= GfrdCfg.SAFETY;

        // Prevent domain creation if there is not enough space
        if (radius - r0() - max_part_radius <= 0.0 || radius < min_radius)
        {
            return false;
        }

        // Calculate radii for center-of-motion vector R, and interparticle vector r
        determine_radii(radius);

        sid_pair_.second = Shell(domainID_, Cylinder(pos, radius, unit_z, height/2), Shell::Code::NORMAL);
        gf_com_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(D_R(), a_R()));
        gf_iv_ = std::make_unique<GreensFunction2DRadAbs>(GreensFunction2DRadAbs(D_r, k_total(), r0(), sigma(), a_r()));
        return true;
    }

    Vector3 create_com_vector(double r, RandomNumberGenerator &rng) const override
    {
        auto plane = structure_.get()->shape();
        auto vec = r * Vector2::random(rng);
        return vec.X() * plane.unit_x() + vec.Y() * plane.unit_y();
    }

    Vector3
    create_interparticle_vector(const PairGreensFunction &gf, double r, RandomNumberGenerator &rng) const override
    {
        auto theta = GreenFunctionHelper::draw_theta_wrapper(rng, gf, r, dt());

        // Note that r0 = length (old_iv)
        auto new_iv = (r/r0()) * Vector3::transformVector(iv_, Matrix4::createRotationA(theta, structure_.get()->shape().unit_z()));

        // Note that unit_z can point two ways rotating the vector clockwise or counterclockwise
        // Since theta is symmetric this doesn't matter.

        // Project the new_iv down on the unit vectors of the surface to prevent the particle from
        // leaving the surface due to numerical problem
        auto shape = structure_.get()->shape();
        auto unit_x = shape.unit_x(), unit_y = shape.unit_y();

        auto new_iv_x = unit_x * (Vector3::dot(new_iv, unit_x));
        auto new_iv_y = unit_y * (Vector3::dot(new_iv, unit_y));

        return new_iv_x + new_iv_y;
    }

    mutable std::unique_ptr<PairGreensFunction> gf_tmp_;
    std::shared_ptr<PlanarSurface> structure_;
    const int LD_MAX_ = 20;
};

// --------------------------------------------------------------------------------------------------------------------------------


class PairMixed2D3D : public Pair
{
public:
    PairMixed2D3D(const DomainID did, const particle_id_pair& pid_pair_2d, const particle_id_pair& pid_pair_3d,
            const std::shared_ptr<Structure>& structure_2d, std::shared_ptr<Structure>& structure_3d, const shell_id_pair& sid_pair, const ReactionRules& reactions, const World &world)
            : Pair(did, pid_pair_2d, pid_pair_3d, sid_pair, reactions),
            pid_pair_2d_(pid_pair_2d), pid_pair_3d_(pid_pair_3d), structure_3d_(structure_3d), world_(world)
    {
//        THROW_UNLESS(illegal_state, sid_pair.second.shape() == Shell::Shape::CYLINDER);
        THROW_UNLESS(illegal_state, structure_2d.get()->id() == pid_pair_2d.second.structure_id());
        THROW_UNLESS(illegal_state, structure_3d.get()->id() == pid_pair_3d.second.structure_id());

        structure_2d_ = std::dynamic_pointer_cast<PlanarSurface>(structure_2d);
        THROW_UNLESS_MSG(illegal_state, structure_2d_ != nullptr, "PairPlanar created with structure that is not a PlanarSurface");


        // The scaling factor needed to make the anisotropic diffusion problem
        // of the IV into an isotropic one. 3D particle is the only contributor
        // to diffusion normal to the plane.
        scaling_factor_ = sqrt(D_tot() / pid_pair_3d.second.D()); // D_r == D_tot
        sqrt_DRDr_ = std::sqrt((2*D_R())/(3*D_tot())); // D_r == D_tot
    }

    const char* type_name() const override { return "PairMixed2D3D"; }

    bool create_updated_shell(const shell_matrix_type &smat, const World &world, ShellID sid1, ShellID sid2) override
    {
        // Create IV and CoM vectors from particles
        do_transform(world);

        THROW_UNLESS(no_space, r0() >= sigma());        // distance_from_sigma (pair gap) between %s and %s = %s < 0' % \(self.single1, self.single2, (self.r0 - self.sigma)))

        auto pos_2d = particle1().position();
        auto pos_3d = world.cyclic_transpose(particle2().position(), pos_2d);
        auto max_part_radius = gsl_max(particle1().radius(), particle2().radius());
        double min_radius = max_part_radius * GfrdCfg.MULTI_SHELL_FACTOR;

        //TODO: cylinder should be centered around CoM, not 2D particle
        double max_radius = std::min(smat.cell_size() / std::sqrt(8.0),                       // any angle cylinder must fit into cell-matrix! 2*sqrt(2)
                                     scaling::dist_to_plane_edge(pos_2d, structure_2d_.get()->id(), world));  // and not exceed its plane edges
        auto radius = max_radius;
        auto height = 2 * min_radius;
        auto D_r = D_tot();
        auto D_2D = particle1().D(), D_3D = particle2().D();
        auto radius2D = particle1().radius(), radius3D = particle2().radius();

        particle_surface_dist_ = structure_2d_->distance(pos_3d);

        std::vector<ShellID> ignored_shells = {sid1, sid2};
        std::vector<StructureID> ignored_structures = {particle1().structure_id(), particle1().structure_id(), world.get_def_structure_id()};

        auto plane = structure_2d_;
        THROW_UNLESS(not_found, plane != nullptr);
        auto unit_z = (plane->shape().unit_z() * Vector3::dot(iv(), plane->shape().unit_z())).normal();

        auto height_through_surface = radius2D * GfrdCfg.SINGLE_SHELL_FACTOR;
        auto height_to_surface = get_distance_to_surface();
        auto static_height = height_through_surface + height_to_surface;
        auto base_pos = pos_2d - height_through_surface * unit_z;
        auto max_dynamic_height = scaling::find_maximal_cylinder_height<shell_matrix_type>(base_pos, unit_z, static_height, scaling_factor_, smat, world, ignored_structures, ignored_shells);

        radius = std::min(max_dynamic_height * scaling_factor_, max_radius);
        radius /= GfrdCfg.SAFETY;

        if (radius < min_radius) return false; // no space for Mixed2D3D domain.. it will be discarded, and a single or multi is created instead!

        if (D_R() == 0.0)
        {
            auto a_r_3D = (radius - particle2().radius());
            a_r_ = a_r_3D;
        }
        else
        {
            auto a_r_2D = (radius - radius2D + r0()*sqrt_DRDr_) / (sqrt_DRDr_ + (D_2D/D_tot()));
            auto a_r_3D = (radius - radius3D + r0()*sqrt_DRDr_) / (sqrt_DRDr_ + (D_3D/D_tot()));
            a_r_ = std::min(a_r_2D, a_r_3D);
        }

        // Height is set to a ratio of the radius, such that diffusion of the IV becomes isotropic
        height = (a_r_ / scaling_factor_) + static_height;

        //TODO: cylinder should be centered around CoM, not 2D particle
        auto cylinder_pos = pos_2d + unit_z * (height/2 - height_through_surface);

        // Calculate radii for center-of-motion vector R, and interparticle vector r
        determine_radii(radius, height/2);

        // We need to check whether the final a_r is large enough to fit r0.
        // If it is larger, but only just so, the next time delta will be in the
        // order of nanoseconds, so it is best to skip creating this domain.
        if (a_r() < r0() * GfrdCfg.SAFETY)
        {
            return false;
        }

        sid_pair_.second = Shell(domainID_, Cylinder(cylinder_pos, radius, unit_z, height/2), Shell::Code::NORMAL);
        gf_com_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(D_R(), a_R()));
        gf_iv_ = std::make_unique<GreensFunction3DRadAbs>(GreensFunction3DRadAbs(D_r, k_total(), r0(), sigma(), a_r()));
        return true;
    }

    Vector3 create_com_vector(double r, RandomNumberGenerator &rng) const override
    {
        auto plane = structure_2d_.get()->shape();
        auto vec = r * Vector2::random(rng);
        return vec.X() * plane.unit_x() + vec.Y() * plane.unit_y();
    }

    Vector3 create_interparticle_vector(const PairGreensFunction &gf, double r, RandomNumberGenerator &rng) const override
    {
        auto old_iv = iv_;
        auto new_iv = iv_;
        auto theta = GreenFunctionHelper::draw_theta_wrapper(rng, gf, r, dt());
        if (theta == 0.0)
        {
            new_iv = old_iv;
        }
        else if (!feq(std::fmod(theta, M_PI), 0.0))
        {
            // Rotate the old iv to the new theta
            auto rotation_axis = Vector3(-old_iv.Y(), old_iv.X(), 0.0); // Crossproduct against z-axis
            new_iv = Vector3::transformVector(old_iv, Matrix4::createRotationA(theta, rotation_axis.normal()));

            // Rotate the new iv around the old iv with angle phi
            auto phi = rng.uniform(0.0, 1.0) * 2 * M_PI;
            rotation_axis = old_iv.normal();
            new_iv = Vector3::transformVector(new_iv, Matrix4::createRotationA(phi, rotation_axis.normal()));
        }
        else
        {
            // theta == pi -> just mirror the old iv
            new_iv = -old_iv;
        }

        // Adjust length of the vector to new r
        // Note that r0 = old_iv.length()
        new_iv = (r/r0()) * new_iv;
        return new_iv;
    }

    std::pair<position_structid_pair, position_structid_pair>
    do_back_transform(Vector3 com, Vector3 iv, double D1, double D2, double radius1, double radius2, StructureID s1,
                      StructureID s2, Vector3 unit, const World &world) const override
    {

        auto weight1 = D1 / D_tot(), weight2 = D2 / D_tot();
        auto z_backtransform = sqrt(D2 / D_tot());

        // Like in create_interparticle_vector(), we assume the 2D surface is a PlanarSurface
        auto surface = dynamic_cast<PlanarSurface*>(structure_2d_.get());
        THROW_UNLESS_MSG(illegal_state, surface != nullptr, "PairMixed2D3D has a non-planar 2D surface, which is not supported.");
        auto shape = surface->shape();

        // Get the coordinates of the iv relative to the system of the surface (or actually the shell)
        auto iv_x = shape.unit_x() * Vector3::dot(iv, shape.unit_x());
        auto iv_y = shape.unit_y() * Vector3::dot(iv, shape.unit_y());

        auto iv_x_length = iv_x.length(), iv_y_length = iv_y.length();

        // Reflect the coordinates in the unit_z direction back to the side of the membrane
        // where the domain is. Note that it's implied that the origin of the coordinate system lies in the
        // plane of the membrane.
        auto iv_z_length =  fabs(Vector3::dot(iv, unit)); // FIXME maybe first project the shell unit_z onto the surface unit vector to prevent numerical problems?

        // Do the reverse scaling in z-direction
        auto iv_z_length_backtransform = iv_z_length * z_backtransform;

        // The following is to avoid overlaps in case that the z-component of the IV is scaled down
        // so much that it would cause a particle overlap.
        // This can happen because the inner boundary in the space of transformed coordinates is not a
        // prolate spheroid, as it should, but approximated by a sphere. This can lead to IV exit points
        // within the prolate spheroid, i.e. within the sphere with radius sigma = radius1 + radius2 in
        // the real space after inverse transform.  FIXME This works, but is still somewhat of a HACK!
        auto min_iv_length = (radius1 + radius2) * GfrdCfg.MINIMAL_SEPARATION_FACTOR;

        if (iv_x_length * iv_x_length + iv_y_length * iv_y_length + iv_z_length_backtransform * iv_z_length_backtransform < min_iv_length * min_iv_length)
        {
            auto z_safety_factor_sq = (min_iv_length * min_iv_length - iv_x_length*iv_x_length - iv_y_length*iv_y_length)
                    / (iv_z_length_backtransform * iv_z_length_backtransform);

            THROW_UNLESS_MSG(illegal_state, z_safety_factor_sq >= 0.0, "Can't have negative z safety factor.");
            Log("GFRD").warn() << "PairMixed2D3D: applying z safety factor to enlarge too small interparticle vector.";

            iv_z_length_backtransform = sqrt(z_safety_factor_sq) * iv_z_length_backtransform;
        }

        // TODO: also check for membrane overlapping?

        // Construct the z-vector-component and the new position vectors
        auto iv_z = unit * iv_z_length_backtransform;
        auto new_2d_pos = com - weight1 * (iv_x + iv_y);
        auto new_3d_pos = com + (weight2 * (iv_x + iv_y)) + iv_z;

        // Make pairs
        auto pair_2d = std::make_pair(new_2d_pos, structure_2d_.get()->id());
        auto pair_3d = std::make_pair(new_3d_pos, structure_3d_.get()->id());
        return std::make_pair(pair_2d, pair_3d);
    }

    const PairGreensFunction &choose_pair_greens_function() const override
    {
        return *gf_iv_;
    }

    void determine_radii(const double shell_radius, const double shell_half_length)
    {
        // Determines the dimensions of the domains used for the Green's functions
        // from the dimensions of the cylindrical shell.
        // Note that the function assumes that the cylinder is dimensioned properly.

        auto particle_2d = pid_pair_2d_.second, particle_3d = pid_pair_3d_.second;
        auto radius_2d = particle_2d.radius(), radius_3d = particle_3d.radius();
        auto D_2d = particle_2d.D(), D_3d = particle_3d.D();

        // The distance the cylinder protrudes through the PlanarSurface. We also call this
        // side of the cylinder the bottom.
        auto z_left = radius_2d * GfrdCfg.SINGLE_SHELL_FACTOR;

        // The distance between the 3D-diffusing particle and the cylinder top cap.
        auto z_right = 2.0 * shell_half_length - get_distance_to_surface() - z_left;

        // Partition the space in the protective domain over the IV and the CoM domains
        // The outer bound of the interparticle vector is set by the particle diffusing in 3D via:
        a_r_ = (z_right - 2*radius_3d) * scaling_factor_;

        // Calculate the maximal displacement of a particle given an interparticle vector bound a_r
        // taking into account the radius for both the 2D and 3D particle. This sets the remaining
        // space for the CoM diffusion.
        auto space_for_iv = gsl_max( (D_2d/D_tot()) * a_r() + radius_2d,
                                     (D_3d/D_tot()) * a_r() + radius_3d);

        // It then determines the CoM vector bound via
        a_R_ = shell_radius - space_for_iv;
        if(feq(a_R(), 0.0, shell_radius)) {
            a_R_ = 0.0;
            Log("GFRD").info() << "(PairMixed2D3D) determine_radii: setting a_R to zero";
        }

        // TODO: Add sanity checks from pair.py
    }

    void do_transform(const World& world) override
    {
        Vector3 pos1 = pid_pair_2d_.second.position();
        Vector3 pos2 = pid_pair_3d_.second.position();
        Vector3 pos2c = world.cyclic_transpose(pos2, pos1);

        // The CoM is calculated in a similar way to a normal 3D pair...
        Vector3 com = D_tot() > 0 ? (pid_pair2_.second.D() * pos1 + pid_pair1_.second.D() * pos2c) / D_tot() : 0.5 * (pos1 + pos2c);

        // ...and then projected onto the plane to make sure the CoM is in the surface
        com = world.cyclic_transpose(com, structure_2d_.get()->position());
        auto projection = structure_2d_.get()->project_point(com);
        com_ = world.apply_boundary(projection.first);

        // IV calculation
        auto iv = pos2c - pos1;

        // Rescale the IV in the axis normal to the plane
        auto plane = structure_2d_;
        auto iv_z_component = Vector3::dot(iv, plane->shape().unit_z());
        auto iv_z = plane->shape().unit_z() * iv_z_component;
        auto new_iv_z = plane->shape().unit_z() * (iv_z_component * scaling_factor_);

        iv_ = iv - iv_z + new_iv_z;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    double get_distance_to_surface() const
    {
        // This is the distance from the particle to the PlanarSurface.
        // Note that this is the distance from the hull of the particle, not its center.
        return particle_surface_dist_;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    double get_distance_to_escape(double half_length) const
    {
        // This calculates the distance from the particle to the flat boundary away from the surface.
        // Note that this is the distance from the hull of the particle, not its center.
        auto cylinder_length = half_length*2;
        auto particle_radius = particle2().radius();
        // Portion of cylinder behind planar surface
        auto cylinder_left = (particle_radius * GfrdCfg.SINGLE_SHELL_FACTOR);
        // Portion of cylinder in front of planar surface
        auto cylinder_right = cylinder_length - cylinder_left;

        auto particle_to_right_boundary = cylinder_right - particle_surface_dist_ - particle_radius;
        return particle_to_right_boundary;
    }

    double sigma() const override {
        // Rescale sigma to correct for the rescaling of the coordinate system.
        // This is the sigma that is used for the evaluation of the Green's function and is in this case
        // slightly different from the sums of the radii of the particles. We rescale sigma in a way
        // that the total flux over the reactive surface of the prolate spheroidal boundary with greater radius
        // sigma is equal to the total flux over a spherical reactive surface with radius = rescaled sigma.
        auto xi = scaling_factor_;
        auto sigma = particle1().radius() + particle2().radius();

        double rescaled_sigma;

        if (feq(xi, 1.0, 1.0)) {
            // If the scaling factor happens to be equal to one (i.e. no rescaling) alpha and sin(alpha)
            // in the rescaling formula for sigma below will be zero and cause zero-division problems.
            // However, if xi=1 we know that we do not have to rescale sigma at all.
            // The case xi=1 usually means that the 2D particle is static (D_2D=0).
            rescaled_sigma = sigma;
        } else {
            auto xi_inverse = 1.0/xi;
            auto alpha = acos(xi_inverse);
            rescaled_sigma = fabs(sigma * sqrt(0.5 + (alpha * xi/(2.0*sin(alpha)))));
        }

        return rescaled_sigma;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    double center_particle_offset(double half_length) const
    {
        auto particle_radius = particle2().radius();
        // Distance between the particle center and cylinder center
        return (half_length - get_distance_to_escape(half_length) - particle_radius);
    }

private:
    friend class Persistence;

    const World &world_;

    const particle_id_pair pid_pair_2d_;
    const particle_id_pair pid_pair_3d_;
    std::shared_ptr<PlanarSurface> structure_2d_;
    std::shared_ptr<Structure> structure_3d_;

    double particle_surface_dist_;
    double scaling_factor_;
    double sqrt_DRDr_;
};