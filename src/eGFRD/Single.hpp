#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <set>
#include "scaling.hpp"
#include "Particle.hpp"
#include "ParticleID.hpp"
#include "DomainID.hpp"
#include "Domain.hpp"
#include "ShellID.hpp"
#include "Shell.hpp"
#include "Cylinder.hpp"
#include "ReactionRule.hpp"
#include "ReactionRuleCollection.hpp"
#include "randomNumberGenerator.hpp"
#include <GreensFunction.hpp>
#include "MatrixSpace.hpp"
#include "World.hpp"
#include "ShellCreateUtils.hpp"
#include "GreenFunctionHelpers.hpp"
#include "exceptions.hpp"
#include "Matrix4.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Single : public Domain
{
public:
    using particle_id_pair = Particle::particle_id_pair;
    using shell_id_pair = Shell::shell_id_pair;

    using ReactionRules = ReactionRuleCollection::reaction_rule_set;
    using shell_matrix_type = MatrixSpace<Shell, ShellID, World::particle_matrix_type::SizeX, World::particle_matrix_type::SizeY, World::particle_matrix_type::SizeZ>;
    using shell_distance_checker = ShellCreateUtils::shell_distance_check<shell_matrix_type>;
    using surface_distance_check = ShellCreateUtils::surface_distance_check;
    using position_structid_pair = World::position_structid_pair;

    // --------------------------------------------------------------------------------------------------------------------------------

    explicit Single(const DomainID did, const particle_id_pair& pid_pair, const shell_id_pair& sid_pair, const ReactionRules& rrules) noexcept
        : Domain(did), pid_pair_(pid_pair), sid_pair_(sid_pair), rrules_(rrules), rrule_(nullptr), gf_(nullptr) {}

    // --------------------------------------------------------------------------------------------------------------------------------

    const particle_id_pair& pip() const { return pid_pair_; }            // this returns a value copy of the Particle, and it should! (because Domain is about to be erased)
    ParticleID particle_id() const { return pid_pair_.first; }
    const Particle& particle() const { return pid_pair_.second; }

    const shell_id_pair& shell_pair() const { return sid_pair_; }        // this returns a value copy of the Shell !
    ShellID shell_id() const { return sid_pair_.first; }
    const Shell& shell() const { return sid_pair_.second; }

    // --------------------------------------------------------------------------------------------------------------------------------

    size_t num_shells() const override { return 1; }

    Multiplicity multiplicity() const override { return Multiplicity::SINGLE; }

    virtual shell_id_pair_ref get_shell() const override { return std::make_pair(sid_pair_.first, std::ref<const Shell>(sid_pair_.second)); }

    virtual std::vector<shell_id_pair_ref> get_shell_list() const override { return std::vector<shell_id_pair_ref >({ get_shell() }); }

    // --------------------------------------------------------------------------------------------------------------------------------

    std::string as_string() const override { return make_string() << ", " << pid_pair_.first << ", " << pid_pair_.second << ", " << sid_pair_.first << ", " << sid_pair_.second; }

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual void determine_next_event(RandomNumberGenerator& rng)
    {
        auto t1 = draw_escape_time(rng);
        auto t2 = draw_reaction_time(rng);

        if (t1 < t2)
        {
            dt_ = t1;
            eventType_ = EventType::SINGLE_ESCAPE;
        }
        else
        {
            dt_ = t2;
            eventType_ = EventType::SINGLE_REACTION;
        }
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual double draw_escape_time(RandomNumberGenerator& rng)
    {
        return GreenFunctionHelper::draw_time_wrapper(rng, *gf_.get());
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual bool create_updated_shell(const shell_matrix_type& smat, const World& world) = 0;

    // --------------------------------------------------------------------------------------------------------------------------------

    double draw_reaction_time(RandomNumberGenerator& rng) const
    {
        return GreenFunctionHelper::draw_reaction_time(k_total(), rng);
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    double k_total() const
    {
        // calculates the total rate for a list of reaction rules
        // The probability for the reaction to happen is proportional to the sum of the rates of all the possible reaction types.
        double k_tot = 0;
        for (auto& rr : rrules_)
            k_tot += rr.getK();
        return k_tot;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    const ReactionRule& draw_reaction_rule(RandomNumberGenerator& rng) const
    {
        // draws a reaction rules out of a list of reaction rules based on their relative rates
        if (rrule_ != nullptr) return *rrule_;

        double rnd = rng.uniform(0, k_total());
        ReactionRules::iterator i = rrules_.begin();
        double k_sum = 0.0;
        while (k_sum + i->getK() < rnd) { k_sum += i->getK();  ++i; }

        rrule_ = &*i;
        return *rrule_;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual double get_inner_a() const = 0;

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual Vector3 create_position_vector(double r, RandomNumberGenerator& rng) const = 0;

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual EventType draw_iv_event_type(RandomNumberGenerator& rng) = 0;

    // --------------------------------------------------------------------------------------------------------------------------------

    position_structid_pair draw_new_position(RandomNumberGenerator& rng) const
    {
        Vector3 newpos, oldpos = particle().position();

        if (particle().D() == 0)
            newpos = oldpos;
        else if (eventType() == EventType::SINGLE_REACTION && draw_reaction_rule(rng).get_products().size() == 0)
            newpos = oldpos;
        else
        {
            double r;
            if (eventType() == EventType::SINGLE_ESCAPE)
                r = get_inner_a();
            else
                // Note that in case of 1D diffusion r has a direction. It is the lateral displacement and can be positive or negative.
                // In other cases r is always positive and denotes radial displacement from the center.
                r = GreenFunctionHelper::draw_r_wrapper(rng, *gf_.get(), dt_, get_inner_a());

            // calculate other coordinates
            Vector3 displacement = create_position_vector(r, rng);

            // Add displacement to shell.shape.position, not to particle.position.
            // This distinction is important only in the case of an asymmetric 1D domain(r0 != 0, or drift), since draw_r always
            // returns a distance relative to the centre of the shell(r = 0), not relative to r0.

            // note that we need to make sure that the shell.shape.position and displacement vector
            // are in the structure to prevent the particle leaving the structure due to numerical errors
            newpos = shell().position() + displacement;
        }

        return std::make_pair(newpos, particle().structure_id());
    }

    // --------------------------------------------------------------------------------------------------------------------------------

protected:
   friend class Persistence;

    particle_id_pair pid_pair_;
    shell_id_pair sid_pair_;
    ReactionRules rrules_;
    mutable const ReactionRule* rrule_;    // once a rule is drawn, remember it.
    std::unique_ptr<GreensFunction> gf_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class SingleSpherical : public Single
{
public:
    SingleSpherical(const DomainID did, const particle_id_pair& pid_pair, const shell_id_pair& sid_pair, const ReactionRules& rrules)
        : Single(did, pid_pair, sid_pair, rrules)
    {
        THROW_UNLESS(illegal_state, sid_pair.second.shape() == Shell::Shape::SPHERE);
    }

    const char* type_name() const override { return "SingleSpherical"; }

    double get_inner_a() const override { return sid_pair_.second.get_sphere().radius() - pid_pair_.second.radius(); }

    EventType draw_iv_event_type(RandomNumberGenerator& rng) override
    {
       UNUSED(rng);
       THROW_EXCEPTION(illegal_state, "SingleSpherical does not support iv events."); 
    }

    Vector3 create_position_vector(double r, RandomNumberGenerator& rng) const override
    {
        return r * Vector3::random(rng);
    }

    GFRD_EXPORT bool create_updated_shell(const shell_matrix_type& smat, const World& world) override;

};

// --------------------------------------------------------------------------------------------------------------------------------

class SingleCylindrical : public Single
{
public:
    SingleCylindrical(const DomainID did, const particle_id_pair& pid_pair, const shell_id_pair& sid_pair, const ReactionRules& rrules)
        : Single(did, pid_pair, sid_pair, rrules)
    {
        if (sid_pair.second.code() == Shell::Code::INIT) THROW_UNLESS(illegal_state, sid_pair.second.shape() == Shell::Shape::SPHERE);
        if (sid_pair.second.code() == Shell::Code::NORMAL) THROW_UNLESS(illegal_state, sid_pair.second.shape() == Shell::Shape::CYLINDER);
    }

    const char* type_name() const override { return "SingleCylindrical"; }

    double get_inner_a() const override { return sid_pair_.second.get_cylinder().radius() - pid_pair_.second.radius(); }

    EventType draw_iv_event_type(RandomNumberGenerator& rng) override
    { 
       UNUSED(rng);
       THROW_EXCEPTION(illegal_state, "SingleCylindrical does not support iv events.");
    }

    Vector3 create_position_vector(double r, RandomNumberGenerator& rng) const override
    {
        // draw random circle and rotate it in the plane
        Vector2 c = r * Vector2::random(rng);
        auto m = Matrix4::createRotationAB(structure_.get()->shape().unit_z(), Vector3::uy);
        return m.multiply(Vector3(c.X(), 0.0, c.Y()));
    }

    GFRD_EXPORT bool create_updated_shell(const shell_matrix_type& smat, const World& world) override;

protected:
    std::shared_ptr<PlanarSurface> structure_;

};

// --------------------------------------------------------------------------------------------------------------------------------

class SingleInteraction : public Single {
public:
    using InteractionRules = ReactionRuleCollection::interaction_rule_set;
    using EventKind = GreensFunction::EventKind;

    SingleInteraction(const DomainID did, const particle_id_pair &pid_pair,
                      const shell_id_pair &sid_pair, const ReactionRules &reaction_rules,
                      const InteractionRules &interaction_rules)
            : Single(did, pid_pair, sid_pair, reaction_rules),
              interaction_rules_(interaction_rules)
    {
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual double draw_iv_event_time(RandomNumberGenerator &rng)
    {
        return GreenFunctionHelper::draw_time_wrapper(rng, *gf_iv_);
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    void determine_next_event(RandomNumberGenerator& rng) override
    {
        auto t_escape = draw_escape_time(rng);
        auto t_reaction = draw_reaction_time(rng);
        auto t_iv = draw_iv_event_time(rng);

        if (t_escape < t_reaction && t_escape < t_iv)
        {
            dt_ = t_escape;
            eventType_ = EventType::SINGLE_ESCAPE;
        }
        else if(t_reaction < t_escape && t_reaction < t_iv)
        {
            dt_ = t_reaction;
            eventType_ = EventType::SINGLE_REACTION;
        }
        else
        {
            dt_ = t_iv;
            eventType_ = EventType::IV_EVENT; // Particular IV event type is determined just-in-time in eGFRDSimulator::process_single_event()
        }
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    EventType draw_iv_event_type(RandomNumberGenerator& rng) override
    {
        auto type = GreenFunctionHelper::draw_eventtype_wrapper(rng, *gf_iv_, dt_);
        iv_event_kind_ = type;

        switch(type)
        {
            case GreensFunction::EventKind::IV_ESCAPE:
                return EventType::SINGLE_ESCAPE;
            case GreensFunction::EventKind::IV_REACTION:
                return EventType::SINGLE_INTERACTION;
            default:
                THROW_EXCEPTION(not_implemented, " Unknown IV event type");
        }
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    double interaction_k_total() const
    {
        // calculates the total rate for a list of interaction rules
        // The probability for the interaction to happen is proportional to the sum of the rates of all the possible interactions.
        double k_tot = 0;
        for (auto& rr : interaction_rules_)
            k_tot += rr.getK();
        return k_tot;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    const InteractionRule& draw_interaction_rule(RandomNumberGenerator& rng) const
    {
        // draws an interaction rule out of a list of reaction rules based on their relative rates
        double rnd = rng.uniform(0, interaction_k_total());
        auto i = interaction_rules_.begin();
        double k_sum = 0.0;
        while (k_sum + i->getK() < rnd) { k_sum += i->getK();  ++i; }

        return *i;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    virtual std::shared_ptr<const Structure*> get_interacting_structure() const = 0;

    // --------------------------------------------------------------------------------------------------------------------------------

protected:
    const InteractionRules& interaction_rules_;
    std::unique_ptr<GreensFunction> gf_iv_;
    EventKind iv_event_kind_;
};

// --------------------------------------------------------------------------------------------------------------------------------


class SinglePlanarInteraction : public SingleInteraction
{
public:

    SinglePlanarInteraction(const DomainID did, const particle_id_pair& pid_pair, const PlanarSurface& structure,
            const shell_id_pair& sid_pair, const ReactionRules& reaction_rules, const InteractionRules& interaction_rules)
            : SingleInteraction(did, pid_pair, sid_pair, reaction_rules, interaction_rules), interacting_structure_(structure),
              structure_pointer_(std::make_shared<const Structure*>(&interacting_structure_))
    {
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    const char* type_name() const override { return "SinglePlanarInteraction"; }

    // --------------------------------------------------------------------------------------------------------------------------------

    Vector3 create_position_vector(double r, RandomNumberGenerator& rng) const override
    {
        // r is drawn from the 2D diffusion Green's function
        Vector2 c = r * Vector2::random(rng);

        // Express the r vector in the unit vectors of the surface to make sure the particle is
        // parallel to the surface (due to numerical problem)
        auto shape = interacting_structure_.shape();
        auto vector_r = c.X() * shape.unit_x() + c.Y() * shape.unit_y();

        // Calculate z vector
        auto shell = sid_pair_.second.get_cylinder();
        auto half_length = shell.half_length();

        // Calculate radiating and absorbing boundary condition distances
        auto z_surface = -get_distance_to_surface(); // (Negative) distance from particle to planar surface
        auto z_absorb = get_distance_to_escape(half_length); // Distance from particle to disk away from planar surface
        double z_displacement = 0;

        if (eventType_ != EventType::IV_EVENT)
        {
            // Normal escape event, z_displacement is drawn from 1D IV Green's function
            z_displacement = GreenFunctionHelper::draw_r_wrapper(rng, *gf_iv_, dt_, z_absorb, z_surface); // last two args= a, sigma
        }
        else if (iv_event_kind_ == EventKind::IV_ESCAPE)
        {
            // IV escape event
            z_displacement = z_absorb;
        }
        else if (iv_event_kind_ == EventKind::IV_REACTION)
        {
            // Interaction event
            z_displacement = z_surface;
        }

        // z_displacement is now relative to the particle center, but draw_new_position() requires it relative to
        // the center of the shell.
        auto offset = center_particle_offset(shell.half_length());
        z_displacement -= offset;

        // Express the vector_z in the unit vectors of the surface to prevent the particle from
        // leaving the surface due to numerical problem.
        auto vector_z = z_displacement * shell.unit_z();

        // The new position is relative to the center of the shell
        // note that we need to make sure that vector_r and vector_z
        // are correct (in the structure) to prevent the particle leaving the structure due to numerical errors
        return vector_r + vector_z;
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
        auto particle_radius = pid_pair_.second.radius();
        // Portion of cylinder behind planar surface
        auto cylinder_left = (particle_radius * GfrdCfg.SINGLE_SHELL_FACTOR);
        // Portion of cylinder in front of planar surface
        auto cylinder_right = cylinder_length - cylinder_left;

        auto particle_to_right_boundary = cylinder_right - particle_surface_dist_ - particle_radius;
        return particle_to_right_boundary;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    double center_particle_offset(double half_length) const
    {
        auto particle_radius = pid_pair_.second.radius();
        // Distance between the particle center and cylinder center
        return (half_length - get_distance_to_escape(half_length));
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    double get_inner_a() const override { return sid_pair_.second.get_cylinder().radius() - pid_pair_.second.radius(); }

    // --------------------------------------------------------------------------------------------------------------------------------

    std::shared_ptr<const Structure*> get_interacting_structure() const override
    {
        return structure_pointer_;
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    GFRD_EXPORT bool create_updated_shell(const shell_matrix_type& smat, const World& world) override;

    // --------------------------------------------------------------------------------------------------------------------------------

protected:
    const PlanarSurface& interacting_structure_;
    const std::shared_ptr<const Structure*> structure_pointer_;
    double particle_surface_dist_ = 0.0;
};