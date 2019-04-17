#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <set>
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

    virtual EventType draw_iv_event_type(RandomNumberGenerator& rng) const = 0;

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

    EventType draw_iv_event_type(RandomNumberGenerator& rng) const override 
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

    EventType draw_iv_event_type(RandomNumberGenerator& rng) const override 
    { 
       UNUSED(rng);
       THROW_EXCEPTION(illegal_state, "SingleCylindrical does not support iv events.");
    }

    Vector3 create_position_vector(double r, RandomNumberGenerator& rng) const override
    {
        // draw random circle and rotate it in the plane
        Vector2 c = r * Vector2::random(rng);
        auto m = Matrix4::createRotationAB(structure_->shape().unit_z(), Vector3::uy);
        return m.multiply(Vector3(c.X(), 0.0, c.Y()));
    }

    GFRD_EXPORT virtual bool create_updated_shell(const shell_matrix_type& smat, const World& world) override;

private:
    std::shared_ptr<PlanarSurface> structure_;

};

// --------------------------------------------------------------------------------------------------------------------------------
