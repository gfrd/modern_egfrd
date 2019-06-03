#include "Single.hpp"
#include "exceptions.hpp"
#include "ShellCreateUtils.hpp"
#include "EGFRDSimulator.hpp"
#include <GreensFunction2DAbsSym.hpp>
#include <GreensFunction3DAbsSym.hpp>
#include <GreensFunction1DRadAbs.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT bool SingleSpherical::create_updated_shell(const shell_matrix_type& smat, const World& world)
{
   auto pos = pid_pair_.second.position();
   double min_radius = particle().radius() * GfrdCfg.MULTI_SHELL_FACTOR;
   double max_radius = smat.cell_size() / 2.0;      // assure spheres cannot overlap
   auto radius = max_radius;

   // when not defusing, use minimal shell (we're not going anywhere, lease space for all others)
   if (pid_pair_.second.D() == 0.0 && pid_pair_.second.v() == 0.0) radius = min_radius;
   
   // check distances to surfaces, ignore def.struct and particle.structure
   for (auto s : world.get_structures())
   {
      if (s->id() != pid_pair_.second.structure_id() && s->id() != world.get_def_structure_id())    // ignore structure that particle is attached to
      {
         auto transposed = world.cyclic_transpose(pos, s->position());
         double distance = s->distance(transposed);
         radius = std::min(radius, distance);
      }
   }

   // Limit radius to space to closest neighbour
   shell_distance_checker sdc(shell_id(), pos, radius, shell_distance_checker::Construct::SHARE5050);
   CompileConfigSimulator::TBoundCondition::each_neighbor(smat, sdc, pos);
   radius = std::min(radius, sdc.distance());

   // Limit radius to space to closest surface
   auto sudc = surface_distance_check(world, std::vector<StructureID>{world.get_def_structure_id()}, pid_pair_.second);
   sudc.measure_distances(world.get_structures());
   radius = std::min(radius, sudc.distance());

   radius /= GfrdCfg.SAFETY;

   if (radius < min_radius) return false;             // no space for Single domain.. it will be discarded, and a Multi is created instead!

   sid_pair_.second = Shell(domainID_, Sphere(pos, radius), Shell::Code::NORMAL);
   gf_ = std::make_unique<GreensFunction3DAbsSym>(GreensFunction3DAbsSym(pid_pair_.second.D(), get_inner_a()));
   return true;
}

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT bool SingleCylindrical::create_updated_shell(const shell_matrix_type& smat, const World& world)
{
   auto structure_id = pid_pair_.second.structure_id();
   auto pos = pid_pair_.second.position();
   double min_radius = particle().radius() * GfrdCfg.MULTI_SHELL_FACTOR;
   double max_radius = smat.cell_size() / std::sqrt(8.0);         // any angle cylinder must fit into cell-matrix! 2*sqrt(2)
   auto radius = max_radius;
   auto height = 2 * min_radius;

   auto search_distance = height + radius;

   // TODO: use actual geometric overlap check
   shell_distance_checker sdc(shell_id(), pos, search_distance, shell_distance_checker::Construct::SHARE5050);
   CompileConfigSimulator::TBoundCondition::each_neighbor(smat, sdc, pos);
   radius = std::min(radius, sdc.distance());

    // Limit radius to space to closest surface
    auto sudc = surface_distance_check(world, std::vector<StructureID>{world.get_def_structure_id(), structure_id}, pid_pair_.second);
    sudc.measure_distances(world.get_structures());
    radius = std::min(radius, sudc.distance());

   radius /= GfrdCfg.SAFETY;

   if (radius < min_radius) return false;             // no space for Single domain.. it will be discarded, and a Multi is created instead!

   auto structure = world.get_structure(pid_pair_.second.structure_id());
   auto *plane = dynamic_cast<PlanarSurface*>(structure.get());
   THROW_UNLESS(not_found, plane != nullptr);
   auto unit_z = plane->shape().unit_z();
   structure_ = std::dynamic_pointer_cast<PlanarSurface>(structure);

   sid_pair_.second = Shell(domainID_, Cylinder(pos, radius, unit_z, height/2), Shell::Code::NORMAL);
   gf_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(pid_pair_.second.D(), get_inner_a()));
   return true;
}

// --------------------------------------------------------------------------------------------------------------------------------


GFRD_EXPORT bool SinglePlanarInteraction::create_updated_shell(const shell_matrix_type& smat, const World& world)
{
    auto structure_id = pid_pair_.second.structure_id();
    auto pos = pid_pair_.second.position();
    auto transposed = world.cyclic_transpose(pos, interacting_structure_.position());
    auto particle_radius = pid_pair_.second.radius();
    particle_surface_dist_ = interacting_structure_.distance(transposed) - particle_radius;

//    THROW_UNLESS_MSG(illegal_state, particle_surface_dist_ >= 0.0, "Particle distance to interacting surface should be positive");

    double min_radius = particle().radius() * GfrdCfg.SINGLE_SHELL_FACTOR;
    double max_radius = smat.cell_size() / std::sqrt(8.0);         // any angle cylinder must fit into cell-matrix! 2*sqrt(2)

    //double max_height = smat.cell_size() / std::sqrt(2.0);
    auto radius = max_radius;

    // TODO: check distances with actual geometric overlap checking
    shell_distance_checker sdc(shell_id(),  pos, radius, shell_distance_checker::Construct::SHARE5050);
    CompileConfigSimulator::TBoundCondition::each_neighbor(smat, sdc, pos);
    radius = std::min(radius, sdc.distance());

    // Limit radius to space to closest surface
    auto sudc = surface_distance_check(world, std::vector<StructureID>{world.get_def_structure_id(), interacting_structure_.id()}, pid_pair_.second);
    sudc.measure_distances(world.get_structures());
    radius = std::min(radius, sudc.distance());

    radius /= GfrdCfg.SAFETY;

    if (radius < min_radius) return false;             // no space for Single domain.. it will be discarded, and a Multi is created instead!

    auto plane = &interacting_structure_;
    THROW_UNLESS(not_found, plane != nullptr);
    auto unit_z = plane->shape().unit_z();

    // orient cylinder on correct side of the plane
    auto orientation = plane->project_point(transposed).second;
    if (orientation.first < 0)
    {
        unit_z *= -1;
    }

    // TODO: scale height proportional to radius with strategy from paper
    double height = std::max(radius, 2 * (particle_radius * GfrdCfg.SINGLE_SHELL_FACTOR) + (get_distance_to_surface() + particle_radius));
    auto cylinder_center_pos = transposed - unit_z * center_particle_offset(height / 2);

    sid_pair_.second = Shell(domainID_, Cylinder(cylinder_center_pos, radius, unit_z, height/2), Shell::Code::NORMAL);
    gf_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(pid_pair_.second.D(), get_inner_a()));
    gf_iv_ = std::make_unique<GreensFunction1DRadAbs>(GreensFunction1DRadAbs(pid_pair_.second.D(), interaction_k_total(), particle_surface_dist_, get_distance_to_surface(), get_distance_to_escape(height/2)));
    return true;
}