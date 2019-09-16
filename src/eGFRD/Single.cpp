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
    auto particle = pid_pair_.second;
    auto structure_id = particle.structure_id();
    auto pos = particle.position();
    double min_radius = particle.radius() * GfrdCfg.MULTI_SHELL_FACTOR;
    double max_radius = std::min(smat.cell_size() / std::sqrt(8.0),                       // any angle cylinder must fit into cell-matrix! 2*sqrt(2)
                                 scaling::dist_to_plane_edge(pos, structure_id, world));  // and not exceed its plane edges
    auto height = 2 * min_radius;

    std::vector<ShellID> ignored_shells = {shell_id()};
    std::vector<StructureID> ignored_structures = {structure_id, world.get_def_structure_id()};

    auto structure = world.get_structure(structure_id);
    auto *plane = dynamic_cast<PlanarSurface*>(structure.get());
    THROW_UNLESS(not_found, plane != nullptr);
    auto unit_z = plane->shape().unit_z();
    structure_ = std::dynamic_pointer_cast<PlanarSurface>(structure);

    auto start1 = pos - unit_z * height / 2;
    auto end1 = pos + unit_z * height / 2;
    auto radius = std::min(max_radius, scaling::find_maximal_cylinder_radius(start1, end1, smat, world,
            ignored_structures, ignored_shells));

    radius /= GfrdCfg.SAFETY;

    if (radius < min_radius)
    {
        return false;             // no space for Single domain.. it will be discarded, and a Multi is created instead!
    }


    sid_pair_.second = Shell(domainID_, Cylinder(pos, radius, unit_z, height/2), Shell::Code::NORMAL);
    gf_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(pid_pair_.second.D(), get_inner_a()));
    return true;
}

// --------------------------------------------------------------------------------------------------------------------------------


GFRD_EXPORT bool SinglePlanarInteraction::create_updated_shell(const shell_matrix_type& smat, const World& world)
{
    auto particle = pid_pair_.second;
    auto structure_id = particle.structure_id();
    auto pos = particle.position();
    auto transposed = world.cyclic_transpose(pos, interacting_structure_.position());
    auto particle_radius = particle.radius();
    particle_surface_dist_ = interacting_structure_.distance(transposed) - particle_radius;

    std::vector<ShellID> ignored_shells = {shell_id()};
    std::vector<StructureID> ignored_structures = {interacting_structure_.id(), particle.structure_id(), world.get_def_structure_id()};

    if(particle_surface_dist_ < 0.0) {
        Log("GFRD").warn() << "Particle " << pid_pair_.first << " is touching a surface it is not bound to, BD motion might have been erroneous.";
        return false;
    }

    double min_radius = particle.radius() * GfrdCfg.SINGLE_SHELL_FACTOR;
    double max_radius = std::min(smat.cell_size() / std::sqrt(8.0),                       // any angle cylinder must fit into cell-matrix! 2*sqrt(2)
                                 scaling::dist_to_plane_edge(pos, interacting_structure_.id(), world));  // and not exceed its plane edges
    double max_height = smat.cell_size() / std::sqrt(8.0);                                // any angle cylinder must fit into cell-matrix! 2*sqrt(2)

    auto plane = &interacting_structure_;
    THROW_UNLESS(not_found, plane != nullptr);
    auto unit_z = plane->shape().unit_z();

    // Orient cylinder on correct side of the plane
    auto orientation = plane->project_point(transposed).second;
    if (orientation.first < 0)
    {
        unit_z *= -1;
    }

    auto height_through_surface = particle_radius * GfrdCfg.SINGLE_SHELL_FACTOR;
    auto height_to_surface = get_distance_to_surface();
    auto static_height = height_through_surface + height_to_surface + particle.radius();
    auto base_pos = transposed - static_height * unit_z;
    auto scaling_factor = 1;
    auto max_dynamic_height = scaling::find_maximal_cylinder_height<shell_matrix_type>(base_pos, unit_z, static_height,
            scaling_factor, smat, world, ignored_structures, ignored_shells);

    // Ensure we don't exceed the matrix cell dimensions
    max_dynamic_height = std::min(max_dynamic_height, max_height - static_height);

    auto radius = max_dynamic_height * scaling_factor;

    // Ensure we don't exceed the matrix cell dimensions or the plane edges
    radius = std::min(radius, max_radius);

    if (radius < min_radius)
    {
        return false;
    }

    if (radius != max_dynamic_height) {
        // Radius was changed to prevent exceeding cell dimensions or plane edge,
        // so we must change the height along with it to ensure the correct scaling factor
        max_dynamic_height = radius / scaling_factor;
    }

    auto height = static_height + max_dynamic_height;
    height /= GfrdCfg.SAFETY;

    // Cylinder height should fit at least the portion through the surface, the distance from particle to surface,
    // and the particle radius.
    if (height < (particle_radius * GfrdCfg.SINGLE_SHELL_FACTOR) + height_to_surface + height_through_surface)
    {
        return false;
    }

    auto cylinder_center_pos = transposed - unit_z * center_particle_offset(height / 2);

    sid_pair_.second = Shell(domainID_, Cylinder(cylinder_center_pos, radius, unit_z, height/2), Shell::Code::NORMAL);
    gf_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(pid_pair_.second.D(), get_inner_a()));
    gf_iv_ = std::make_unique<GreensFunction1DRadAbs>(GreensFunction1DRadAbs(pid_pair_.second.D(), interaction_k_total(), particle_surface_dist_, get_distance_to_surface(), get_distance_to_escape(height/2)));
    return true;
}