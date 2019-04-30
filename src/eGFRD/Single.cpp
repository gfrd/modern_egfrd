#include "Single.hpp"
#include "exceptions.hpp"
#include "ShellCreateUtils.hpp"
#include "EGFRDSimulator.hpp"
#include <GreensFunction2DAbsSym.hpp>
#include <GreensFunction3DAbsSym.hpp>

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
         double distance = s->distance(pos);    // fix cyclic world !
         radius = std::min(radius, distance);
      }
   }

   shell_distance_checker sdc(shell_id(), pos, radius, shell_distance_checker::Construct::SHARE5050);
   CompileConfigSimulator::TBoundCondition::each_neighbor(smat, sdc, pos);
   radius = std::min(radius, sdc.distance());

   if (radius < min_radius) return false;             // no space for Single domain.. it will be discarded, and a Multi is created instead!

   sid_pair_.second = Shell(domainID_, Sphere(pos, radius), Shell::Code::NORMAL);
   gf_ = std::make_unique<GreensFunction3DAbsSym>(GreensFunction3DAbsSym(pid_pair_.second.D(), get_inner_a()));
   return true;
}

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT bool SingleCylindrical::create_updated_shell(const shell_matrix_type& smat, const World& world)
{
   auto pos = pid_pair_.second.position();
   double min_radius = particle().radius() * GfrdCfg.MULTI_SHELL_FACTOR;
   double max_radius = smat.cell_size() / std::sqrt(8.0);         // any angle cylinder must fit into cell-matrix! 2*sqrt(2)
   //double max_heigth = smat.cell_size() / std::sqrt(2.0);
   auto radius = max_radius;

   // check distances to surfaces, ignore def.struct and particle.structure
   for (auto s : world.get_structures())
   {
      if (s->id() != pid_pair_.second.structure_id() && s->id() != world.get_def_structure_id())    // ignore structure that particle is attached to
      {
         double distance = s->distance(pos);    // fix cyclic world !
         radius = std::min(radius, distance);
      }
   }
   // TODO 
   shell_distance_checker sdc(shell_id(), pos, radius, shell_distance_checker::Construct::SHARE5050);
   CompileConfigSimulator::TBoundCondition::each_neighbor(smat, sdc, pos);
   radius = std::min(radius, sdc.distance()) / GfrdCfg.SAFETY;

   if (radius < min_radius) return false;             // no space for Single domain.. it will be discarded, and a Multi is created instead!

   auto structure = world.get_structure(pid_pair_.second.structure_id());
   auto *plane = dynamic_cast<PlanarSurface*>(structure.get());
   THROW_UNLESS(not_found, plane != nullptr);
   auto unit_z = plane->shape().unit_z();
   structure_ = std::dynamic_pointer_cast<PlanarSurface>(structure);

   sid_pair_.second = Shell(domainID_, Cylinder(pos, radius, unit_z, radius), Shell::Code::NORMAL);
   gf_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(pid_pair_.second.D(), get_inner_a()));
   return true;
}

// --------------------------------------------------------------------------------------------------------------------------------


GFRD_EXPORT bool SinglePlanarInteraction::create_updated_shell(const shell_matrix_type& smat, const World& world)
{
    auto pos = pid_pair_.second.position();
    double min_radius = particle().radius() * GfrdCfg.MULTI_SHELL_FACTOR;
    double max_radius = smat.cell_size() / std::sqrt(8.0);         // any angle cylinder must fit into cell-matrix! 2*sqrt(2)
    //double max_heigth = smat.cell_size() / std::sqrt(2.0);
    auto radius = max_radius;

    // check distances to surfaces, ignore def.struct and particle.structure
    for (auto s : world.get_structures())
    {
        if (s->id() != interacting_structure_.id() && s->id() != world.get_def_structure_id())    // ignore structure that particle is attached to
        {
            double distance = s->distance(pos);    // fix cyclic world !
            radius = std::min(radius, distance);
        }
    }
    // TODO
    shell_distance_checker sdc(shell_id(),  pos, radius, shell_distance_checker::Construct::SHARE5050);
    CompileConfigSimulator::TBoundCondition::each_neighbor(smat, sdc, pos);
    radius = std::min(radius, sdc.distance()) / GfrdCfg.SAFETY;

    if (radius < min_radius) return false;             // no space for Single domain.. it will be discarded, and a Multi is created instead!

    auto plane = &interacting_structure_;
    THROW_UNLESS(not_found, plane != nullptr);
    auto unit_z = plane->shape().unit_z();

    sid_pair_.second = Shell(domainID_, Cylinder(pos, radius, unit_z, radius), Shell::Code::NORMAL);
    gf_ = std::make_unique<GreensFunction2DAbsSym>(GreensFunction2DAbsSym(pid_pair_.second.D(), get_inner_a()));
    return true;
}