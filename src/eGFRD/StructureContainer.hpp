#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <set>
#include <map>
#include <memory>
#include "exceptions.hpp"
#include "ConnectivityContainer.hpp"
#include "CuboidalRegion.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "DiskSurface.hpp"
#include "SphericalSurface.hpp"
#include "makeString.hpp"
#include "SerialIDGenerator.hpp"
#include "genericIterator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class StructureContainer
{
public:
   typedef std::set<StructureID>                                      structure_id_set;
   typedef std::pair<const StructureID, std::shared_ptr<Structure>>  structure_id_pair;
   typedef std::map<StructureID, std::shared_ptr<Structure>>         structure_map;

protected:
   typedef std::map<StructureID, structure_id_set>                       per_structure_substructure_id_set;
   typedef std::pair<StructureID, Vector3>                               neighbor_id_vector; // FIXME this is actually a datatype of the ConnectivityContainer.

public:
   typedef std::pair<StructureID, std::shared_ptr<SphericalSurface>>    spherical_surface_id_pair_type;
   typedef std::pair<StructureID, std::shared_ptr<DiskSurface>>         disk_surface_id_pair_type;
   typedef std::pair<StructureID, std::shared_ptr<CylindricalSurface>>  cylindrical_surface_id_pair_type;
   typedef std::pair<StructureID, std::shared_ptr<PlanarSurface>>       planar_surface_id_pair_type;
   typedef std::pair<StructureID, std::shared_ptr<CuboidalRegion>>      cuboidal_region_id_pair_type;
   typedef ConnectivityContainer<2>                                      cylindrical_surface_bc_type;  // FIXME number of neighbors should be the length of the side enumerator
   typedef ConnectivityContainer<4>                                      planar_surface_bc_type;
   typedef ConnectivityContainer<6>                                      cuboidal_region_bc_type;
   typedef std::pair<Vector3, StructureID>                               position_structid_pair;
   typedef gi::iteratorRange<structure_map, std::shared_ptr<Structure>, gi::SelectSecond<std::shared_ptr<Structure>>>                 structures_range;

   StructureContainer() noexcept { }

   template <typename TStructure>
   StructureID add_structure(std::shared_ptr<TStructure> structure)
   {
      const StructureID sid(structidgen_());
      structure->set_id(sid);
      VERIFY(update_structure(std::make_pair(sid, structure)));
      return sid;
   }

   // Update_structure methods for the different kind of structures supported by the StructureContainer
   bool update_structure(const spherical_surface_id_pair_type& structid_sphere)
   {
      return update_structure_base(structid_sphere);  // No connecting needs to be done
   }

   bool update_structure(const disk_surface_id_pair_type& structid_disk)
   {
      return update_structure_base(structid_disk);
   }

   bool update_structure(const cuboidal_region_id_pair_type& structid_cube)
   {
      if (update_structure_base(structid_cube))
      {
         // add to Connectivity container for cuboidal regions
         // NOTE this information is never queried!! Maybe leave the neighbor info out?
         //cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 0, std::make_pair(structid_cube.first, structid_cube.second->shape().unit_z()));
         //cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 1, std::make_pair(structid_cube.first, -structid_cube.second->shape().unit_z()));
         //cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 2, std::make_pair(structid_cube.first, structid_cube.second->shape().unit_y()));
         //cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 3, std::make_pair(structid_cube.first, -structid_cube.second->shape().unit_y()));
         //cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 4, std::make_pair(structid_cube.first, structid_cube.second->shape().unit_x()));
         //cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 5, std::make_pair(structid_cube.first, -structid_cube.second->shape().unit_x()));
         return true;
      }
      return false;
   }

   bool update_structure(const cylindrical_surface_id_pair_type& structid_cylinder)
   {
      if (update_structure_base(structid_cylinder))
      {
         // add to Connectivity container for cylindrical surfaces
         //cylindrical_structs_bc_.set_neighbor_info(structid_cylinder.first, 0, std::make_pair(structid_cylinder.first, -structid_cylinder.second->shape().unit_z()));
         //cylindrical_structs_bc_.set_neighbor_info(structid_cylinder.first, 1, std::make_pair(structid_cylinder.first, structid_cylinder.second->shape().unit_z()));
         return true;
      }
      return false;
   }

   bool update_structure(const planar_surface_id_pair_type& structid_plane)
   {
      // call regular update method (this also checks if the structure was already present)
      if (update_structure_base(structid_plane))
      {
         // We now assume that the structure was not already in the boundary_conditions_thing
         // add to Connectivity container for planar surfaces and set *reflective* boundary conditions. // TODO This may cause many problems in single reactions!!! 
         planar_structs_bc_.set_neighbor_info(structid_plane.first, 0, std::make_pair(structid_plane.first, -structid_plane.second->shape().unit_y()));
         planar_structs_bc_.set_neighbor_info(structid_plane.first, 1, std::make_pair(structid_plane.first, structid_plane.second->shape().unit_y()));
         planar_structs_bc_.set_neighbor_info(structid_plane.first, 2, std::make_pair(structid_plane.first, structid_plane.second->shape().unit_x()));
         planar_structs_bc_.set_neighbor_info(structid_plane.first, 3, std::make_pair(structid_plane.first, -structid_plane.second->shape().unit_x()));
         return true;
      }
      return false;
   }

   // Methods to get the 'n-th' neighbor structure of a particular structure from the container.
   // We pass the structure and not just the ID because with just the ID we can't differentiate between the various types of structures.
   neighbor_id_vector get_neighbor_info(const PlanarSurface& structure, int n) const
   {
      return planar_structs_bc_.get_neighbor_info(structure.id(), n);
   }

   //neighbor_id_vector get_neighbor_info(const CylindricalSurface& structure, int n) const
   //{
   //   return cylindrical_structs_bc_.get_neighbor_info(structure.id(), n);
   //}

   ////bool connect_structures(PlanarSurface const& structure1, planar_surface_side_type const& side1, PlanarSurface const& structure2, planar_surface_side_type const& side2)
   //bool connect_structures(const PlanarSurface& structure1, int side1, const PlanarSurface& structure2, int side2)
   //{
   //   // TODO if the structure are not the same, check that the planes actually touch at the right edge.

   //   // First check whether the planes are orthogonal or parallel, or none of this.
   //   bool planes_are_orthogonal(false);
   //   bool planes_are_parallel(false);
   //   Vector3 structure1_unit_z(structure1.shape().unit_z());
   //   Vector3 structure2_unit_z(structure2.shape().unit_z());

   //   if (feq(Vector3::dot(structure1_unit_z, structure2_unit_z), 0.0, 1.0))
   //      planes_are_orthogonal = true;
   //   if (feq(Vector3::dot(structure1_unit_z, structure2_unit_z), 1.0, 1.0))
   //      planes_are_parallel = true;
   //   // If they are neither parallel or orthogonal we cannot treat this situation and stop here
   //   if (!(planes_are_orthogonal || planes_are_parallel))
   //      throw unsupported("StructureContainer: only orthogonal or parallel planes can be connected.");

   //   // Now determine which of the unit vectors of the neighboring plane is the one that shall be
   //   // used for transforming the position that reaches out of the origin plane. This is relevant
   //   // only for the orthogonal case. If the planes are parallel, the position does not have
   //   // to be transformed at all (only the StructureID changes).
   //   // For parallel planes we therefore return a zero vector below. This at once serves as 
   //   // a flag indicating that the planes are indeed parallel.
   //   //Vector3 zero_vector( zero_vector() );
   //   Vector3 structure1_vector(Vector3(0.0, 0.0, 0.0));
   //   Vector3 structure2_vector(Vector3(0.0, 0.0, 0.0));
   //   // this is already fine for parallel planes

   //   if (planes_are_orthogonal)
   //   {
   //      // Side convention for coding neighbors:
   //      //
   //      //  |-------------|
   //      //  |      0      |
   //      //  |             |
   //      //  | 2         3 |
   //      //  |             |
   //      //  |      1      |
   //      //  |-------------|
   //      //
   //      // FIXME cleanup code below            
   //      switch (side1)
   //      {
   //      case 0: structure1_vector = -structure1.shape().unit_y();
   //         break;
   //      case 1: structure1_vector = structure1.shape().unit_y();
   //         break;
   //      case 2: structure1_vector = structure1.shape().unit_x();
   //         break;
   //      case 3: structure1_vector = -structure1.shape().unit_x();
   //         break;
   //      }
   //      switch (side2)
   //      {
   //      case 0: structure2_vector = -structure2.shape().unit_y();
   //         break;
   //      case 1: structure2_vector = structure2.shape().unit_y();
   //         break;
   //      case 2: structure2_vector = structure2.shape().unit_x();
   //         break;
   //      case 3: structure2_vector = -structure2.shape().unit_x();
   //         break;
   //      }
   //   }

   //   planar_structs_bc_.set_neighbor_info(structure1.id(), side1, std::make_pair(structure2.id(), structure2_vector));
   //   planar_structs_bc_.set_neighbor_info(structure2.id(), side2, std::make_pair(structure1.id(), structure1_vector));

   //   return true;
   //}

   ////bool connect_structures(CylindricalSurface const& structure1, cylindrical_surface_side_type const& side1,CylindricalSurface const& structure2, cylindrical_surface_side_type const& side2)

   //bool connect_structures(const CylindricalSurface& structure1, int side1, const CylindricalSurface& structure2, int side2)
   //{
   //   if (structure1.id() != structure2.id())
   //   {
   //      throw unsupported("Two different cylinders cannot be connected (yet...).");
   //   }
   //   const double factor((side1 == 1) ? 1.0 : -1.0);
   //   if (side1 == side2)
   //   {
   //      cylindrical_structs_bc_.set_neighbor_info(structure1.id(), side1, std::make_pair(structure2.id(), structure2.shape().unit_z() * -factor));
   //      cylindrical_structs_bc_.set_neighbor_info(structure2.id(), side2, std::make_pair(structure1.id(), structure1.shape().unit_z() * -factor));
   //   }
   //   else
   //   {
   //      cylindrical_structs_bc_.set_neighbor_info(structure1.id(), side1, std::make_pair(structure2.id(), structure2.shape().unit_z() *  factor));
   //      cylindrical_structs_bc_.set_neighbor_info(structure2.id(), side2, std::make_pair(structure1.id(), structure1.shape().unit_z() * -factor));
   //   }
   //   return true;
   //}

   //// Template member function to call the right 'apply_boundary' function for the various types of structures.
   //template <typename Tstructure_>
   //position_structid_pair apply_boundary(const Tstructure_& structure, const position_structid_pair& pos_struct_id) const
   //{
   //   return ::apply_boundary(pos_struct_id, structure, *this);
   //}
   //template <typename Tstructure_>
   //position_structid_pair cyclic_transpose(const Tstructure_& structure, const position_structid_pair& pos_struct_id) const
   //{
   //   return ::cyclic_transpose(pos_struct_id, structure, *this);
   //}


   // --------------------------------------------------------------------------------------------------------------------------------

   bool has_structure(const StructureID id) const
   {
      structure_map::const_iterator i(structure_map_.find(id));
      return i != structure_map_.end();
   }

   std::shared_ptr<Structure> get_structure(const StructureID id) const
   {
      structure_map::const_iterator i(structure_map_.find(id));
      if (structure_map_.end() == i) throw not_found(make_string() << "Unknown structure (id=" << id << ")");
      return (*i).second;
   }

   structures_range get_structures() const
   {
      return gi::iteratorRange<structure_map, std::shared_ptr<Structure>, gi::SelectSecond<std::shared_ptr<Structure>>>(structure_map_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   //std::shared_ptr<Structure> get_some_structure_of_type(const StructureTypeID sid) const
   //{
   //   // This returns a structure of the structure type given as an argument, if it exists
   //   // in the container. Otherwise it throws an exception.
   //   structure_map::const_iterator i;

   //   for (i = structure_map_.begin(); i != structure_map_.end(); ++i)

   //      if ((*i).second->sid() == sid)        break;

   //   if (structure_map_.end() == i)
   //      throw not_found(make_string() << "Could not find any structure of specified structure type (sid=" << sid << ")");
   //   return (*i).second;
   //}

   //bool remove_structure(const StructureID id)
   //{
   //   // TODO
   //   //  -find StructureID in all ConnectivityContainers -> remove references
   //   //  -get all the substructures of structure
   //   //  -only remove if no substructures
   //   return false;
   //}

   structure_id_set get_visible_structures(const StructureID id) const
   {
      // This function returns all structures that are visible form the structure with the ID.
      // A structure is visible if it is a substructure or if it has the same parent structure.
      auto structure = get_structure(id);

      // First check the "sibling" structures
      structure_id_set visible_structures;
      if (structure->parent_id() == id || !structure->parent_id())      // self-ref or zero-id (means top)
      {
         visible_structures = structure_id_set();
      }
      else
      {
         visible_structures = get_visible_structures(structure->parent_id());
         visible_structures.erase(id);
      }

      // Now the substructures
      structure_id_set substructures = get_substructure_ids(id);
      visible_structures.insert(substructures.begin(), substructures.end());
      return visible_structures;
   }

private:

   bool update_structure_base(const structure_id_pair& structid_pair)
   {
      structure_map::const_iterator i(structure_map_.find(structid_pair.first));
      if (i != structure_map_.end())
      {
         if ((*i).second->parent_id() != structid_pair.second->parent_id())
            // If the structure had a different parent structure
            // Note that all the substructures of the structures are kept (they are also moved moved through the hierarchy)
         {
            structure_substructures_map_[(*i).second->parent_id()].erase((*i).first);
            structure_substructures_map_[structid_pair.second->parent_id()].insert(structid_pair.first);
         }
         structure_map_[(*i).first] = structid_pair.second;
         return false;
      }

      // The structure was not yet in the world. Create a new item in the structure mapping
      structure_map_[structid_pair.first] = structid_pair.second;
      // add the id the mapping 'super_structure_id -> set of substructures'
      structure_substructures_map_[structid_pair.second->parent_id()].insert(structid_pair.first);
      // create a new mapping from structure id -> set of substructures
      structure_substructures_map_[structid_pair.first] = structure_id_set();
      return true;
   }

   structure_id_set get_substructure_ids(const StructureID id) const
   {
      per_structure_substructure_id_set::const_iterator i(structure_substructures_map_.find(id));
      if (i == structure_substructures_map_.end()) throw not_found(make_string() << "Unknown structure (id=" << id << ")");
      return (*i).second;
   }

protected:
   friend class Persistence;

   SerialIDGenerator<StructureID>      structidgen_;       // structureID generator
   structure_map                       structure_map_;     // mapping: structure_id -> structure
   per_structure_substructure_id_set   structure_substructures_map_;

   // The next objects describe the boundary conditions that are applied to the different types of structures.
   planar_surface_bc_type      planar_structs_bc_;

   //cylindrical_surface_bc_type cylindrical_structs_bc_;
   //cuboidal_region_bc_type     cuboidal_structs_bc_;
};

/*

inline std::pair<Vector3, StructureID> cyclic_transpose(const std::pair<Vector3, StructureID>& pos_structure_id, const PlanarSurface& planar_surface, const StructureContainer& sc)
{
}


// --------------------------------------------------------------------------------------------------------------------------------

// The template needs to be parameterized with the appropriate shape (which then parameterizes the type
// of the ConnectivityContainer that we use.
// We supply the structure container as an argument so that we can get the structures that we need and
// to query the boundary conditions and connectivity.

inline std::pair<Vector3, StructureID> apply_boundary(const std::pair<Vector3, StructureID>& pos_structure_id, const CylindricalSurface& cylindrical_surface, const StructureContainer& sc)
{
   typedef CylindricalSurface::shape_type   cylinder_type;
   typedef std::pair<StructureID, Vector3>         neighbor_id_vector;


   // Note that we assume that the new position is in the cylinder (dot(pos, unit_r)==0)
   const cylinder_type origin_cylinder(cylindrical_surface.shape());
   const double half_length(origin_cylinder.half_length());

   const Vector3 pos_vector((pos_structure_id.first - origin_cylinder.position()));
   const double component_z(Vector3::dot(pos_vector, origin_cylinder.unit_z()));

   // declare the variables that will be written
   StructureID new_id(pos_structure_id.second);
   Vector3 neighbor_cylinder_inl;

   if (half_length < component_z)
   {
      // we are at the 'right' side of the cylinder (side nr. 0, side in the direction of the normal vector)
      // get the unit vector pointing from the edge between the two planes to the center of the plane
      const neighbor_id_vector niv(sc.get_neighbor_info(cylindrical_surface, 0));

      new_id = niv.first;
      neighbor_cylinder_inl = origin_cylinder.unit_z()* half_length + niv.second* (component_z - half_length);
   }
   else if (component_z < -half_length)
   {
      // we are at the 'right' side of the cylinder (side nr. 0, side in the direction of the normal vector)
      // get the unit vector pointing from the edge between the two planes to the center of the plane
      const neighbor_id_vector niv(sc.get_neighbor_info(cylindrical_surface, 1));

      new_id = niv.first;
      neighbor_cylinder_inl = origin_cylinder.unit_z() * -half_length + niv.second *  (-component_z - half_length);
   }
   else
   {
      // we are still in the cylinder (did not pass any of the boundaries)
      // don't have to do anything
      return pos_structure_id;
   }

   const Vector3 new_pos(origin_cylinder.position() + neighbor_cylinder_inl);
   return std::make_pair(new_pos, new_id);
}

// --------------------------------------------------------------------------------------------------------------------------------

inline std::pair<Vector3, StructureID> cyclic_transpose(const std::pair<Vector3, StructureID>& pos_structure_id, const CylindricalSurface& cylindrical_surface, const StructureContainer& sc)
{
   // This is not used and not implemented yet.
}
*/
// --------------------------------------------------------------------------------------------------------------------------------
