#ifndef CONNECTIVITY_CONTAINER_HPP
#define CONNECTIVITY_CONTAINER_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "exceptions.hpp"
#include "StructureID.hpp"
#include "Vector3.hpp"
#include <array>
#include "makeString.hpp"

// --------------------------------------------------------------------------------------------------------------------------------


// implements the datastructure of objects that each have exactly 'num_neighbors' neighbors and have Tdata_ data associated
// with each edge from the one of the object to the other.
// NOTE the edges to neighbors are MANDATORY and have a specific ORDER.

// It allows you to check what structure you enter into when you leave the current one.

template<size_t NNeighbors>
class ConnectivityContainer
{
public:
   typedef std::pair<StructureID, Vector3>                                          obj_data_pair_type;

protected:
   // the structureID of the neighbor and the normal vector pointing from the
   // intersection in the direction of the center of the neighbor
   typedef std::array<obj_data_pair_type, NNeighbors>                         obj_data_pair_array_type;
   // the mapping from the structure and all its neighbors and how they are connected.
   typedef std::map<StructureID, obj_data_pair_array_type>                            obj_neighbor_objdata_array_map;

public:
   static const size_t max = NNeighbors;

   void set_neighbor_info(const StructureID obj, const size_t n, const obj_data_pair_type& obj_data_pair)
   {
      if (n >= max) throw illegal_argument(make_string() << "Index out of range for neighbor (n= " << n << "), max= " << NNeighbors);
      neighbor_mapping_[obj][n] = obj_data_pair;
   }

   obj_data_pair_type get_neighbor_info(const StructureID id, const size_t n) const
   {
      if (n >= max) throw illegal_argument(make_string() << "Index out of range for neighbor(n = " << n << "), max = " << NNeighbors);
      auto it(neighbor_mapping_.find(id));
      if (it == neighbor_mapping_.end()) throw not_found(make_string() << "Object not found(id=" << id << ")");
      return (*it).second[n];
   }

private:
   friend class Persistence;

   obj_neighbor_objdata_array_map  neighbor_mapping_;       // mapping StructureID -> array<(StructureID, vector3D), num_neighbors>
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* CONNECTIVITY_CONTAINER_HPP */

