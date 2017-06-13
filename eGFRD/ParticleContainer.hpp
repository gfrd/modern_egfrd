#ifndef PARTICLE_CONTAINER_HPP
#define PARTICLE_CONTAINER_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <utility>
#include <vector>
#include "Particle.hpp"
#include "ParticleID.hpp"
#include "genericIterator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Transaction;

class ParticleContainer
{
public:

   using particle_id_pair = Particle::particle_id_pair;
   using particle_id_pair_and_distance = std::pair<particle_id_pair, double>;
   using particle_id_pair_and_distance_list = std::vector<particle_id_pair_and_distance>;
   using particle_id_pair_generator = agi::iteratorRange<particle_id_pair>;
   using particle_id_list = const std::vector<ParticleID>&;

   virtual ~ParticleContainer() = default;

   // matrix and size properties

   virtual size_t num_particles() const = 0;

   virtual Vector3 world_size() const = 0;

   double world_volume() const { auto ws = world_size();  return ws.X() * ws.Y() * ws.Z(); }

   virtual double cell_size() const = 0;

   virtual std::array<uint, 3> matrix_size() const = 0;


   // particle methods

   virtual particle_id_pair new_particle(const Particle&& p) = 0;

   virtual bool update_particle(const particle_id_pair& pip) = 0;

   virtual bool remove_particle(const ParticleID pid) = 0;

   virtual particle_id_pair get_particle(const ParticleID pid) const = 0;

   virtual bool has_particle(const ParticleID pid) const = 0;

   virtual particle_id_pair_generator get_particles() const = 0;

   // create transaction object (to undo changes in particlecontainer)     /* NOT USED */

   virtual Transaction* create_transaction() = 0;




   // not migrated!!

   // virtual const SpeciesInfo& get_species(const SpeciesTypeID& id) const = 0;                                                       NOT VIRTUAL, because species are not is ParticleContainer

   //virtual bool add_structure_type(const StructureType& Structure) = 0;   // TODO


   // structure methods

   //virtual bool has_structure(const StructureID& id) const = 0;

   //virtual std::shared_ptr<Structure> get_structure(const StructureID& id) const = 0;

   //virtual std::shared_ptr<Structure> get_some_structure_of_type(const SpeciesTypeID& sid) const = 0;

   //virtual bool remove_structure(const StructureID& id) = 0;

   // virtual StructureID get_def_structure_id() const = 0;

   // virtual bool update_structure(const structure_id_pair& structid_pair) = 0; <- not virtual because of template calls to StructureContianer

   //virtual structure_id_set get_structure_ids(const SpeciesTypeID& sid) const = 0;

   //virtual structures_range get_structures() const = 0;     // TODO? or python binding only function?

   //virtual StructureType get_structure_type(const SpeciesTypeID& sid) const = 0; // TODO? or python binding only function?

   // virtual structure_types_range get_structure_types() const = 0; 

   // virtual SpeciesTypeID get_def_structure_type_id() const = 0;

   // virtual StructureID add_structure(const Structure& structure) = 0;   // TODO


   //virtual structure_id_pair_and_distance_list* get_close_structures(const Vector3& pos, const StructureID& current_struct_id, const StructureID& ignore) const = 0;

   //virtual structure_id_pair_and_distance_list* check_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID& current, double sigma) const = 0;

   //virtual structure_id_pair_and_distance_list* check_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID& current, double sigma, const StructureID& ignore) const = 0;

   //virtual structure_id_pair_and_distance_list* check_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID& current, double sigma, const StructureID& ignore1, const StructureID& ignore2) const = 0;

   //virtual particle_id_pair_and_distance_list check_overlap(const Sphere& s) const = 0;

   //virtual particle_id_pair_and_distance_list* check_overlap(const Sphere& s, const ParticleID& ignore) const = 0;

   //virtual particle_id_pair_and_distance_list* check_overlap(const Sphere& s, const ParticleID& ignore1, const ParticleID& ignore2) const = 0;

   //virtual double distance(const Vector3& lhs, const Vector3& rhs) const = 0;

   //virtual Vector3 apply_boundary(const Vector3& v) const = 0;

   //virtual double apply_boundary(double v) const = 0;

   //virtual position_structid_pair apply_boundary(const position_structid_pair& pos_struct_id) const = 0;

   //virtual Vector3 cyclic_transpose(const Vector3& p0, const Vector3& p1) const = 0;

   //virtual double cyclic_transpose(double p0, double p1) const = 0;

   //virtual position_structid_pair cyclic_transpose(const position_structid_pair& pos_struct_id, const Structure& structure) const = 0;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* PARTICLE_CONTAINER_HPP */
