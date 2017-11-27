#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <cassert>
#include <unordered_set>
#include "DefsEgfrd.hpp"
#include "ParticleContainer.hpp"
#include "exceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Transaction : public ParticleContainer
{
public:
   using particle_id_generator = gi::iteratorRange<std::unordered_set<ParticleID>, const ParticleID>;

   virtual particle_id_generator get_added_particles() const = 0;

   virtual particle_id_generator get_removed_particles() const = 0;

   virtual particle_id_generator get_modified_particles() const = 0;

   virtual void rollback() = 0;
};

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TParticleContainer>
class TransactionImpl : public Transaction
{
private:
   using particleid_pair_map_type = std::unordered_map<ParticleID, Particle>;
   using particleid_set_type = std::unordered_set<ParticleID>;

public:

   // matrix and size properties, route to underlying ParticleContainer

   std::size_t num_particles() const override { return pc_.num_particles(); }

   Vector3 world_size() const override { return pc_.world_size(); }

   double cell_size() const override { return pc_.cell_size(); }

   std::array<uint, 3> matrix_size() const override { return pc_.matrix_size(); }

   // particle methods

   particle_id_pair new_particle(const Particle&& p) override
   {
      particle_id_pair retval = pc_.new_particle(std::move(p));
      added_particles_.emplace(retval.first);
      return retval;
   }

   bool update_particle(const particle_id_pair& pip) override
   {
      ASSERT(removed_particles_.end() == removed_particles_.find(pip.first));

      auto r = orig_particles_.insert(particle_id_pair(pip.first, Particle()));
      if (r.second && added_particles_.end() == added_particles_.find(pip.first))
      {
         modified_particles_.emplace(pip.first);
         Particle _v(pc_.get_particle(pip.first).second);
         std::swap((*r.first).second, _v);
      }
      return pc_.update_particle(pip);
   }

   bool remove_particle(const ParticleID pid) override
   {
      auto r = orig_particles_.insert(particle_id_pair(pid, Particle()));
      if (r.second)
      {
         Particle v = pc_.get_particle(pid).second;
         std::swap((*r.first).second, v);
      }

      if (added_particles_.erase(pid) == 0)
      {
         modified_particles_.erase(pid);
         removed_particles_.emplace(pid);
      }
      else
      {
         orig_particles_.erase(pid);
      }
      return pc_.remove_particle(pid);
   }

   particle_id_pair get_particle(const ParticleID pid) const override
   {
      return pc_.get_particle(pid);
   }

   bool has_particle(const ParticleID pid) const override
   {
      return pc_.has_particle(pid);
   }

   particle_id_pair_generator get_particles() const override
   {
      return pc_.get_particles();
   }

   Transaction* create_transaction() override
   {
      return pc_.create_transaction();
   }



   //virtual particle_id_pair_and_distance_list* check_overlap(Sphere const& s) const
   //{
   //   return pc_.check_overlap(s);
   //}

   //virtual particle_id_pair_and_distance_list* check_overlap(Sphere const& s, ParticleID const& ignore) const
   //{
   //   return pc_.check_overlap(s, ignore);
   //}

   //virtual particle_id_pair_and_distance_list* check_overlap(Sphere const& s, ParticleID const& ignore1, ParticleID const& ignore2) const
   //{
   //   return pc_.check_overlap(s, ignore1, ignore2);
   //}

   //virtual structure_id_pair_and_distance_list* check_surface_overlap(Sphere const& s, Vector3 const& old_pos, StructureID const& current,
   //   double const& sigma) const
   //{
   //   return pc_.check_surface_overlap(s, old_pos, current, sigma);
   //}

   //virtual structure_id_pair_and_distance_list* check_surface_overlap(Sphere const& s, Vector3 const& old_pos, StructureID const& current,
   //   double const& sigma, StructureID const& ignore) const
   //{
   //   return pc_.check_surface_overlap(s, old_pos, current, sigma, ignore);
   //}

   //virtual structure_id_pair_and_distance_list* check_surface_overlap(Sphere const& s, Vector3 const& old_pos, StructureID const& current,
   //   double const& sigma, StructureID const& ignore1, StructureID const& ignore2) const
   //{
   //   return pc_.check_surface_overlap(s, old_pos, current, sigma, ignore1, ignore2);
   //}

   //virtual const SpeciesInfo& get_species(const SpeciesTypeID& id) const override
   //{
   //   return pc_.get_species(id);
   //}

   // StructureType stuff
//    virtual bool add_structure_type(StructureType const& Structure);   // TODO

//   virtual StructureType get_structure_type(SpeciesTypeID const& sid) const
//   {
//      return pc_.get_structure_type(sid);
//   }
//   virtual structure_types_range get_structure_types() const
//   {
//      return pc_.get_structure_types();
//   }
//   virtual SpeciesTypeID get_def_structure_type_id() const
//   {
//      return pc_.get_def_structure_type_id();
//   }
//
//   // Start Structure stuff
////    virtual StructureID add_structure(Structure const& structure);   // TODO
//
//   virtual std::shared_ptr<Structure> get_structure(StructureID const& id) const
//   {
//      return pc_.get_structure(id);
//   }
//   virtual structures_range get_structures() const
//   {
//      return pc_.get_structures();
//   }

   //virtual std::shared_ptr<Structure> get_some_structure_of_type(SpeciesTypeID const& sid) const override
   //{
   //   return pc_.get_some_structure_of_type(sid);
   //}

   //template <typename Tstructid_pair_>
   //bool update_structure(Tstructid_pair_ const& structid_pair)
   //{
   //   return pc_.update_structure(structid_pair);
   //}
   //virtual bool remove_structure(StructureID const& id)
   //{
   //   return pc_.remove_structure(id);
   //}
   //virtual structure_id_set get_structure_ids(SpeciesTypeID const& sid) const
   //{
   //   return pc_.get_structure_ids(sid);
   //}
   //virtual StructureID get_def_structure_id() const
   //{
   //   return pc_.get_def_structure_id();
   //}
   //virtual structure_id_pair_and_distance_list* get_close_structures(Vector3 const& pos, StructureID const& current_struct_id,
   //   StructureID const& ignore) const
   //{
   //   return pc_.get_close_structures(pos, current_struct_id, ignore);
   //}

   //virtual particle_id_pair_generator get_particles() const override
   //{
   //   return pc_.get_particles();
   //}


   particle_id_generator get_added_particles() const override
   {
      return particle_id_generator(added_particles_ );
   }

   particle_id_generator get_removed_particles() const override
   {
      return particle_id_generator(removed_particles_);
   }

   particle_id_generator get_modified_particles() const override
   {
      return particle_id_generator(modified_particles_);
   }

   void rollback() override
   {
      for (auto i : orig_particles_)
         pc_.update_particle(i);

      for (auto i : added_particles_)
         pc_.remove_particle(i);

      added_particles_.clear();
      modified_particles_.clear();
      removed_particles_.clear();
      orig_particles_.clear();
   }

   //virtual double distance(Vector3 const& lhs, Vector3 const& rhs) const
   //{
   //   return pc_.distance(lhs, rhs);
   //}

   //virtual Vector3 apply_boundary(Vector3 const& v) const
   //{
   //   return pc_.apply_boundary(v);
   //}

   //virtual double apply_boundary(double const& v) const
   //{
   //   return pc_.apply_boundary(v);
   //}

   //virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id) const
   //{
   //   return pc_.apply_boundary(pos_struct_id);
   //}

   //virtual Vector3 cyclic_transpose(Vector3 const& p0, Vector3 const& p1) const
   //{
   //   return pc_.cyclic_transpose(p0, p1);
   //}

   //virtual double cyclic_transpose(double const& p0, double const& p1) const
   //{
   //   return pc_.cyclic_transpose(p0, p1);
   //}

   //virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
   //   Structure const& structure) const
   //{
   //   return pc_.cyclic_transpose(pos_struct_id, structure);
   //}

   TransactionImpl(TParticleContainer& pc) : pc_(pc) { }

private:
   particle_id_pair get_original_particle(const ParticleID id) const
   {
      auto i(orig_particles_.find(id));
      if (orig_particles_.end() == i) throw not_found(make_string() << "No such particle: id=" << id);
      return *i;
   }

private:
   TParticleContainer&         pc_;                    // the associated particle container
   particleid_set_type         added_particles_;
   particleid_set_type         modified_particles_;
   particleid_set_type         removed_particles_;
   particleid_pair_map_type    orig_particles_;
};

// --------------------------------------------------------------------------------------------------------------------------------
