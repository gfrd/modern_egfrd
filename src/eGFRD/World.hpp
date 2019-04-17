#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <map>
#include "helperFunctions.hpp"
#include "genericIterator.hpp"
#include "ParticleContainerImpl.hpp"
#include "exceptions.hpp"
#include "Particle.hpp"
#include "ParticleID.hpp"
#include "SpeciesTypeID.hpp"
#include "Model.hpp"
#include "Structure.hpp"
#include "StructureType.hpp"
#include "StructureID.hpp"
#include "SpeciesType.hpp"
#include "Logger.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class World : public ParticleContainerImpl
{
public:
   using base_type = ParticleContainerImpl;
   using boundary_type = base_type::boundary_type;

protected:
   using particle_id_set = std::set<ParticleID>;
   using structure_id_set = std::set<StructureID>;
   using species_map = std::map<SpeciesTypeID, SpeciesType>;
   using structure_map = std::map<StructureID, std::shared_ptr<Structure>>;
   using structure_type_map = std::map<StructureTypeID, StructureType>;
   using per_species_particle_id_set = std::map<SpeciesTypeID, particle_id_set>;
   using per_structure_type_structure_id_set = std::map<StructureTypeID, structure_id_set>;
   using per_structure_particle_id_set = std::map<StructureID, particle_id_set>;

public:
   using SpeciesIteratorRange = gi::iteratorRange<species_map>;
   using StructureIteratorRange = gi::iteratorRange<structure_map>;
   using position_structid_pair = std::pair<Vector3, StructureID>;
   using structure_id_pair = std::pair<StructureID, std::shared_ptr<Structure>>;
   using structure_id_pair_and_distance = std::pair<structure_id_pair, double>;
   using structure_id_pair_and_distance_list = std::vector<structure_id_pair_and_distance>;


   World() : ParticleContainerImpl(), default_structure_type_id_(0), default_structure_id_(0) { }

   particle_id_pair add_particle(const SpeciesTypeID sid, const StructureID structure_id, const Vector3& pos)
   {
      return new_particle(Particle(get_species(sid), structure_id, pos));        // lookup speciesTypeID and construct Particle
   }

   particle_id_pair new_particle(const Particle&& p) override
   {
      // check speciesType matches structureType
      ASSERT(structures_.get_structure(p.structure_id())->sid() == get_species(p.sid()).structure_type_id());

      // insert into particle container, indirectly calls update_particle
      return base_type::new_particle(std::move(p));
   }


   bool update_particle(const particle_id_pair& pi_pair) override
   {
      particle_matrix_type::iterator i(pmat_.find(pi_pair.first));
      if (i != pmat_.end())       // The particle was already in the matrix
      {
         if ((*i).second.sid() != pi_pair.second.sid())                          // If the species changed we need to change the 'species_id -> particle ids' mapping
         {
            particle_pool_[(*i).second.sid()].erase((*i).first);
            particle_pool_[pi_pair.second.sid()].insert(pi_pair.first);
         }
         if ((*i).second.structure_id() != pi_pair.second.structure_id())        // If the structure changed we need to change the 'structure_id -> particle ids' mapping
         {
            particleonstruct_pool_[(*i).second.structure_id()].erase((*i).first);
            particleonstruct_pool_[pi_pair.second.structure_id()].insert(pi_pair.first);
         }
         pmat_.update(i, pi_pair);
         return false;
      }

      // The particle didn't exist yet
      VERIFY(base_type::update_particle(pi_pair));                                     // update particle
      particle_pool_[pi_pair.second.sid()].insert(pi_pair.first);                      // update species->particles map
      particleonstruct_pool_[pi_pair.second.structure_id()].insert(pi_pair.first);     // update structure->particles map
      return true;
   }

   virtual bool remove_particle(const ParticleID id) override
   {
      bool found(false);
      particle_id_pair pp(get_particle(id, found));                                     // call method in ParticleContainerBase
      if (!found) return false;
      particle_pool_[pp.second.sid()].erase(id);                                        // remove particle from species->particles
      particleonstruct_pool_[pp.second.structure_id()].erase(id);                       // remove particle from structure->particles
      base_type::remove_particle(id);                                                   // remove particle
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   position_structid_pair apply_boundary(const position_structid_pair& pos_struct_id) const
   {
      //         // cyclic transpose the position with the structure
      //         const std::shared_ptr<const Structure> structure (get_structure(pos_struct_id.second));
      //         const position_structid_pair cyc_pos_struct_id (std::make_pair(cyclic_transpose(pos_struct_id.first, structure->position()),pos_struct_id.second));

      // TESTING: cyclic transpose should not be applied when going around an edge!
      // TODO Clean this up when it is clear that cyclic_transpose is unnecessary here.

      const std::shared_ptr<Structure> structure(structures_.get_structure(pos_struct_id.second));
      // 1. applies the boundary of the structure/surface (go around the corner etc)
      auto new_pos_struct_id(structure->apply_boundary(pos_struct_id, structures_));
      // 2. Then also apply the boundary condition of the world
      return std::make_pair(apply_boundary(new_pos_struct_id.first), new_pos_struct_id.second);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double distance(const Vector3& lhs, const Vector3& rhs) const
   {
      return boundary_type::distance(lhs, rhs, world_size());
   }

   Vector3 apply_boundary(const Vector3& v) const
   {
      return boundary_type::apply_boundary(v, world_size());
   }

   Vector3 cyclic_transpose(const Vector3& p0, const Vector3& p1) const
   {
      return boundary_type::cyclic_transpose(p0, p1, world_size());
   }

   //virtual double cyclic_transpose(double p0, double p1) const override
   //{
   //   return boundary_type::cyclic_transpose(p0, p1, world_size());
   //}

   //virtual position_structid_pair cyclic_transpose(const position_structid_pair& pos_struct_id, const Structure& structure) const override
   //{
   //   const position_structid_pair cyc_pos_struct_id(std::make_pair(cyclic_transpose(pos_struct_id.first, structure.position()), pos_struct_id.second));
   //   return structure.cyclic_transpose(cyc_pos_struct_id, structures_);
   //}

   // THIS SEEMS STRANGE TO PUT THIS HERE. 
   //template<typename T1_>
   //T1_ calculate_pair_CoM(const T1_& p1, const T1_& p2,  const typename element_type_of<T1_>::type& D1, const typename element_type_of<T1_>::type& D2)
   //{
   //   //typedef typename element_type_of< T1_ >::type element_type;   

   //   T1_ retval;

   //   const T1_ p2t(cyclic_transpose(p2, p1));
   //   return modulo(divide(add(multiply(p1, D2), multiply(p2t, D1)), add(D1, D2)), world_size());
   //}




   // --------------------------------------------------------------------------------------------------------------------------------

   bool test_particle_overlap(const Sphere& s, const ParticleID ignore) const
   {
      return test_particle_overlap(s, std::array<ParticleID, 1>({ ignore }));
   }

   bool test_particle_overlap(const Sphere& s, const ParticleID ignore1, const ParticleID ignore2) const
   {
      return test_particle_overlap(s, std::array<ParticleID, 2>({ ignore1, ignore2 }));
   }

   template<typename TSet>
   bool test_particle_overlap(const Sphere& s, const TSet& ignore) const
   {
      ParticleContainerUtils::flag_collector<TSet> fc(ignore);
      boundary_type::overlap_check(pmat_, fc, s);
      return fc.result();
   }

   bool test_particle_overlap(const Sphere& s) const
   {
      ParticleContainerUtils::flag_collector<std::array<ParticleID, 0>> fc;
      boundary_type::overlap_check(pmat_, fc, s);
      return fc.result();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   template<typename TShape>
   particle_id_pair_and_distance_list check_particle_overlap(const TShape& s, const ParticleID ignore) const
   {
      return check_particle_overlap(s, std::array<ParticleID, 1>({ ignore }));
   }

   template<typename TShape>
   particle_id_pair_and_distance_list check_particle_overlap(const TShape& s, const ParticleID ignore1, const ParticleID ignore2) const
   {
      return check_particle_overlap(s, std::array<ParticleID, 2>({ ignore1, ignore2 }));
   }

   template<typename TSet, typename TShape>
   particle_id_pair_and_distance_list check_particle_overlap(const TShape& s, const TSet& ignore) const
   {
      ParticleContainerUtils::neighbor_collector<particle_id_pair_and_distance_list, TSet> nc(ignore);
      boundary_type::overlap_check(pmat_, nc, s);
      return nc.result();
   }

   template<typename TShape>
   particle_id_pair_and_distance_list check_particle_overlap(const TShape& s) const
   {
      ParticleContainerUtils::neighbor_collector<particle_id_pair_and_distance_list, std::array<ParticleID, 0>> nc;
      boundary_type::overlap_check(pmat_, nc, s);
      return nc.result();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_species_type(const SpeciesType& species)
   {
      get_structure_type(species.structure_type_id());      // does StructureType defined in the SpeciesType exists
      species_map_[species.id()] = species;
      particle_pool_[species.id()] = particle_id_set();
   }

   std::size_t num_species() const { return species_map_.size(); }

   const SpeciesType& get_species(const SpeciesTypeID id) const
   {
      species_map::const_iterator i(species_map_.find(id));
      if (species_map_.end() == i) throw not_found(make_string() << "Unknown species (id=" << id << ")");
      return (*i).second;
   }

   SpeciesIteratorRange get_species() const { return SpeciesIteratorRange(species_map_); }



   particle_id_set get_particle_ids(const SpeciesTypeID sid) const
   {
      per_species_particle_id_set::const_iterator i = particle_pool_.find(sid);
      if (i == particle_pool_.end()) throw not_found(make_string() << "Unknown species (id=" << sid << ")");
      return (*i).second;
   }


   template <typename TStructure>
   StructureID add_structure(std::shared_ptr<TStructure> structure)
   {
      get_structure_type(structure->sid());     // does StructureType defined in the Structure

      // add to structure container
      const StructureID sid = structures_.add_structure(structure);
      // update the mapping 'structure_type_id -> set of structure ids'
      structure_pool_[structure->sid()].insert(sid);
      // create a new mapping from structure id -> set of particles
      particleonstruct_pool_[sid] = particle_id_set();
      return sid;
   }

   template <typename TStructIDPair>
   void update_structure(const TStructIDPair& structid_pair)
   {
      if (structures_.has_structure(structid_pair.first))                   // The item was already found
      {
         // FIXME what to do when the structure has particles and moves or changes Structure?!?!
         //  get all particles on the structure
         //  assert that no particles if structure type changes.

         const std::shared_ptr<const Structure> structure(structures_.get_structure(structid_pair.first));
         if (structure->sid() != structid_pair.second->sid())           // If the structuretype changed we need to update the 'structure_type_id->structure ids' mapping
         {
            structure_pool_[structure->sid()].erase(structure->id());
            structure_pool_[structid_pair.second->sid()].insert(structid_pair.first);
         }
         structures_.update_structure(structid_pair);
      }

      // The structure was not yet in the world.
      throw not_found(make_string() << "structure does not exist in the world, add_strucuter first!");
   }

   // Remove structure
   //virtual bool remove_structure(const StructureID id) override
   //{
   //   // TODO
   //   //  -get all particles on structure
   //   //  -only remove if no particles
   //   return base_type::remove_structure(id);
   //}

   StructureContainer::structures_range get_structures() const
   {
      return structures_.get_structures();
   }

   //   std::set<StructureID> get_visible_structures(StructureID sid)
   //   {
   //      return structures_.get_visible_structures(sid);
   //   }

   std::shared_ptr<Structure> get_structure(const StructureID id) const
   {
      return structures_.get_structure(id);
   }


   //// Connectivity related stuff
   //template <typename Tstructure_, typename Tside_enum_>
   //bool connect_structures(Tstructure_ const& structure1, Tside_enum_ const& side1, Tstructure_ const& structure2, Tside_enum_ const& side2)
   //{
   //   return structures_.connect_structures(structure1, side1, structure2, side2);
   //}

   //template <typename Tstructure_, typename Tside_enum_>
   //neighbor_id_vector get_neighbor_info(Tstructure_ const& structure, Tside_enum_ side)
   //{
   //   return structures_.get_neighbor_info(structure, side);
   //}

   //template <typename Tstructure_, typename Tside_enum_>
   //StructureID get_neighbor_id(Tstructure_ const& structure, Tside_enum_ side)
   //{
   //   return structures_.get_neighbor_info(structure, side).first;
   //}

   // Get all the structure ids by Structure id
   structure_id_set get_structure_ids(const StructureTypeID sid) const
   {
      auto i(structure_pool_.find(sid));
      if (i == structure_pool_.end()) throw not_found(make_string() << "Unknown StructureType (id=" << sid << ")");
      return (*i).second;
   }

   //particle_id_set get_particle_ids_on_struct(const StructureID struct_id) const
   //{
   //   auto i(particleonstruct_pool_.find(struct_id));
   //   if (i == particleonstruct_pool_.end()) throw not_found(make_string() << "Unknown structure (id=" << struct_id << ")");
   //   return (*i).second;
   //}

   void add_structure_type(const StructureType& structure)
   {
      structure_type_map_[structure.id()] = structure;
      structure_pool_[structure.id()] = structure_id_set();
   }

   const StructureType& get_structure_type(const StructureTypeID sid) const
   {
      structure_type_map::const_iterator i(structure_type_map_.find(sid));
      if (structure_type_map_.end() == i) throw not_found(make_string() << "Unknown Structure (id=" << sid << ")");
      return (*i).second;
   }

   // TODO
   //virtual structure_types_range get_structure_types() const
   //{
   //   return structure_types_range(structure_type_iterator(structure_type_map_.begin(), structure_types_second_selector_type()),structure_type_iterator(structure_type_map_.end(), structure_types_second_selector_type()),structure_type_map_.size());
   //}




   // --------------------------------------------------------------------------------------------------------------------------------

   structure_id_pair_and_distance_list check_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID sid, double sigma) const
   {
      ParticleContainerUtils::neighbor_collector<structure_id_pair_and_distance_list, std::array<StructureID, 0>> nc;
      surface_overlap_checker(s, old_pos, sid, sigma, nc);
      return nc.result();
   }

   structure_id_pair_and_distance_list check_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID sid, double sigma, const StructureID ignore) const
   {
      ParticleContainerUtils::neighbor_collector<structure_id_pair_and_distance_list, std::array<StructureID, 1>> nc(std::array<StructureID, 1>({ ignore }));
      surface_overlap_checker(s, old_pos, sid, sigma, nc);
      return nc.result();
   }

   structure_id_pair_and_distance_list check_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID sid, double sigma, const StructureID ignore1, const StructureID ignore2) const
   {
      ParticleContainerUtils::neighbor_collector<structure_id_pair_and_distance_list, std::array<StructureID, 2>> nc(std::array<StructureID, 2>({ ignore1, ignore2 }));
      surface_overlap_checker(s, old_pos, sid, sigma, nc);
      return nc.result();
   }

   bool test_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID sid, double sigma) const
   {
      ParticleContainerUtils::flag_collector<std::array<StructureID, 0>> nc;
      surface_overlap_checker(s, old_pos, sid, sigma, nc);
      return nc.result();
   }

   bool test_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID sid, double sigma, const StructureID ignore) const
   {
      ParticleContainerUtils::flag_collector<std::array<StructureID, 1>> nc(std::array<StructureID, 1>({ ignore }));
      surface_overlap_checker(s, old_pos, sid, sigma, nc);
      return nc.result();
   }

   bool test_surface_overlap(const Sphere& s, const Vector3& old_pos, const StructureID sid, double sigma, const StructureID ignore1, const StructureID ignore2) const
   {
      ParticleContainerUtils::flag_collector<std::array<StructureID, 2>> nc(std::array<StructureID, 2>({ ignore1, ignore2 }));
      surface_overlap_checker(s, old_pos, sid, sigma, nc);
      return nc.result();
   }

   template<typename TFun>
   void surface_overlap_checker(const Sphere& s, const Vector3& old_pos, const StructureID sid, double sigma, TFun& nc) const
   {
      structure_map visible_structures;
      for (auto i : structures_.get_visible_structures(sid))
         visible_structures[i] = structures_.get_structure(i);

      for (structure_map::const_iterator i(visible_structures.begin()), e(visible_structures.end()); i != e; ++i)
      {
         Vector3 cyc_old_pos = cyclic_transpose(old_pos, s.position());       // The old position transposed towards the new position (which may also be modified by periodic BC's)
         Vector3 displacement = s.position() - cyc_old_pos;                   // the relative displacement from the 'old' position towards the real new position
         Vector3 cyc_pos = cyclic_transpose(s.position(), (*i).second->position());    // new position transposed to the structure in question
         Vector3 cyc_old_pos2 = cyc_pos - displacement;                                // calculate the old_pos relative to the transposed new position.

         double distance = (*i).second->newBD_distance(cyc_pos, s.radius(), cyc_old_pos2, sigma);
         if (distance < s.radius()) nc(i, distance);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   //virtual std::pair<Vector3, StructureID> cyclic_transpose(const std::pair<Vector3, StructureID>& pos_struct_id, const Structure& structure) const
   //{
   //   const std::pair<Vector3, StructureID> cyc_pos_struct_id(std::make_pair(cyclic_transpose(pos_struct_id.first, structure.position()), pos_struct_id.second));
   //   return structure.cyclic_transpose(cyc_pos_struct_id, structures_);
   //}

   // --------------------------------------------------------------------------------------------------------------------------------

   StructureTypeID get_def_structure_type_id() const
   {
      if (!default_structure_type_id_) throw not_found("Default StructureType is not defined.");
      return default_structure_type_id_;
   }

   void set_def_structure_type_id(const StructureTypeID sid)
   {
      default_structure_type_id_ = sid;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Get the default StructureID of the World 
   StructureID get_def_structure_id() const
   {
      if (!default_structure_id_) throw not_found("Default Structure is not defined.");         // use set_def_structure first!
      return default_structure_id_;
   }

   // Set the default Structure of the World (generates its StructureID)
   void set_def_structure(const std::shared_ptr<CuboidalRegion> cuboidal_region)
   {
      if (!default_structure_id_)
      {
         if (!default_structure_type_id_ || cuboidal_region->sid() != default_structure_type_id_)
            throw illegal_state("Default structure is not of default StructureType");
      }
      else
         throw illegal_state("Default structure is already defined!");

      // add cuboidal_region to structures
      default_structure_id_ = add_structure(cuboidal_region);

      if (cuboidal_region->parent_id() != StructureID(0))
         throw illegal_state("Default structure should have no parent.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Initialize world, copy types from the model
   void initialize(double worldsize, const Model& m)
   {
      // cell_size is scaled according to world_size and cells in X, other dimensions follow.
      pmat_.initialize(worldsize / CompileConfigSimulator::MatrixCellsX);

      for (auto s : m.structures())                            // first copy StructureTypes from Model to World
         add_structure_type(s.second);

      for (auto s : m.get_species())                           // then copy SpeciesTypes from Model to World
         add_species_type(s.second);

      StructureTypeID sid = m.get_def_structure_type_id();     // default StructureTypeID

      Vector3 x = world_size() / 2;                            // create a cuboidal box with size of world
      auto region = std::make_shared<CuboidalRegion>(CuboidalRegion("world", sid, StructureID(0), Box(x, x)));

      set_def_structure_type_id(sid);                          // set StructureTypeID first
      set_def_structure(region);                               // then set region as default structure
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void throwInParticles(const SpeciesTypeID sid, uint numOfParticles, RandomNumberGenerator& rng, bool ignore_structure_overlaps = false, Vector3 bound1 = Vector3(), Vector3 bound2 = Vector3())
   {
      auto species_type = get_species(sid);
      auto stid = species_type.structure_type_id();
      auto structure_type = get_structure_type(stid);
      auto structure_set = get_structure_ids(stid);

      Vector3 ws(world_size());
      if ((bound1 - bound2).length() <= 0.0)         // empty bounds? set to full space
      {
         bound1 = Vector3(0, 0, 0);
         bound2 = ws;
      }

      {
         auto b21 = bound2 - bound1;
         auto bw2 = ws - bound2;
         bool correct = bound1.abs() == bound1 && b21.abs() == b21 && b21.length() > 0 && bw2.abs() == bw2;       // if abs(x)==x : x is positive or zero, quick for testing all elements >= 0
         if (!correct) throw illegal_argument(make_string() << "bounding box is invalid, should be inside world (0 <= b <= ws) and (bound1 < bound2).");
      }

      bool warn = false;
      uint retrycount = GfrdCfg.ThrowInRetryCount;
      while (numOfParticles)
      {
         int idx = structure_set.size() > 1 ? rng.uniform_int(0, static_cast<int>(structure_set.size()) - 1) : 0;

         auto it = structure_set.begin();
         std::advance(it, idx);
         auto structure_id = *it;
         auto structure = structures_.get_structure(structure_id);
         Vector3 pos = structure->random_position(rng);

         std::tie(pos, structure_id) = apply_boundary(position_structid_pair(pos, structure_id));        // does not seem useful, position is generated by the structure, only edges might corner the particle

         bool particle_overlaps = test_particle_overlap(Sphere(pos, species_type.radius()));
         bool surface_overlaps = test_surface_overlap(Sphere(pos, species_type.radius()*GfrdCfg.MINIMAL_SEPARATION_FACTOR), pos, structure_id, species_type.radius());
         bool out_of_bounds = (pos - bound1).abs() != (pos - bound1) || (bound2 - pos).abs() != (bound2 - pos);

         if (!particle_overlaps && (!surface_overlaps || ignore_structure_overlaps) && !out_of_bounds)
         {
            new_particle(Particle(get_species(sid), structure_id, pos));
            numOfParticles--;
            retrycount = GfrdCfg.ThrowInRetryCount;
         }
         else
         {
            retrycount--;
            THROW_UNLESS_MSG(no_space, retrycount > 0, "Add particle failed after " << GfrdCfg.ThrowInRetryCount << " attempts (overlaps or out-of-bounds) " << sid);
            if (!warn && retrycount < 4*GfrdCfg.ThrowInRetryCount/5)
            {
               Log("World").warn() << "Particle placement for " << sid << " may take some time due to restricted available space...";
               warn = true;
            }
         }
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void check() const
   {
      // TODO
      /*
      def check_particle_matrix(self):
         total = sum(len(self.world.get_particle_ids(s.id)) for s in self.world.species)

         if total != self.world.num_particles:
            raise RuntimeError('total number of particles %d != ''self.world.num_particles %d' % (total, self.world.num_particles))

      def check_particles(self):
         for pid_particle_pair in self.world:
            pid = pid_particle_pair[0]
            pos = pid_particle_pair[1].position
            if (pos >= self.world.world_size).any() or (pos < 0.0).any():
               raise RuntimeError('%s at position %s out of the world ''(world size=%g).' % (pid, pos, self.world.world_size))
      */
   };


   // --------------------------------------------------------------------------------------------------------------------------------

   World(const World&) = delete;                                     // don't allow value construct of the world
   World& operator=(const World& w) = delete;                        // don't allow value assignment of the world
   World(World&&) = default;                                         // allow move construction of a world (note that it implementation may be a copy!)
   World& operator=(World&& w) = default;                            // allow move assignment of a world

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   friend class Persistence;

   StructureContainer                  structures_;                     // containing for the structures
   species_map                         species_map_;                    // mapping: SpeciesTypeID -> SpeciesType
   structure_type_map                  structure_type_map_;             // mapping: StructureTypeID -> StructureType
   per_structure_type_structure_id_set structure_pool_;                 // mapping: StructureTypeID -> set of StructureIDs of that StructureType
   per_species_particle_id_set         particle_pool_;                  // mapping: SpeciesTypeID -> set of particle ids of that species
   per_structure_particle_id_set       particleonstruct_pool_;          // mapping: StructureID -> set of particle ids of that structure

   StructureTypeID                     default_structure_type_id_;      // The ID of the default StructureType of the World (does not contain any useful info)
   StructureID                         default_structure_id_;           // The ID of the default Structure of the World (usually this is a CuboidlRegion with size of world)
};

// --------------------------------------------------------------------------------------------------------------------------------
