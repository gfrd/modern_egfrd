#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "MatrixSpace.hpp"
#include "exceptions.hpp"
#include "ParticleContainer.hpp"
#include "StructureContainer.hpp"
#include "SerialIDGenerator.hpp"

#include "Transaction.hpp"
#include "Shell.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ParticleContainerUtils
{
   // The flag_collector is class for the overlap_check that flags when an overlap occurs (not who or where)
   template<typename TIgnoreSet>
   struct flag_collector
   {
      flag_collector(const TIgnoreSet& ignore = TIgnoreSet()) : ignore_(ignore), result_() { }

      template<typename TIter>
      void operator()(TIter& i, double)
      {
         if (result_) return; // already set, quick exit
         bool found = std::find(ignore_.cbegin(), ignore_.cend(), (*i).first) != ignore_.cend();
         if (!found) result_ = true;
      }

      bool result() const { return result_; }

   private:
      const TIgnoreSet&      ignore_;
      bool                   result_;
   };



   // The neighbor_collector is list of pairs (object, distance) which is sorted on the second part (the distance) when requesting results.
   template<typename TObjectIDPairDistanceList, typename TIgnoreSet>
   struct neighbor_collector
   {
      typedef TObjectIDPairDistanceList  list_type;

      neighbor_collector(const TIgnoreSet& ignore = TIgnoreSet()) : ignore_(ignore), result_() { }

      template<typename TIter>
      void operator()(TIter& i, double distance)
      {
         bool found = std::find(ignore_.cbegin(), ignore_.cend(), (*i).first) != ignore_.cend();
         if (!found) result_.emplace_back(std::make_pair(*i, distance));
      }

      const list_type&& result(bool nosort = false)
      {
         if (!nosort) std::sort(result_.begin(), result_.end(),             // sort with lambda on distance!
            [](const typename list_type::value_type &lhs, const typename list_type::value_type &rhs) { return lhs.second < rhs.second; });
         return std::move(result_);
      }

   private:
      const TIgnoreSet&                    ignore_;
      list_type                            result_;
   };


   // The overlap_check is a MatrixSpace each_neighbor_xx collector, that selects those 'neighbors' that overlap with the given shape (usally sphere or cylinder).
   // results are stored in a xx_collector object, with distances defined as 'neighbor-edge' to centre of shape (not the overlap or radius)
   template<typename TMatrixSpace, typename TNeighborCollector, typename TShape>
   class overlap_check
   {
      typedef typename TMatrixSpace::const_iterator         const_iterator;

   public:
      overlap_check(TNeighborCollector& c, const TShape& shape) : c_(c), shape_(shape) { }

      void operator()(const_iterator i, const Vector3& offset) const
      {
         typename const_iterator::reference item(*i);

         //const double dist(item.second.shape().offset(offset).distance(shape_.position()));
         //if (dist < shape_.radius())
         //{
         //   c_(i, dist);     // put item in the collector object
         //}

         // item.second is always a sphere, calculate distance from (edge of shape_) to center of that sphere, minus its radius gives:
         // distance < 0 : penetrating, distance = 0 : touching, distance > 0 : gap
         double distance = shape_.distance(item.second.shape().offset(offset).position()) - item.second.shape().radius();
         if (distance < 0.0)
         {
            c_(i, distance);     // put item in the collector object
         }
      }

   private:
      TNeighborCollector&        c_;        // structure (overlap checker) storing the overlapping particles.
      const TShape&              shape_;    // The spherical particle whose neighbors are being checked.
   };
};

// --------------------------------------------------------------------------------------------------------------------------------

struct WorldNoBounds
{
   // This is the normal world (without periodic/cyclic boundary conditions)

   static Vector3 apply_boundary(const Vector3& v, const Vector3& world_size)
   {
      UNUSED(world_size);
      return v;
   }

   static Vector3 cyclic_transpose(const Vector3& p0, const Vector3& p1, const Vector3& world_size)
   {
      UNUSED(p1, world_size);
      return p0;
   }

   static double distance(const Vector3& p0, const Vector3& p1, const Vector3& world_size)
   {
      UNUSED(world_size);
      return (p0 - p1).length();
   }

   template<typename TMatrixSpace, typename TNeighborCollector, typename TShape>
   static void overlap_check(TMatrixSpace& ms, TNeighborCollector& col, const TShape& shape)
   {
      auto oc = ParticleContainerUtils::overlap_check<TMatrixSpace, TNeighborCollector, TShape>(col, shape);
      ms.each_neighbor(ms.index(shape.position()), oc);
   }

   template<typename TMatrixSpace, typename TNeighborCollector>
   static void each_neighbor(const TMatrixSpace& ms, TNeighborCollector& col, const Vector3& pos)
   {
      ms.each_neighbor(ms.index(pos), col);
   }
};

// --------------------------------------------------------------------------------------------------------------------------------

struct WorldCyclic
{
   // This is the default world (with periodic/cyclic boundary conditions)

   static Vector3 apply_boundary(const Vector3& v, const Vector3& world_size)
   {
      return v.modulo(world_size);
   }

   static Vector3 cyclic_transpose(Vector3  p0, Vector3  p1, const Vector3& world_size)
      // selects the copy of p0 (over the periodic boundaries) such that the distance between p0 and p1 is minimized over over the periodic boundary conditions.
   {
      return cyclic::cyclic_transpose(p0, p1, world_size);
   }

   static double distance(const Vector3& p0, const Vector3& p1, const Vector3& world_size)
   {
      return cyclic::distance_cyclic(p0, p1, world_size);
   }

   template<typename TMatrixSpace, typename TNeighborCollector, typename TShape>
   static void overlap_check(TMatrixSpace& ms, TNeighborCollector& col, const TShape& shape)
   {
      auto oc = ParticleContainerUtils::overlap_check<TMatrixSpace, TNeighborCollector, TShape>(col, shape);
      ms.each_neighbor_cyclic(ms.index(shape.position()), oc);
   }

   template<typename TMatrixSpace, typename TNeighborCollector>
   static void each_neighbor(const TMatrixSpace& ms, TNeighborCollector& col, const Vector3& pos)
   {
      ms.each_neighbor_cyclic(ms.index(pos), col);
   }

};

// --------------------------------------------------------------------------------------------------------------------------------

class ParticleContainerImpl : public ParticleContainer
{
public:
   using boundary_type = CompileConfigSimulator::TBoundCondition;
   using particle_matrix_type = MatrixSpace<Particle, ParticleID, CompileConfigSimulator::MatrixCellsX, CompileConfigSimulator::MatrixCellsY, CompileConfigSimulator::MatrixCellsZ>;
   using particle_id_pair_and_distance = std::pair<particle_id_pair, double>;

   ParticleContainerImpl() : pmat_() { }

   void initialize(double cell_size) { pmat_.initialize(cell_size); }

   std::size_t num_particles() const override { return pmat_.size(); }

   Vector3 world_size() const override { return pmat_.world_size(); }

   double cell_size() const override { return pmat_.cell_size(); }

   std::array<uint, 3> matrix_size() const override { return pmat_.matrix_size(); }

   // --------------------------------------------------------------------------------------------------------------------------------

   particle_id_pair get_particle(const ParticleID id, bool& found) const
   {
      particle_matrix_type::const_iterator i(pmat_.find(id));
      if (pmat_.end() == i)
      {
         found = false;
         return particle_id_pair();
      }

      found = true;
      return *i;
   }

   virtual particle_id_pair get_particle(const ParticleID pid) const override
   {
      particle_matrix_type::const_iterator i(pmat_.find(pid));
      if (pmat_.end() == i) throw not_found(make_string() << "No such particle: id=" << pid);
      return *i;
   }

   virtual bool has_particle(const ParticleID pid) const override
   {
      return pmat_.end() != pmat_.find(pid);
   }

   virtual Transaction* create_transaction() override
   {
      return new TransactionImpl<ParticleContainerImpl>(*this);
   }

   virtual particle_id_pair_generator get_particles() const override
   {
      return agi::iteratorRange<particle_id_pair>(pmat_);
   }

   virtual particle_id_pair new_particle(const Particle&& p) override
   {
      particle_id_pair pip(std::make_pair(pidgen_(), std::move(p)));
      VERIFY(update_particle(pip));    // check really new
      return pip;
   }
   virtual bool update_particle(const particle_id_pair& pip) override
   {
      return pmat_.update(pip).second;
   }

   virtual bool remove_particle(const ParticleID pid) override
   {
      return pmat_.erase(pid);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   ///// Structure methods
   ///// The following are mostly wrappers of the methods defined in StructureContainer
   ///// and handle structure connectivity and boundary application.

   //// Check for whether a container contains a structure with a certain ID
   //virtual bool has_structure(const StructureID& id) const override
   //{
   //   return structures_.has_structure(id);
   //}

   //// Getter structure_id -> structure
   //virtual std::shared_ptr<Structure> get_structure(const StructureID& id) const override
   //{
   //   return structures_.get_structure(id);
   //}

   //// Get range of structures in container
   ////virtual structures_range get_structures() const
   ////{
   ////   return structures_.get_structures_range();
   ////}

   //// Getter Structure_id -> some structure of that type
   //virtual std::shared_ptr<Structure> get_some_structure_of_type(const SpeciesTypeID& sid) const override
   //{
   //   return structures_.get_some_structure_of_type(sid);
   //}

   //// Update structure wrapper
   //template <typename Tstructid_pair_>
   //bool update_structure(Tstructid_pair_ const& structid_pair)
   //{
   //   return structures_.update_structure(structid_pair);
   //}

   //// Remove structure wrapper
   //virtual bool remove_structure(const StructureID& id) override
   //{
   //   return structures_.remove_structure(id);
   //}

   // Get all structures close to a position pos, taking care of structure types
   // and an ignore parameter (defining a set of ignored structures)
   // The actual distance measurement is performed by method structure->distance(cyc_pos) below,
   // which is implemented in the respective structure (sub-) classes.
   //virtual structure_id_pair_and_distance_list* get_close_structures(const Vector3& pos, const StructureID& current_struct_id, const StructureID& ignore) const override
   //{
   //   return nullptr;
   //   //const structure_id_set visible_structure_IDs(structures_.get_visible_structures(current_struct_id));

   //   //// Get and temporarily store all the visible structures (upto now we only had their IDs)
   //   //structure_map visible_structures;
   //   //for (auto i(visible_structure_IDs.begin()), e(visible_structure_IDs.end()); i != e; ++i)
   //   //{
   //   //   visible_structures[(*i)] = get_structure(*i);
   //   //}

   //   //// Calculate the distances and store the surface,distance tuple in the overlap checker (in a list is sorted by distance).
   //   //ParticleContainerUtils::neighbor_collector<structure_id_pair_and_distance_list, std::array<StructureID, 1>> checker(std::array<StructureID, 1>({ ignore }));
   //   //for (structure_map::const_iterator i(visible_structures.begin()), e(visible_structures.end()); i != e; ++i)
   //   //{
   //   //   const Vector3 cyc_pos(cyclic_transpose(pos, ((*i).second)->position()));
   //   //   // Here we perform the actual distance measurement
   //   //   const double dist((*i).second->distance(cyc_pos));
   //   //   checker(i, dist);
   //   //}
   //   //return checker.result();
   //}

protected:
   friend class Persistence;

   SerialIDGenerator<ParticleID>       pidgen_;                         // SerialIDGenerator used to produce the unique ids for the particles
   particle_matrix_type                pmat_;          // the structure (MatrixSpace) containing the particles.
};

// --------------------------------------------------------------------------------------------------------------------------------
