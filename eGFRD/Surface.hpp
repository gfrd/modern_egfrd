#ifndef SURFACE_HPP
#define SURFACE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "Structure.hpp"
#include "SpeciesType.hpp"
#include <functional>

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TShape>
class Surface : public Structure              // ParticleSimulationStructure<Ttraits_> <- removed class in hierarchy but may be needed to stuff it back in again
{
public:
   typedef TShape                                   shape_type;
   typedef std::pair<double, double>                 DoublePair;
   typedef std::pair<Vector3, DoublePair>            projected_type;
   typedef std::pair<Vector3, bool>                  position_flag_pair_type;


   Surface(const std::string& name, const StructureTypeID sid, const StructureID pid, const shape_type& shape)
      : Structure(name, sid, pid), shape_(shape) {}

   const shape_type& shape() const { return shape_; }

   virtual bool operator==(const Structure& rhs) const override
   {
      Surface const* _rhs(dynamic_cast<Surface const*>(&rhs));
      return _rhs && Structure::operator==(rhs) && shape_ == _rhs->shape_;
   }

   projected_type project_point(const Vector3& pos) const override
   {
      return shape_.project_point(pos);
   }

   //virtual projected_type project_point_on_surface(const Vector3& pos) const override
   //{
   //   return shape_.project_point_on_surface(pos);
   //}

   double distance(const Vector3& pos) const override
   {
      return shape_.distance(pos);
   }

   const Vector3& position() const override
   {
      return shape_.position();
   }

   Vector3 random_position(RandomNumberGenerator& rng) const override
   {
      return shape_.random_position(rng);
   }

   // Vector used to determine whether a particle has crossed the structure
   // Here we return the zero-vector because there is no "sides" to cross
   Vector3 side_comparison_vector() const override
   {
      return Vector3();
   }



   // MS: implement these here for the time being, if left un-implemented) Structure types are abstract, and cannot be instantiated.
   position_structid_pair get_pos_sid_pair(const Structure& target_structure, const Vector3& position, double offset, double rl, RandomNumberGenerator& rng) const override
   {
      UNUSED(rng, rl, offset, position, target_structure);
      // TODO, rethink this whole two level dispatch pattern for template and OOP mix
      return position_structid_pair();
   }

   position_structid_pair_pair get_pos_sid_pair_pair(const Structure& target_structure, const Vector3& position, const SpeciesType& s_orig, const SpeciesType& s_targ, double rl, RandomNumberGenerator& rng) const override
   {
      // TODO, rethink this whole two level dispatch pattern for template and OOP mix

      if (id() == target_structure.id())
      {

         //auto new_positions = geminate_dissociation_positions(rng, s_orig, s_targ, position, rl);
         // geminate_dissociation_positions will produce two new positions close to old_pos taking into account
         // the type of origin_structure and the properties of the two product species
         // (the displacements from old_pos are weighted by the diffusion constants)

         // TODO copied in place, for having reactions with two products work in cubic world
         double r01 = s_orig.radius() + s_targ.radius();
         double D01 = s_orig.D() + s_targ.D();

         double r01l = r01 + rl;
         double r01l_cb = r01l * r01l * r01l;
         double r01_cb = r01 * r01 * r01;
         double diss_vec_length = std::cbrt(rng() * (r01l_cb - r01_cb) + r01_cb);
         auto m = random_vector(diss_vec_length, rng);

         auto pos1 = position - m * s_orig.D() / D01;
         auto pos2 = position + m * s_targ.D() / D01;

         return std::make_pair(std::make_pair(pos1, id()), std::make_pair(pos2, id()));
      }

      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition (CuboidalRegion->CuboidalRegion/CuboidalRegion).");
   }

   position_structid_pair get_pos_sid_pair_2o(Structure const& origin_structure2, const SpeciesTypeID target_sid, const Vector3& CoM, double offset, double reaction_length, RandomNumberGenerator& rng) const override
   {
      UNUSED(rng, reaction_length, offset, CoM, target_sid, origin_structure2);
      // TODO, rethink this whole two level dispatch pattern for template and OOP mix
      return position_structid_pair();
   }






   //virtual position_flag_pair_type deflect(const Vector3& pos0, const Vector3& displacement) const override
   //{
   //   return shape_.deflect(pos0, displacement);
   //}

   //virtual Vector3 deflect_back(const Vector3& pos, const Vector3& u_z) const
   //{
   //   return shape_.deflect_back(pos, u_z);
   //}

   std::size_t hash() const override
   {
      return Structure::hash() ^ std::hash<shape_type>()(shape_);
   }

   std::string as_string() const override
   {
      return make_string() << ", " << shape_;
   }

   const char* type_name() const override { return "Surface"; }

protected:
   friend class Persistence;

   shape_type shape_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* SURFACE_HPP */
