#ifndef SPHERICAL_SURFACE_HPP
#define SPHERICAL_SURFACE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "Surface.hpp"
#include "Sphere.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class StructureContainer;

class SphericalSurface : public Surface < Sphere >
{
public:
   typedef Surface<Sphere> base_type;

   typedef std::pair<Vector3, Vector3>                               position_pair_type;
   typedef std::pair<Vector3, StructureID>                           position_structid_pair;
   typedef std::pair<position_structid_pair, position_structid_pair>   position_structid_pair_pair;

   SphericalSurface(const std::string& name, const StructureTypeID sid, const StructureID pid, const shape_type& shape) noexcept
      : base_type(name, sid, pid, shape) {}

   const char* type_name() const override { return "SphericalSurface"; }

   // *** Simple structure-specific sampling functions 
   Vector3 random_position(RandomNumberGenerator& rng) const override
   {
      UNUSED(rng);
      return Vector3(); // TODO
   }

   Vector3 random_vector(double r, RandomNumberGenerator& rng) const override
   {
      UNUSED(r,rng);
      return Vector3(); // TODO
   }

   Vector3 bd_displacement(double mean, double r, RandomNumberGenerator& rng) const override
   {
      UNUSED(mean,r,rng);
      return Vector3(); // TODO
   }

   //virtual double get_1D_rate_geminate(double k, double r01) const override
   //{
   //   return double(); //TODO
   //}

   double get_1D_rate_surface(double k, double r0) const override
   {
      UNUSED(k,r0);
      return double(); //TODO
   }

   double particle_reaction_volume(double r01, double rl) const override
   {
      UNUSED(r01,rl);
      return double(); //TODO
   }

   double surface_reaction_volume(double r0, double rl) const override
   {
      UNUSED(r0,rl);
      return double(); //TODO
   }

   //virtual Vector3 surface_dissociation_vector(RandomNumberGenerator& rng, double r0, double rl) const override
   //{
   //   return Vector3(); //TODO  
   //}

   //virtual Vector3 surface_dissociation_unit_vector(RandomNumberGenerator& rng) const override
   //{
   //   return Vector3(); //TODO  
   //}

   //// Vector used to determine whether a particle has crossed the structure
   //// Here we return the zero-vector because there is no "sides" to cross
   //virtual Vector3 side_comparison_vector() const override         // MS implemented in Surface (common implementation)
   //{
   //   return Vector3();
   //}

   /* TODO TODO TODO TODO TODO TODO


   virtual position_pair_type geminate_dissociation_positions(RandomNumberGenerator& rng, SpeciesType const& s0, SpeciesType const& s1, Vector3 const& op, double const& rl) const
   {
   return position_pair_type(); //TODO
   }

   virtual position_pair_type special_geminate_dissociation_positions(RandomNumberGenerator& rng, SpeciesType const& s_surf, SpeciesType const& s_bulk, Vector3 const& op_surf, double const& rl) const
   {
   return position_pair_type(); //TODO
   }


   //virtual double minimal_distance(double const& radius) const
   //{
   //    return 0.; // TODO
   //}*/


   double newBD_distance(const Vector3& new_pos, double radius, const Vector3& old_pos, double sigma) const override
   {
      UNUSED(sigma, old_pos, radius);
      return distance(new_pos);
   }

   // *** Boundary condition handling ***
   position_structid_pair apply_boundary(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const override
   {
      UNUSED(structure_container);
      return pos_struct_id;      // The apply_boundary for a spherical surface is trivial (since there is no boundary!)
   }

   position_structid_pair cyclic_transpose(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const override
   {
      UNUSED(structure_container);
      return pos_struct_id;       // Two spherical surface cannot be connected (there is no boundary!)
   }

   /*
   // *** Dynamic dispatch for the structure functions *** //
   // *** 1 *** - One new position
   // This requires a double dynamic dispatch.
   // First dispatch
   virtual position_structid_pair get_pos_sid_pair(Structure const& target_structure, Vector3 const& position,
   double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   return target_structure.get_pos_sid_pair_helper(*this, position, offset, reaction_length, rng);
   }
   // Second dispatch
   virtual position_structid_pair get_pos_sid_pair_helper(CuboidalRegion const& origin_structure, Vector3 const& position,
   double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<CuboidalRegion >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(SphericalSurface const& origin_structure, Vector3 const& position,
   double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<SphericalSurface >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(CylindricalSurface const& origin_structure, Vector3 const& position,
   double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<CylindricalSurface >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(DiskSurface const& origin_structure, Vector3 const& position,
   double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<DiskSurface >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(PlanarSurface const& origin_structure, Vector3 const& position,
   double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<PlanarSurface >(origin_structure, position, offset, rl, rng);
   }
   // The template function that defines the actual final dispatch procedure.
   template<typename Tstruct_>
   position_structid_pair get_pos_sid_pair_helper_any(Tstruct_ const& origin_structure, Vector3 const& position,
   double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   // redirect to structure function with well-defined typing
   return ::get_pos_sid_pair(origin_structure, *this, position, offset, rl, rng);
   };

   // *** 2 *** - Two new positions
   // Same principle as above, but different return type
   // First dispatch
   virtual position_structid_pair_pair get_pos_sid_pair_pair(Structure const& target_structure, Vector3 const& position,
   SpeciesType const& s1, SpeciesType const& s2, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   return target_structure.get_pos_sid_pair_pair_helper(*this, position, s1, s2, reaction_length, rng);
   }
   // Second dispatch
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(CuboidalRegion const& origin_structure, Vector3 const& position,
   SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<CuboidalRegion >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(SphericalSurface const& origin_structure, Vector3 const& position,
   SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<SphericalSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(CylindricalSurface const& origin_structure, Vector3 const& position,
   SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<CylindricalSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(DiskSurface const& origin_structure, Vector3 const& position,
   SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<DiskSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(PlanarSurface const& origin_structure, Vector3 const& position,
   SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<PlanarSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   // The template function that defines the actual final dispatch procedure.
   template<typename Tstruct_>
   position_structid_pair_pair get_pos_sid_pair_pair_helper_any(Tstruct_ const& origin_structure, Vector3 const& position,
   SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   // redirect to structure function with well-defined typing
   return ::get_pos_sid_pair_pair(origin_structure, *this, position, s_orig, s_targ, rl, rng);
   };

   // *** 3 *** - Pair reactions => two origin structures
   // First dispatch
   //     // Overloading get_pos_sid_pair with signature (origin_structure2, target_Structure_id, ...)
   //    virtual position_structid_pair get_pos_sid_pair(Structure const& origin_structure2, SpeciesTypeID const& target_sid, Vector3 const& CoM,
   //                                                         double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   //     {
   //         // this just redirects
   //         return this->get_pos_sid_pair_2o(origin_structure2, target_sid, CoM, offset, reaction_length, rng);
   //     }
   //     // The actual implementation of the first dispatch
   virtual position_structid_pair get_pos_sid_pair_2o(Structure const& origin_structure2, SpeciesTypeID const& target_sid,
   Vector3 const& CoM, double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   return origin_structure2.get_pos_sid_pair_2o_helper(*this, target_sid, CoM, offset, reaction_length, rng);
   }
   // Second dispatch
   virtual position_structid_pair get_pos_sid_pair_2o_helper(CuboidalRegion const& origin_structure1, SpeciesTypeID const& target_sid,
   Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<CuboidalRegion >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(SphericalSurface const& origin_structure1, SpeciesTypeID const& target_sid,
   Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<SphericalSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(CylindricalSurface const& origin_structure1, SpeciesTypeID const& target_sid,
   Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<CylindricalSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(DiskSurface const& origin_structure1, SpeciesTypeID const& target_sid,
   Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<DiskSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(PlanarSurface const& origin_structure1, SpeciesTypeID const& target_sid,
   Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<PlanarSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   // The template function that defines the actual final dispatch procedure.
   template<typename Tstruct_>
   position_structid_pair get_pos_sid_pair_2o_helper_any(Tstruct_ const& origin_structure1, SpeciesTypeID const& target_sid, Vector3 const& CoM,
   double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   // This method has to figure out where the product will be placed in case of a bimolecular reaction.
   // As a default, we place particles on the substructure or the lower-dimensional structure. If the structures
   // have the same structure type (=> same dimensionality) it does not matter on which structure we put the product,
   // as long as it has the structure type id of the product species. This is handled in cases '1' below.

   // 1 - Check whether one of the structures is the parent of the other. If yes, the daughter structure is the target.
   if (this->is_parent_of_or_has_same_sid_as(origin_structure1) && origin_structure1.has_valid_target_sid(target_sid))
   // origin_structure1 is target
   return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);

   else if (origin_structure1.is_parent_of_or_has_same_sid_as(*this) && this->has_valid_target_sid(target_sid))
   // this structure is target
   return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);

   // 2 - Check which structures has the lower dimensionality / particle degrees of freedom, and put the product there.
   else if (origin_structure1.shape().dof() < this->shape().dof() && origin_structure1.has_valid_target_sid(target_sid))
   // origin_structure1 is target
   return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);

   else if (this->shape().dof() < origin_structure1.shape().dof() && this->has_valid_target_sid(target_sid))
   // this structure is target
   return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);

   else throw propagation_error("Invalid target structure type: does not match product species structure type or has wrong hierarchy or dimensionality.");
   }

   //     // *** 4 *** - Generalized functions for pair reactions with two origin structures and one target structure
   //     // NOTE: This is yet unused, but possibly useful in the future.
   //     // Overloading get_pos_sid_pair again with signature (origin_structure2, target_structure, ...) and introducing
   //     // a triple dynamic dispatch.
   //     virtual position_structid_pair get_pos_sid_pair(Structure const& origin_structure2, Structure const& target_structure, Vector3 const& position,
   //                                                          double const& offset, double const& reaction_length, RandomNumberGenerator const& rng) const
   //     {
   //         return origin_structure2.get_pos_sid_pair_helper1(*this, target_structure, position, offset, reaction_length, rng);
   //     }
   */

   ///*** Formerly used functions of the Morelli scheme ***/
   //// DEPRECATED
   //virtual double drawR_gbd(double const& rnd, double const& r01, double const& dt, double const& D01, double const& v) const
   //{    
   //    return double(); // TODO
   //}
   //// DEPRECATED
   //virtual double p_acceptance(double const& k_a, double const& dt, double const& r01, Vector3 const& ipv, 
   //                            double const& D0, double const& D1, double const& v0, double const& v1) const
   //{    
   //    return double(); //TODO
   //}
   //// DEPRECATED
   //virtual Vector3 dissociation_vector( RandomNumberGenerator& rng, double const& r01, double const& dt, 
   //                                            double const& D01, double const& v ) const
   //{
   //    return Vector3(); //TODO
   //}

   //virtual void accept(ImmutativeStructureVisitor const& visitor) const
   //{
   //    visitor(*this);
   //}

   //virtual void accept(MutativeStructureVisitor const& visitor)
   //{
   //    visitor(*this);
   //}
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* SPHERICAL_SURFACE_HPP */
