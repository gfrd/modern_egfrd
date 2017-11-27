#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "Surface.hpp"
#include "Disk.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class StructureContainer;

class DiskSurface : public Surface < Disk >
{
public:
   typedef Surface<Disk> base_type;

   typedef std::pair<Vector3, Vector3>                               position_pair_type;
   typedef std::pair<Vector3, StructureID>                           position_structid_pair;
   typedef std::pair<position_structid_pair, position_structid_pair>   position_structid_pair_pair;

   DiskSurface(const std::string& name, const StructureTypeID sid, const StructureID pid, const shape_type& shape) noexcept
      : base_type(name, sid, pid, shape) {}


   const char* type_name() const override { return "DiskSurface"; }

   // Produce a "random vector" on the disk; returns the same as random_position()
   Vector3 random_vector(double r, RandomNumberGenerator& rng) const override
   {
      UNUSED(r, rng);
      return shape_.position();
   }

   // BD displacement = zero vector because there is only one legal position on the disk
   Vector3 bd_displacement(double mean, double r, RandomNumberGenerator& rng) const override
   {
      UNUSED(mean,r,rng);
      return Vector3();  // multiply(shape_.unit_z(), 0.0);  // TODO is there not cheaper way to pass a zero vector?
   }

   //// Rate for binding to particle on the structure
   //virtual double get_1D_rate_geminate(double k, double r01) const override
   //{
   //   // Same as for particle-particle reactions on the cylinder
   //   return k;
   //}

   // Rate for binding to the structure
   double get_1D_rate_surface(double k, double r0) const override
   {
      UNUSED(r0);
      return k;
   }

   // Reaction volume for binding to particle in the structure
   double particle_reaction_volume(double r01, double rl) const override
   {
      UNUSED(r01);
      // The disk can only hold 1 particle; this function therefore never should be called.        
      return rl;

   //   // FIXME: This reaction volume is only correct for particles coming from a rod.
   //   // If the interaction partner of the disk particle comes from a plane or from
   //   // the bulk, a different volume factor shall be used. Then the return value of
   //   // this function would depend on properties of the asker -> how to do???
   }

   // Reaction volume for binding to the structure
   double surface_reaction_volume(double r0, double rl) const override
   {
      UNUSED(r0);
      // The reaction volume for a particle on the rod interacting with a disk;
      // should be the same as for two particles interacting with each other on the rod.
      return rl;
   }

   //// Vector of dissociation from the structure into the bulk
   //virtual Vector3 surface_dissociation_vector(RandomNumberGenerator& rng, double offset, double rl) const override
   //{
   //   // This function produces a position for a particle unbinding from the disk.
   //   // It should lie within the reaction volume around the disk.
   //   // Note that this is the same code as for the CylindricalSurface.
   //   double X(rng.uniform(0., 1.));
   //   double   const disk_radius = shape_.radius();
   //   Vector3 const unit_z = shape_.unit_z();

   //   // Calculate the length of the vector first
   //   double const rrl(disk_radius + offset + rl);
   //   double const rrl_sq(gsl_pow_2(rrl));
   //   double const rr_sq(gsl_pow_2(disk_radius + offset));
   //   // Create a random length between rr_sq and rrl_sq
   //   double const diss_vec_length(sqrt(rr_sq + X * (rrl_sq - rr_sq)));

   //   // Create a 3D vector with totally random orientation
   //   Vector3 v(rng.uniform(0., 1.) - .5, rng.uniform(0., 1.) - .5, rng.uniform(0., 1.) - .5);
   //   // Subtract the part parallel to the axis to get the orthogonal components and normalize
   //   // This creates a normed random vector in the disk plane
   //   v = (v - unit_z * Vector3::dot(unit_z, v)).normal();

   //   // Return the created vector with the right length
   //   return v * MINIMAL_SEPARATION_FACTOR * diss_vec_length;
   //   // TODO define a global MINIMAL_SEPARATION_FACTOR also for BD mode
   //}

   //// Normed direction of dissociation from the structure to parent structure
   //virtual Vector3 surface_dissociation_unit_vector(RandomNumberGenerator& rng) const override
   //{
   //   return shape_.unit_z(); // FIXME
   //}

   //// Vector used to determine whether a particle has crossed the structure
   //// Here we return the zero-vector because there is no "sides" to cross
   //virtual Vector3 side_comparison_vector() const override         // MS implemented in Surface (common implementation)
   //{
   //   return Vector3();
   //}



   /* TODO TODO TODO TODO TODO TODO



   // Positions created at dissociation of one particle on the structure into two particles on the structure
   virtual position_pair_type geminate_dissociation_positions( RandomNumberGenerator& rng, SpeciesType const& s0, SpeciesType const& s1, Vector3 const& op,
   double const& rl ) const
   {
   // The positions of a particle dissociating into two new ones on the disk; should never happen,
   // therefore this function just returns a dummy positions pair (2x the disk center)
   return position_pair_type( shape_.position(), shape_.position() );
   }

   // Positions created at dissociation of one particle on the structure into two particles, one of which ends up in the bulk
   virtual position_pair_type special_geminate_dissociation_positions( RandomNumberGenerator& rng, SpeciesType const& s_disk, SpeciesType const& s_diss,
   Vector3 const& reactant_pos, double const& rl ) const
   {
   // This function produces two new positions for a dissociating particle in the case
   // that one stays on the surface and the other one changes to the parent structure.
   // TODO We have to distinguish between sink and cap here and between dissociation onto
   // the rod and into the bulk/plane!

   // Note: s_disk = disk-bound species, s_diss = dissociating species (may go to bulk or plane)

   double const disk_radius( shape_.radius() );
   double const r01( s_disk.radius() + s_diss.radius() );

   // The following is the additional distance that we have to pass to surface_dissociation_vector() below
   // to place the unbinding particle in contact with the disk particle or the disk, whatever has the
   // larger radius. surface_dissociation_vector() will add it to the disk_radius.
   double offset( s_disk.radius() > disk_radius ? r01 - disk_radius : s_diss.radius());

   position_pair_type pp01;
   // Particle 0 is the one that stays on the origin structure. It does not move.
   pp01.first  = reactant_pos;
   // Particle 1 is the one that unbinds from the disk and is placed in the reaction volume around it.
   pp01.second = add(reactant_pos, surface_dissociation_vector(rng, offset, rl));

   return pp01;
   }*/

   // Used by newBDPropagator
   double newBD_distance(const Vector3& new_pos, double radius, const Vector3& old_pos, double sigma) const override
   {
      UNUSED(radius);

      double disk_radius = shape_.radius();
      auto new_pos_rz = shape().to_internal(new_pos);
      auto old_pos_rz = shape().to_internal(old_pos);
      if (new_pos_rz.Y() * old_pos_rz.Y() < 0 && ((new_pos_rz.X() < disk_radius) || (old_pos_rz.X() < disk_radius)))
         return -1.0 * distance(new_pos) + sigma;
      return distance(new_pos) + sigma;
   }

   //virtual double minimal_distance(double const& radius) const
   //{
   //    // TODO
   //    double cylinder_radius = shape_.radius();
   //    // Return minimal distance *to* surface.
   //    return (cylinder_radius + radius) * traits_type::MINIMAL_SEPARATION_FACTOR - cylinder_radius;
   //}

   // *** Boundary condition handling ***
   // FIXME This is a mess but it works. See ParticleContainerBase.hpp for explanation.
   position_structid_pair apply_boundary(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const override
   {
      UNUSED(structure_container);
      return pos_struct_id;   // This seems a little strange, but we assume that particles are immobile on the disk
   }

   position_structid_pair cyclic_transpose(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const override
   {
      UNUSED(structure_container);
      return pos_struct_id;   // Disks can also not be connected, so no cyclic transpose.
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
   if( this->is_parent_of_or_has_same_sid_as(origin_structure1) && origin_structure1.has_valid_target_sid(target_sid) )
   // origin_structure1 is target
   return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);

   else if( origin_structure1.is_parent_of_or_has_same_sid_as(*this) && this->has_valid_target_sid(target_sid) )
   // this structure is target
   return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);

   // 2 - Check which structures has the lower dimensionality / particle degrees of freedom, and put the product there.
   else if( origin_structure1.shape().dof() < this->shape().dof() && origin_structure1.has_valid_target_sid(target_sid) )
   // origin_structure1 is target
   return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);

   else if( this->shape().dof() < origin_structure1.shape().dof() && this->has_valid_target_sid(target_sid) )
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

   /*** Formerly used functions of the Morelli scheme ***/
   //// DEPRECATED
   //virtual double drawR_gbd(double const& rnd, double const& r01, double const& dt, double const& D01, double const& v) const
   //{
   //     // TODO: This is part of the old BD scheme and should be removed at some point
   //    return drawR_gbd_1D(rnd, r01, dt, D01, v);
   //}
   //// DEPRECATED
   //virtual double p_acceptance(double const& k_a, double const& dt, double const& r01, Vector3 const& ipv, 
   //                            double const& D0, double const& D1, double const& v0, double const& v1) const
   //{
   //    // TODO: This is part of the old BD scheme and should be removed at some point
   //    /*
   //        The I_bd factors used for calculating the acceptance probability are dependent on the direction 
   //        of the overlap step (r = r_1 - r_0), compared to the direction of the drift. 
   //        The I_bd factors are defined for a particle creating an overlap comming from the right (r < 0).
   //        Since the I_bd terms calulated here are for the backward move, we have to invert their drifts.
   //        When the particle comes from the left (r > 0) we have to invert its drift again.

   //        ---Code below is used for drift dependent backstep.

   //        double numerator = g_bd_1D(ipv, r01, dt, D0, -v0);
   //        double denominator = g_bd_1D(ipv, r01, dt, D0, v0)*exp( ipv/abs_ipv*(abs_ipv - r01)*v/D01 );
   //        double correction = numerator/denominator;
   //        if( ipv < 0 )
   //            return correction*( k_a * dt / ( I_bd_1D(r01, dt, D0, -v0) + I_bd_1D(r01, dt, D1, v1) ) );
   //        else
   //            return correction*( k_a * dt / ( I_bd_1D(r01, dt, D0, v0) + I_bd_1D(r01, dt, D1, -v1) ) );

   //        Also change v -> -v in drawR for the dissociation move.
   //    */

   //    return 0.5*( k_a * dt / ( I_bd_1D(r01, dt, D0, v0) + I_bd_1D(r01, dt, D1, v1) ) );

   //}
   //// DEPRECATED
   //virtual Vector3 dissociation_vector( RandomNumberGenerator& rng, double const& r01, double const& dt, 
   //                                            double const& D01, double const& v ) const
   //{
   //    // TODO: This is part of the old BD scheme and should be removed at some point
   //    return random_vector(drawR_gbd(rng.uniform(0., 1.), r01, dt, D01, v), rng);
   //}
   //
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
