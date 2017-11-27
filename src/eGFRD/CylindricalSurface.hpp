#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "Surface.hpp"
#include "Cylinder.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class StructureContainer;

class CylindricalSurface : public Surface < Cylinder >
{
public:
   typedef Surface<Cylinder> base_type;

   typedef std::pair<Vector3, Vector3>                                           position_pair_type;
   typedef std::pair<Vector3, StructureID>                                       position_structid_pair;
   typedef std::pair<position_structid_pair, position_structid_pair>   position_structid_pair_pair;

   CylindricalSurface(const std::string& name, const StructureTypeID sid, const StructureID pid, const shape_type& shape) noexcept
      : base_type(name, sid, pid, shape) {}

   const char* type_name() const override { return "CylindricalSurface"; }

   // Produce a random vector along the axis of the cylinder
   virtual Vector3 random_vector(double r, RandomNumberGenerator& rng) const override
   {
      return shape_.unit_z() * (rng.uniform_int(0, 1) * 2 - 1) * r;
   }

   // BD displacement for BD on the axis of the cylinder
   virtual Vector3 bd_displacement(double mean, double r, RandomNumberGenerator& rng) const override
   {
      return shape_.unit_z() * rng.normal(mean, r);
   }

   //// Rate for binding to particle on the structure
   //virtual double get_1D_rate_geminate(double k, double r01) const override
   //{
   //   return k;
   //}

   // Rate for binding to the structure
   virtual double get_1D_rate_surface(double k, double r0) const override
   {
      return k / (2 * M_PI * (shape_.radius() + r0));
   }

   // Reaction volume for binding to particle in the structure
   double particle_reaction_volume(double r01, double rl) const override
   {
      UNUSED(r01);
      return rl;
   }

   // Reaction volume for binding to the structure
   double surface_reaction_volume(double r0, double rl) const override
   {
      double rc(shape_.radius() + r0);
      double rcl(rc + rl);
      double rcl_sq(rcl * rcl);
      double rc_sq(rc * rc);
      return M_PI * (rcl_sq - rc_sq);
   }

   //// Vector of dissociation from the structure into the bulk
   //virtual Vector3 surface_dissociation_vector(RandomNumberGenerator& rng, double r0, double rl) const override
   //{
   //   double X(rng.uniform(0., 1.));
   //   double const rod_radius = shape_.radius();
   //   Vector3 const unit_z = shape_.unit_z();

   //   // Calculate the length of the vector first
   //   double const rrl(rod_radius + r0 + rl);
   //   double const rrl_sq(gsl_pow_2(rrl));
   //   double const rr_sq(gsl_pow_2(rod_radius + r0));
   //   // Create a random length between rr_sq and rrl_sq
   //   double const diss_vec_length(sqrt(rr_sq + X * (rrl_sq - rr_sq)));

   //   // Create a 3D vector with totally random orientation
   //   Vector3 v(rng.uniform(0., 1.) - .5, rng.uniform(0., 1.) - .5, rng.uniform(0., 1.) - .5);
   //   // Subtract the part parallel to the axis to get the orthogonal components and normalize
   //   // This creates a normed random vector orthogonal to the cylinder axis
   //   v = (v - unit_z * Vector3::dot(unit_z, v)).normal();

   //   // Return the created vector with the right length
   //   return v * diss_vec_length;
   //}

   //// Normed direction of dissociation from the structure to parent structure
   //virtual Vector3 surface_dissociation_unit_vector(RandomNumberGenerator& rng) const override
   //{
   //   Vector3 const unit_z = shape_.unit_z();
   //   // Create a 3D vector with totally random orientation
   //   Vector3 v(rng.uniform(0., 1.) - .5, rng.uniform(0., 1.) - .5, rng.uniform(0., 1.) - .5);
   //   // Subtract the part parallel to the axis to get the orthogonal components and normalize
   //   // This creates a normed random vector orthogonal to the cylinder axis
   //   return (v - unit_z * Vector3::dot(unit_z, v)).normal();
   //}

   //// Vector used to determine whether a particle has crossed the structure
   //// Here we return the zero-vector because there is no "sides" to cross
   //virtual Vector3 side_comparison_vector() const override         // MS implemented in Surface (common implementation)
   //{
   //   return Vector3();
   //}


   /* TODO TODO TODO TODO TODO TODO


   // Positions created at dissociation of one particle on the structure into two particles on the structure
   virtual position_pair_type geminate_dissociation_positions( RandomNumberGenerator& rng, SpeciesType const& s0, SpeciesType const& s1, Vector3 const& op, double const& rl ) const
   {
   double const r01( s0.radius() + s1.radius() );
   double const D01( s0.D() + s1.D() );

   double const X( rng.uniform(0.,1.) );

   double const diss_vec_length( X*rl + r01 );

   Vector3 const m( random_vector( diss_vec_length, rng ) );

   return position_pair_type( op - m * s0.D() / D01,
   op + m * s1.D() / D01 );
   }

   // Positions created at dissociation of one particle on the structure into two particles, one of which ends up in the bulk
   virtual position_pair_type special_geminate_dissociation_positions( RandomNumberGenerator& rng, SpeciesType const& s_surf, SpeciesType const& s_bulk,  Vector3 const& op_surf, double const& rl ) const
   {
   // FIXME FIXME FIXME What the hell is all of that calculation? What are these angles / transformations?
   // This should just place the dissociating particle in the reaction volume close to the cylinder!
   // NOTE Probably just copied from PlanarSurface ?

   double const rod_radius( shape_.radius() );

   //Species living on the rod should have a larger radius than the rod.
   assert( rod_radius < s_surf.radius() ); // FIXME This is dangerous to put here! Revise!

   double const r01( s_bulk.radius() + s_surf.radius() );
   double const D01( s_bulk.D() + s_surf.D() );
   double const D_bulk_D01( s_bulk.D() / D01 );
   double const D_surf_D01( s_surf.D() / D01 );

   //Commented code for direct binding case with c.o.m. reaction.
   //double const theta_min( asin(rod_radius / r01) );
   double const theta_min( asin( (rod_radius + s_bulk.radius()) / r01) );
   double const theta( theta_min + rng.uniform(0.,1.) * (M_PI - 2 * theta_min) );
   double const phi( rng.uniform(0.,1.) * 2 * M_PI );

   double const X( rng.uniform(0.,1.) );
   double const r01l( r01 + rl );
   double const r01l_cb( r01l * r01l * r01l );
   double const r01_cb( r01 * r01 * r01 );

   double const diss_vec_length( cbrt( X * (r01l_cb - r01_cb) + r01_cb ) );

   Vector3 v;
   v[0] = 1.; v[1] = 1.; v[2] = 1.;

   Vector3 const unit_z( shape_.unit_z() );
   Vector3 const unit_x( normalize( subtract( v,
   multiply( unit_z, dot_product(v, unit_z) ) ) ) );
   Vector3 const unit_y( normalize( cross_product( unit_x, unit_z ) ) );

   double const x( diss_vec_length * sin( theta ) * cos( phi ) );
   double const y( diss_vec_length * sin( theta ) * sin( phi ) );
   double const z( diss_vec_length * cos( theta ) );

   position_pair_type pp01;

   // This is the particle that ends up on the structure
   pp01.first = subtract( op_surf, unit_z * (z * D_surf_D01) );

   // This is the particle that ends up in the bulk
   pp01.second = add( op_surf,
   add( unit_x * x,
   add( unit_y * y,
   unit_z * (z * D_bulk_D01) ) ) );

   return pp01;
   } */

   // Used by newBDPropagator
   double newBD_distance(const Vector3& new_pos, double radius, const Vector3& old_pos, double sigma) const override
   {
      UNUSED(radius,old_pos,sigma);
      return distance(new_pos);
   }

   //TODO what is that?
   //virtual double minimal_distance(double const& radius) const
   //{
   //double cylinder_radius = shape_.radius();
   //// Return minimal distance *to* surface.
   //return (cylinder_radius + radius) * traits_type::MINIMAL_SEPARATION_FACTOR - cylinder_radius;
   //}


   // Boundary condition handling
   // FIXME This is a mess but it works. See ParticleContainerBase.hpp for explanation.
   position_structid_pair apply_boundary(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const override
   {
      UNUSED(structure_container);
      return pos_struct_id;         //TODO  structure_container.apply_boundary(*this, pos_struct_id);
   }

   position_structid_pair cyclic_transpose(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const override
   {
      UNUSED(structure_container);
      return pos_struct_id;       // for now we do not support connected cylindrical surfaces.
   }

   /*
   // *** Dynamic dispatch for the structure functions *** //
   // *** 1 *** - One new position
   // This requires a double dynamic dispatch.
   // First dispatch
   virtual position_structid_pair get_pos_sid_pair(Structure const& target_structure, Vector3 const& position, double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   return target_structure.get_pos_sid_pair_helper(*this, position, offset, reaction_length, rng);
   }
   // Second dispatch
   virtual position_structid_pair get_pos_sid_pair_helper(CuboidalRegion const& origin_structure, Vector3 const& position, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<CuboidalRegion >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(SphericalSurface const& origin_structure, Vector3 const& position, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<SphericalSurface >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(CylindricalSurface const& origin_structure, Vector3 const& position, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<CylindricalSurface >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(DiskSurface const& origin_structure, Vector3 const& position, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<DiskSurface >(origin_structure, position, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_helper(PlanarSurface const& origin_structure, Vector3 const& position, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_helper_any<PlanarSurface >(origin_structure, position, offset, rl, rng);
   }
   // The template function that defines the actual final dispatch procedure.
   template<typename Tstruct_>
   position_structid_pair get_pos_sid_pair_helper_any(Tstruct_ const& origin_structure, Vector3 const& position, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   // redirect to structure function with well-defined typing
   return ::get_pos_sid_pair(origin_structure, *this, position, offset, rl, rng);
   };

   // *** 2 *** - Two new positions
   // Same principle as above, but different return type
   // First dispatch
   virtual position_structid_pair_pair get_pos_sid_pair_pair(Structure const& target_structure, Vector3 const& position, SpeciesType const& s1, SpeciesType const& s2, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   return target_structure.get_pos_sid_pair_pair_helper(*this, position, s1, s2, reaction_length, rng);
   }
   // Second dispatch
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(CuboidalRegion const& origin_structure, Vector3 const& position, SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<CuboidalRegion >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(SphericalSurface const& origin_structure, Vector3 const& position, SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<SphericalSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(CylindricalSurface const& origin_structure, Vector3 const& position, SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<CylindricalSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(DiskSurface const& origin_structure, Vector3 const& position, SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<DiskSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(PlanarSurface const& origin_structure, Vector3 const& position, SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_pair_helper_any<PlanarSurface >(origin_structure, position, s_orig, s_targ, rl, rng);
   }
   // The template function that defines the actual final dispatch procedure.
   template<typename Tstruct_>
   position_structid_pair_pair get_pos_sid_pair_pair_helper_any(Tstruct_ const& origin_structure, Vector3 const& position, SpeciesType const& s_orig, SpeciesType const& s_targ, double const& rl, RandomNumberGenerator& rng) const
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
   virtual position_structid_pair get_pos_sid_pair_2o(Structure const& origin_structure2, SpeciesTypeID const& target_sid, Vector3 const& CoM, double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   return origin_structure2.get_pos_sid_pair_2o_helper(*this, target_sid, CoM, offset, reaction_length, rng);
   }
   // Second dispatch
   virtual position_structid_pair get_pos_sid_pair_2o_helper(CuboidalRegion const& origin_structure1, SpeciesTypeID const& target_sid, Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<CuboidalRegion >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(SphericalSurface const& origin_structure1, SpeciesTypeID const& target_sid, Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<SphericalSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(CylindricalSurface const& origin_structure1, SpeciesTypeID const& target_sid, Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<CylindricalSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(DiskSurface const& origin_structure1, SpeciesTypeID const& target_sid, Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<DiskSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   virtual position_structid_pair get_pos_sid_pair_2o_helper(PlanarSurface const& origin_structure1, SpeciesTypeID const& target_sid, Vector3 const& CoM, double const& offset, double const& rl, RandomNumberGenerator& rng) const
   {
   return this->get_pos_sid_pair_2o_helper_any<PlanarSurface >(origin_structure1, target_sid, CoM, offset, rl, rng);
   }
   // The template function that defines the actual final dispatch procedure.
   template<typename Tstruct_>
   position_structid_pair get_pos_sid_pair_2o_helper_any(Tstruct_ const& origin_structure1, SpeciesTypeID const& target_sid, Vector3 const& CoM, double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   {
   // This method has to figure out where the product will be placed in case of a bimolecular reaction.
   // As a default, we place particles on the substructure or the lower-dimensional structure. If the structures
   // have the same structure type (=> same dimensionality) it does not matter on which structure we put the product,
   // as long as it has the structure type id of the product species. This is handled in cases '1' below.

   // 1 - Check whether one of the structures is the parent of the other. If yes, the daughter structure is the target.
   if (this->is_parent_of_or_has_same_sid_as(origin_structure1) && origin_structure1.has_valid_target_sid(target_sid) )
   // origin_structure1 is target
   return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);

   else if (origin_structure1.is_parent_of_or_has_same_sid_as(*this) && this->has_valid_target_sid(target_sid) )
   // this structure is target
   return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);

   // 2 - Check which structures has the lower dimensionality / particle degrees of freedom, and put the product there.
   else if (origin_structure1.shape().dof() < this->shape().dof() && origin_structure1.has_valid_target_sid(target_sid) )
   // origin_structure1 is target
   return ::get_pos_sid_pair(*this, origin_structure1, CoM, offset, reaction_length, rng);

   else if (this->shape().dof() < origin_structure1.shape().dof() && this->has_valid_target_sid(target_sid) )
   // this structure is target
   return ::get_pos_sid_pair(origin_structure1, *this, CoM, offset, reaction_length, rng);

   else throw propagation_error("Invalid target structure type: does not match product species structure type or has wrong hierarchy or dimensionality.");
   }

   //     // *** 4 *** - Generalized functions for pair reactions with two origin structures and one target structure
   //     // NOTE: This is yet unused, but possibly useful in the future.
   //     // Overloading get_pos_sid_pair again with signature (origin_structure2, target_structure, ...) and introducing
   //     // a triple dynamic dispatch.
   //     virtual position_structid_pair get_pos_sid_pair(Structure const& origin_structure2, Structure const& target_structure, Vector3 const& position,
   //                                                          double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const
   //     {
   //         return origin_structure2.get_pos_sid_pair_helper1(*this, target_structure, position, offset, reaction_length, rng);
   //     }


   */

   ///*** Formerly used functions of the Morelli scheme ***/
   //// DEPRECATED
   //virtual double drawR_gbd(double const& rnd, double const& r01, double const& dt, double const& D01, double const& v) const
   //{
   //    return drawR_gbd_1D(rnd, r01, dt, D01, v);
   //}
   //// DEPRECATED
   //virtual double p_acceptance(double const& k_a, double const& dt, double const& r01, Vector3 const& ipv, 
   //                            double const& D0, double const& D1, double const& v0, double const& v1) const
   //{
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
   //        if (ipv < 0 )
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
   //    return random_vector(drawR_gbd(rng.uniform(0., 1.), r01, dt, D01, v), rng);
   //}

   // TODO
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
