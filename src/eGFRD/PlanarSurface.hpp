#ifndef PLANAR_SURFACE_HPP
#define PLANAR_SURFACE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "Surface.hpp"
#include "Plane.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class StructureContainer;

class PlanarSurface : public Surface<Plane>
{
public:
   typedef Surface<Plane>   base_type;

   typedef std::pair<Vector3, Vector3>                               position_pair_type;
   typedef std::pair<Vector3, StructureID>                           position_structid_pair;
   typedef std::pair<position_structid_pair, position_structid_pair>   position_structid_pair_pair;

   PlanarSurface(const std::string& name, const StructureTypeID sid, const StructureID pid, const shape_type& shape) noexcept
      : base_type(name, sid, pid, shape) {}

   const char* type_name() const override { return "PlanarSurface"; }

   // Create a random vector in the plane
   virtual Vector3 random_vector(double r, RandomNumberGenerator& rng) const override
   {
      return (shape_.unit_x() * rng.uniform(-1., 1.) + shape_.unit_y() * rng.uniform(-1., 1.)).normal() * r;
   }

   // Calculate a bd displacement for particles diffusing in the plane
   Vector3 bd_displacement(double mean, double r, RandomNumberGenerator& rng) const override
   {
      UNUSED(mean);
      double x = rng.normal(0., r), y = rng.normal(0., r);
      return shape_.unit_x() * x + shape_.unit_y() * y;
   }

   //// Rate for binding to particle on the structure
   //virtual double get_1D_rate_geminate(double k, double r01) const override
   //{
   //   return k / (2 * M_PI * r01);
   //}

   // Rate for binding to the structure
   double get_1D_rate_surface(double k, double r0) const override
   {
      UNUSED(r0);
      return k;
   }

   // Reaction volume for binding to particle in the structure
   virtual double particle_reaction_volume(double r01, double rl) const override
   {
      double r01l(r01 + rl);
      double r01l_sq(r01l * r01l);
      double r01_sq(r01 * r01);
      return M_PI * (r01l_sq - r01_sq);       // TODO Is there a correction factor if binding is allowed from both sides of the planes?
   }

   // Reaction volume for binding to the structure
   double surface_reaction_volume(double r0, double rl) const override
   {
      UNUSED(r0);
      return shape_.is_one_sided() ? rl : 2.0*rl;
   }

   // Vector of dissociation from the structure to parent structure
   //virtual Vector3 surface_dissociation_vector(RandomNumberGenerator& rng, double r0, double rl) const override
   //{
   //   double X(rng.uniform(0., 1.));
   //   double diss_vec_length(X*rl);

   //   if (shape_.is_one_sided())
   //   {
   //      return shape_.unit_z() * diss_vec_length;
   //   }
   //   else
   //   {
   //      double sign(rng.uniform_int(0, 1) * 2 - 1);
   //      diss_vec_length *= sign;
   //      return shape_.unit_z() * diss_vec_length;
   //   }
   //}

   //// Normed direction of dissociation from the structure to parent structure
   //virtual Vector3 surface_dissociation_unit_vector(RandomNumberGenerator& rng) const override
   //{
   //   return shape_.unit_z();
   //}

   //// Vector used to determine whether a particle has crossed the structure
   //// For the plane the normal vector is the natural choice
   virtual Vector3 side_comparison_vector() const override
   {
      return shape_.unit_z();
   }


   /* TODO TODO TODO TODO TODO TODO


   // Positions created at dissociation of one particle on the structure into two particles on the structure
   virtual position_pair_type geminate_dissociation_positions(RandomNumberGenerator& rng, SpeciesType const& s0, SpeciesType const& s1,
   Vector3 const& op, double const& rl) const
   {
   double const r01(s0.radius() + s1.radius());
   double const D01(s0.D() + s1.D());

   double X(rng.uniform(0., 1.));

   double const r01l(r01 + rl);
   double const r01l_sq(r01l * r01l);
   double const r01_sq(r01 * r01);

   // The square takes into account the radial character of the sampled length
   double const diss_vec_length(sqrt(r01_sq + X * (r01l_sq - r01_sq)));
   assert(diss_vec_length >= r01);

   Vector3 const m(random_vector(diss_vec_length, rng));

   return position_pair_type(op - m * s0.D() / D01,
   op + m * s1.D() / D01);
   }

   // Positions created at dissociation of one particle on the structure into two particles, one of which ends up in the bulk
   virtual position_pair_type special_geminate_dissociation_positions(RandomNumberGenerator& rng, SpeciesType const& s_surf, SpeciesType const& s_bulk,
   Vector3 const& op_surf, double const& rl) const
   {
   double const r01(s_bulk.radius() + s_surf.radius());
   double const D01(s_bulk.D() + s_surf.D());
   double const D_bulk_D01(s_bulk.D() / D01);
   double const D_surf_D01(s_surf.D() / D01);

   //double theta_max( M_PI/2 - asin(s_bulk.radius() / ( r01 )) );
   //theta_max = theta_max < 0 ? 0 : theta_max;

   // Create random angles for later construction of randomly oriented vector
   double const theta(rng.uniform(0., 1.) * M_PI);
   double const phi(rng.uniform(0., 1.) * 2 * M_PI);

   // Construct dissociation length
   double const X(rng.uniform(0., 1.));
   double const r01l(r01 + rl);
   double const r01l_cb(r01l * r01l * r01l);
   double const r01_cb(r01 * r01 * r01);

   double const diss_vec_length(cbrt(X * (r01l_cb - r01_cb) + r01_cb));


   // Determine direction of dissociation for the particle that ends up in the bulk
   // As a standard it is the plane's unit_z vector
   Vector3 unit_z(shape_.unit_z());
   // If the plane however is two-sided: randomize!
   if (!(shape_.is_one_sided()))
   unit_z = multiply(unit_z, rng.uniform_int(0, 1) * 2 - 1);

   // Construct randomly oriented dissociation vector
   double const x(diss_vec_length * sin(theta) * cos(phi));
   double const y(diss_vec_length * sin(theta) * sin(phi));
   double const z(diss_vec_length * cos(theta));

   position_pair_type pp01;

   // This is the particle that will end up on the surface
   pp01.first = subtract(op_surf,
   add(shape_.unit_x() * (x * D_surf_D01),
   shape_.unit_y() * (y * D_surf_D01)));

   // This is the particle that will end up in the bulk
   pp01.second = add(op_surf,
   add(shape_.unit_x() * (x * D_bulk_D01),
   add(shape_.unit_y() * (y * D_bulk_D01),
   unit_z * abs(z))));

   return pp01;

   }
*/
// Used by newBDPropagator
   double newBD_distance(const Vector3& new_pos, double radius, const Vector3& old_pos, double sigma) const override
   {
      UNUSED(radius);
      auto half_lengths(shape_.half_extent());
      auto new_pos_xyz(shape_.to_internal(new_pos));
      auto old_pos_xyz(shape_.to_internal(old_pos));
      if (new_pos_xyz.Z() * old_pos_xyz.Z() < 0 && ((abs(new_pos_xyz.X()) < half_lengths.X() && abs(new_pos_xyz.Y()) < half_lengths.Y()) ||
         (abs(old_pos_xyz.X()) < half_lengths.X() && abs(old_pos_xyz.Y()) < half_lengths.Y()))) // If the new and old positions lie on different sides of the plane
         return -1.0 * distance(new_pos) + sigma;
      return distance(new_pos) + sigma;
   }

   //virtual double minimal_distance(double const& radius) const
   //{
   //// PlanarSurface has thickness of 0.
   //return radius * traits_type::MINIMAL_SEPARATION_FACTOR;
   //}


   //  *** Boundary condition handling ***
   GFRD_EXPORT virtual position_structid_pair apply_boundary(const position_structid_pair& pos_structure_id, const StructureContainer& sc) const override;
   GFRD_EXPORT virtual position_structid_pair cyclic_transpose(const position_structid_pair& pos_structure_id, const StructureContainer& sc) const override;



   // *** Dynamic dispatch for the structure functions *** //
   // *** 1 *** - One new position
   // This requires a double dynamic dispatch.
   // First dispatch
   /*
   virtual position_structid_pair get_pos_sid_pair(Structure const& target_structure, Vector3 const& position,
   double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const override
   {
      return position_structid_pair(); //target_structure.get_pos_sid_pair_helper(*this, position, offset, reaction_length, rng);
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
   SpeciesType const& s1, SpeciesType const& s2, double const& reaction_length, RandomNumberGenerator& rng) const override
   {
   return position_structid_pair_pair();// target_structure.get_pos_sid_pair_pair_helper(*this, position, s1, s2, reaction_length, rng);
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
   Vector3 const& CoM, double const& offset, double const& reaction_length, RandomNumberGenerator& rng) const override
   {
   return position_structid_pair(); //origin_structure2.get_pos_sid_pair_2o_helper(*this, target_sid, CoM, offset, reaction_length, rng);
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


   */

   //     // *** 4 *** - Generalized functions for pair reactions with two origin structures and one target structure
   //     // NOTE: This is yet unused, but possibly useful in the future.
   //     // Overloading get_pos_sid_pair again with signature (origin_structure2, target_structure, ...) and introducing
   //     // a triple dynamic dispatch.
   //     virtual position_structid_pair get_pos_sid_pair(Structure const& origin_structure2, Structure const& target_structure, Vector3 const& position,
   //                                                          double const& offset, double const& reaction_length, RandomNumberGenerator const& rng) const
   //     {
   //         return origin_structure2.get_pos_sid_pair_helper1(*this, target_structure, position, offset, reaction_length, rng);
   //     }


   ///*** Formerly used functions of the Morelli scheme ***/
   //// DEPRECATED
   //virtual double drawR_gbd(double const& rnd, double const& r01, double const& dt, double const& D01, double const& v) const
   //{
   //    //TODO: use the 2D BD function instead of the 3D one - failed on very hard integral.
   //    return drawR_gbd_3D(rnd, r01, dt, D01);
   //}
   //// DEPRECATED
   //virtual double p_acceptance(double const& k_a, double const& dt, double const& r01, Vector3 const& ipv,
   //    double const& D0, double const& D1, double const& v0, double const& v1) const
   //{
   //    //TODO: use the 2D BD function instead of the 3D one. - Solution known
   //    return k_a * dt / ((I_bd_3D(r01, dt, D0) + I_bd_3D(r01, dt, D1)) * 4.0 * M_PI);
   //}
   //// DEPRECATED
   //virtual Vector3 dissociation_vector(RandomNumberGenerator& rng, double const& r01, double const& dt,
   //    double const& D01, double const& v) const
   //{
   //    return random_vector(drawR_gbd(rng(), r01, dt, D01, v), rng);
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

#endif /* PLANAR_SURFACE_HPP */
