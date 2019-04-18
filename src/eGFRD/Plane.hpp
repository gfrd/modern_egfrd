#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <algorithm>
#include "DefsEgfrd.hpp"
#include "exceptions.hpp"
#include "Vector3.hpp"
#include "Vector2.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Plane
{
public:

   // constructor(s)

   Plane() noexcept : Plane(Vector3(), Vector3::ux, Vector3::uy, 0.5, 0.5, false) {}

   Plane(const Vector3& position, const Vector3& vx, const Vector3& vy, double half_lx, double half_ly, bool is_one_sided)
      : position_(position), half_extent_(Vector2(half_lx, half_ly)), unitX_(vx.normal()), unitY_(vy.normal()), unitZ_(Vector3::cross(vx, vy).normal()), is_one_sided_(is_one_sided)
   {
      auto a = std::abs(Vector3::dot(unitX_, unitY_));
      THROW_UNLESS(illegal_argument, a < 1E-8);    // unitX,unitY should be at 90' angle
   }

   // member functions

   const Vector3& position() const { return position_; }

   const Vector3& unit_x() const { return unitX_; }

   const Vector3& unit_y() const { return unitY_; }

   const Vector3& unit_z() const { return unitZ_; }

   const Vector2& half_extent() const { return half_extent_; }

   static int dof() { return 2; }     // degrees of freedom for particle movement

   bool is_one_sided() const { return is_one_sided_; }

   bool operator==(const Plane& rhs) const
   {
      return position_ == rhs.position_ && half_extent_ == rhs.half_extent_
         && unitX_ == rhs.unitX_ && unitY_ == rhs.unitY_ && unitZ_ == rhs.unitZ_ && is_one_sided_ == rhs.is_one_sided_;
   }

   bool operator!=(const Plane& rhs) const
   {
      return !operator==(rhs);
   }

   // some math utility functions

   Vector3 to_internal(const Vector3& pos) const
      // The function calculates the coefficients to express 'pos' into the base of the plane 'obj'
   {
      Vector3 intern(pos - position_);
      return Vector3(Vector3::dot(intern, unitX_), Vector3::dot(intern, unitY_), Vector3::dot(intern, unitZ_));
   }

   std::pair<Vector3, std::pair<double, double>> project_point(const Vector3& pos) const
      // Calculates the projection of 'pos' onto the plane 'obj' and also returns the coefficient
      // for the normal component (z) of 'pos' in the basis of the plane (pair entry .second.first)
      // and the negative distance to the closest edge (pair entry .second.second). Note that the
      // latter also indicates whether the projected point is in the plane.
   {

      Vector3 intern(to_internal(pos));
      const double dx(abs(intern.X()) - half_extent_.X());
      const double dy(abs(intern.Y()) - half_extent_.Y());
      const double min_dist(dx <= 0 && dy <= 0 ? std::max(dx, dy) : 1.0);      // TODO make this is proper distance if we need it
      return std::make_pair(position_ + unitX_ * intern.X() + unitY_ * intern.Y(), std::make_pair(intern.Z(), min_dist));
   }

   std::pair<Vector3, std::pair<double, double>> project_point_on_surface(const Vector3& pos) const
      // Since the projected point on the plane, is already on the surface of the plane,
      // this function is just a wrapper of projected point.
   {
      return project_point(pos);
   }

   double distance(const Vector3& pos) const
      // Calculates the distance from 'pos' to plane 'obj' Note that when the plane is finite,
      // it also calculates the distance to the edge of the plane if necessary
   {
      Vector3 intern(to_internal(pos));

      // The (negative) distances to the plane edges (in x/y direction)
      double const dx(fabs(intern.X()) - half_extent_.X());
      double const dy(fabs(intern.Y()) - half_extent_.Y());
      // The z-distance (may be positive or negative, depending on the side)
      double const dz(intern.Z());

      if (dx < 0 && dy < 0)
      {
         // pos is positioned over the plane (projected point is in the plane, not next to it).
         return fabs(dz);
      }

      if (dx > 0) // outside the plane in the x direction
         return std::sqrt((dy > 0 ? std::pow(dx,2) + std::pow(dy,2) : std::pow(dx,2)) + std::pow(dz,2));
      // inside the plane in x direction
      return dy > 0 ? std::sqrt(std::pow(dy,2) + std::pow(dz,2)) : fabs(dz);
   }

   std::pair<Vector3, bool> deflect(const Vector3& r0, const Vector3& d) const
      // This routine deflects a displacement on a plane.
      // If the displacement vector intersects with the plane (starting from r0) the part
      // ranging out of the plane will be deflected into the plane, always into the direction towards the plane center.
      // ATTENTION! As for now, this only works for displacement vectors perpendicular to the plane
   {
      double intersect_parameter;

      // Calculate the intersection parameter and intersection point.
      // r0 is the origin position, d is a displacement vector.
      // If intersect_parameter <= 1 we have an intersection.
      if (Vector3::dot(d, unitZ_) != 0.0)
         intersect_parameter = Vector3::dot(position_ - r0, unitZ_) / Vector3::dot(d, unitZ_);
      else
         intersect_parameter = 1e100;       // infinity, displacement is parallel to edge

      // Check whether the displacement actually crosses the plane;
      // If not, just return original position plus displacement;
      // if yes, calculate the deflection.
      if (intersect_parameter > 1.0 || intersect_parameter < 0.0)
      {
         // No intersection; the new position is just r0+displacement
         return std::make_pair(r0 + d, false);
      }

      // Calculate the intersection point and the part of the displacement
         // that ranges out of the target plane.
         // Project all points that are supposed to be in the plane into it,
         // just to be sure.
         Vector3 intersect_pt = project_point(r0 + d * intersect_parameter).first;
         Vector3 d_out = d * (1.0 - intersect_parameter);

         // Calculate the length of the component of d_out perpendicular to the edge
         // and the vector of the component parallel to the edge
         double l_perp = Vector3::dot(d_out, unitZ_);       // note that this can be positive or negative,
         // depending on whether u_z points in or out!
         Vector3 d_edge = d_out - unitZ_ * l_perp;
         Vector3 d_edge_n = d_edge.normal();

         // Find the vector pointing from the edge towards the center of the plane, which is
         // the component of (center_pt - intersect_pt) perpendicular to the edge.
         // First we calculate the component parallel to the edge
         Vector3 icv = position_ - intersect_pt;
         Vector3 icv_edge = d_edge_n * Vector3::dot(icv, d_edge_n);
         Vector3 icv_perp = icv - icv_edge;

         // Calculate the component perpendicular to the edge in the plane.
         // Note that l_perp can be pos. or neg. while icv_perp always points towards
         // the center of the plane; therefore use abs(l_perp).
         Vector3 d_perp = icv_perp.normal() * abs(l_perp);

         // Construct the new position vector, make sure it's in the plane to
         // avoid trouble with periodic boundary conditions
         Vector3 new_pos = project_point(intersect_pt + (d_edge + d_perp)).first;

      // for now this returns the new position without changes
      return std::make_pair(new_pos, true);
   }

   Vector3 random_position(RandomNumberGenerator& rng) const
   {
      // -1 < rng() < 1. See for example PlanarSurface.hpp.
      return position_ + unitX_ * half_extent_.X() * rng.uniform(-1, 1) + unitY_ * half_extent_.Y() * rng.uniform(-1, 1);
   }

protected:

   Vector3 position_;
   Vector2 half_extent_;
   Vector3 unitX_, unitY_, unitZ_;
   bool is_one_sided_;
};


// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Plane& p)
{
   stream << "Plane{P=" << p.position() << ", HE=" << p.half_extent() << ", U=x" << p.unit_x() << ",y" << p.unit_y() << ",z" << p.unit_z() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Plane >
   {
      std::size_t operator()(const Plane& p) const
      {
         return hash<Vector3>()(p.position()) ^
            hash<Vector3>()(p.unit_x()) ^ hash<Vector3>()(p.unit_y()) ^ hash<Vector3>()(p.unit_z()) ^
            hash<Vector2>()(p.half_extent());
      }
   };
}

// --------------------------------------------------------------------------------------------------------------------------------
