#ifndef CYLINDER_HPP
#define CYLINDER_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <algorithm>
#include "Vector2.hpp"
#include "Vector3.hpp"
#include <functional>

// --------------------------------------------------------------------------------------------------------------------------------

class Cylinder
{
public:

   // constructors

   Cylinder() noexcept : Cylinder(Vector3(), double(1.0), Vector3::uz, double(0.5)) {}

   Cylinder(const Vector3& position, double radius, const Vector3& uz, double half_length) noexcept : position_(position), radius_(radius), half_length_(half_length), unitZ_(uz) { }

   // member functions

   const Vector3& position() const { return position_; }

   double radius() const { return radius_; }

   const Vector3& unit_z() const { return unitZ_; }

   double half_length() const { return half_length_; }

   static int dof() { return 1; }     // degrees of freedom for particle movement

   bool operator==(const Cylinder& rhs) const
   {
      return position_ == rhs.position_ && half_length_ == rhs.half_length_ && radius_ == rhs.radius_ && unitZ_ == rhs.unitZ_;
   }

   bool operator!=(const Cylinder& rhs) const
   {
      return !operator==(rhs);
   }

   // some math utility functions

   Cylinder offset(Vector3 offset) const
   {
      return Cylinder(position_ + offset, radius_, unit_z(), half_length_);
   }


   Vector2 to_internal(const Vector3& pos) const
   {
      const Vector3 intern(pos - position_);
      const double z(Vector3::dot(intern, unitZ_));                 // z can be < 0
      const double r((intern - (unitZ_ * z)).length());              // r is always >= 0
      return Vector2(r, z);                                          // z = internal distance along cyl, r = internal distance from center line
   }

   // The following projects 'pos' on the cylinder
   // It returns a pair of which the first component is the projected position.
   // The second is again a pair of which the first entry is the distance of 'pos'
   // from the cylinder axis, the second a length l which indicates whether the
   // projected position is in the cylinder (true if l negative; l is the negative
   // distance of the projected position to the cylinder edge).
   std::pair<Vector3, std::pair<double, double>> project_point(const Vector3& pos) const
   {
      // The projection lies on the z-axis.
      Vector2 rz(to_internal(pos));
      return std::make_pair(position_ + unitZ_ * rz.X(), std::make_pair(rz.X(), std::abs(rz.Y()) - half_length_));
   }


   //Almost equal to projected point method, but for the subtraction of the cylinder radius from the radial distance r.
   //And projected point now lies on the surface, not on the central axis.
   std::pair<Vector3, std::pair<double, double>>  project_point_on_surface(const Vector3& pos) const
   {
      // Here we do not call 'to_internal' for efficiency
      const Vector3 intern(pos - position_);
      const double zlength(Vector3::dot(intern, unitZ_));
      const Vector3 z(unitZ_ * zlength);
      const Vector3 r(intern - z);
      const double rlength(r.length());
      const Vector3 projected_point(position_ + z);
      return std::make_pair(projected_point + r.normal() * radius_, std::make_pair(rlength - radius_, std::abs(zlength) - half_length_));
   }

   double distance(const Vector3& pos) const
   {
      // First compute the (r,z) components of pos in a coordinate system
      // defined by the vectors unitR and unit_z, where unitR is
      // chosen such that unitR and unit_z define a plane in which pos lies.

      Vector2 rz(to_internal(pos));
      const double dz(std::abs(rz.Y()) - half_length_);
      const double dr(rz.X() - radius_);

      if (dz > 0) // pos is (either) to the right or to the left of the cylinder.
         return dr > 0.0 ? std::sqrt(dz * dz + dr * dr) : dz;
      return dr >= 0 ? dr : std::max(dr, dz);
   }

   std::pair<Vector3, bool> deflect(const Vector3& r0, const Vector3& d) const
   {
      // Displacements are not deflected on cylinders (yet),
      // but this function has to be defined for every shape to be used in structure.
      // For now it just returns original pos. + displacement. The changeflage = false.
      return std::make_pair(r0 + d, false);
   }

   Vector3 random_position(RandomNumberGenerator& rng) const
   {
      // -1 < rng() < 1. See for example CylindricalSurface.hpp.
      return position_ + unitZ_ * rng.uniform(-1, 1) * half_length_;
   }

private:
   Vector3 position_;
   double radius_;
   double half_length_;
   Vector3 unitZ_;                  // Z unit vectors defining orientation in space, should be normalized, and for now it always one of the ortho-normal base (x,y,z)
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Cylinder& c)
{
   stream << "Cylinder{P=" << c.position() << ", R=" << c.radius() << ", HL=" << c.half_length() << ", U=z" << c.unit_z() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Cylinder >
   {
      std::size_t operator()(const Cylinder& c) const
      {
         return hash<Vector3>()(c.position()) ^ hash<double>()(c.radius()) ^
            hash<Vector3>()(c.unit_z()) ^ hash<double>()(c.half_length());
      }
   };
} // namespace std

#endif /* CYLINDER_HPP */
