#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "DefsEgfrd.hpp"
#include "Vector3.hpp"
#include "Vector2.hpp"
#include "helperFunctions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Disk
{
public:

   // constructors

   Disk() noexcept : position_(), radius_(0), unitZ_(Vector3::uz) {}

   Disk(const Vector3& position, double radius, const Vector3& unitZ) noexcept : position_(position), radius_(radius), unitZ_(unitZ) {}

   // member functions

   const Vector3& position() const { return position_; }

   double radius() const { return radius_; }

   const Vector3& unit_z() const { return unitZ_; }

   static int dof() { return 0; }        // degrees of freedom for particle movement

   bool operator==(const Disk& rhs) const
   {
      return position_ == rhs.position_ && radius_ == rhs.radius_ && unitZ_ == rhs.unitZ_;
   }

   bool operator!=(const Disk& rhs) const
   {
      return !operator==(rhs);
   }

   // some math utility functions

   Vector2 to_internal(const Vector3& pos) const
   {
      const Vector3 intern(pos - position_);
      const double z(Vector3::dot(intern, unitZ_));               // z can be < 0
      const double r((intern - unitZ_* z).length());              // r is always >= 0
      return Vector2(r, z);
   }

   std::pair<Vector3, std::pair<double, double>>  project_point(const Vector3& pos) const
      // Calculates the projection of 'pos' onto the disk 'obj' and also returns the coefficient
      // for the normal component (z) of 'pos' in the basis of the disk and the distance of the
      // projected point to the 'edge' of the disk. Here a positive number means the projected
      // point is outside the disk, and a negative numbers means it is 'inside' the disk.
      // As a special case, we calculate the projection differently for a position that is in
      // the plane of the disk. Then we imagine that the particle is always 'above' the structure
      // and return -1 as a standard, instead of the distance to the disk center.
   {
      // Here we do not call 'to_internal' for efficiency
      const Vector3 intern(pos - position_);
      const double z(Vector3::dot(intern, unitZ_));
      const Vector3 r(intern - unitZ_ * z);
      const double rlength(r.length());

      // The quantities that will be return, with default values
      // for the standard case (pos is not in plane of disk)
      Vector3 proj_pos(position_ + r);
      double normal_comp(z);
      double dist_to_edge(rlength - radius_);

      // Special case: pos is in the same plane as the disk
      if (feq(z, 0.0, radius_)) // third argument is typical scale
      {
         proj_pos = position_; // projected position = disk center
         normal_comp = rlength;
         dist_to_edge = -1.0;
      }

      return std::make_pair(proj_pos, std::make_pair(normal_comp, dist_to_edge));
   }

   std::pair<Vector3, std::pair<double, double>> project_point_on_surface(const Vector3& pos) const
   {
      return project_point(pos);
   }

   double distance(const Vector3& pos) const
   {
      // First compute the (r,z) components of pos in a coordinate system
      // defined by the vectors unitR and unit_z, where unitR is
      // chosen such that unitR and unit_z define a plane in which pos lies.
      Vector2 rz(to_internal(pos));
      const double dz(std::abs(rz.Y()));
      const double dr(rz.X() - radius_);
      return dr > 0 ? std::sqrt(dz*dz + dr*dr) : dz;
   }

   std::pair<Vector3, bool> deflect(const Vector3& r0, const Vector3& d) const
   {
      // Displacements are not deflected on disks (yet),
      // but this function has to be defined for every shape to be used in structure.
      // For now it just returns original pos. + displacement. The change flag = false.
      return std::make_pair(r0 + d, false);
   }

   Vector3 random_position(RandomNumberGenerator& rng) const
   {
      UNUSED(rng);
      // The disk has only one "legal" position = its center
      return position_;
   }

private:
   Vector3   position_;    // centre.
   double    radius_;
   Vector3   unitZ_;      // Z-unit_z. should be normalized.
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Disk& d)
{
   stream << "Disk{P=" << d.position() << ", R=" << d.radius() << ", U=z" << d.unit_z() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Disk >
   {
      std::size_t operator()(const Disk& d) const
      {
         return hash<Vector3>()(d.position()) ^ hash<double>()(d.radius()) ^ hash<Vector3>()(d.unit_z());
      }
   };
}

// --------------------------------------------------------------------------------------------------------------------------------
