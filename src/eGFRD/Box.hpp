#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "DefsEgfrd.hpp"
#include "Vector3.hpp"
#include "helperFunctions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Box
{
public:

   // constructor(s)

   Box() noexcept : Box(Vector3(), Vector3(0.5, 0.5, 0.5), Vector3::ux, Vector3::uy, Vector3::uz) {}

   Box(const Vector3& position) noexcept : Box(position, Vector3(0.5, 0.5, 0.5), Vector3::ux, Vector3::uy, Vector3::uz) {}

   Box(const Vector3& position, const Vector3& half_extent) noexcept : Box(position, half_extent, Vector3::ux, Vector3::uy, Vector3::uz) {}

   Box(const Vector3& position, const Vector3& half_extent, const Vector3& ux, const Vector3& uy, const Vector3& uz) noexcept : position_(position), half_extent_(half_extent), unitX_(ux), unitY_(uy), unitZ_(uz) { }

   // member functions

   const Vector3& position() const { return position_; }

   const Vector3& unit_x() const { return unitX_; }

   const Vector3& unit_y() const { return unitY_; }

   const Vector3& unit_z() const { return unitZ_; }

   const Vector3& half_extent() const { return half_extent_; }

   static int dof() { return 3; }     // degrees of freedom for particle movement

   bool operator==(const Box& rhs) const
   {
      return position_ == rhs.position_ && half_extent_ == rhs.half_extent_
         && unitX_ == rhs.unitX_ && unitY_ == rhs.unitY_ && unitZ_ == rhs.unitZ_;
   }

   bool operator!=(const Box& rhs) const
   {
      return !operator==(rhs);
   }

   // some math utility functions

   Vector3 to_internal(const Vector3& pos) const
   {
      Vector3 intern(pos - position_);
      //return Vector3(Vector3::dot(intern, obj.unitX_), Vector3::dot(intern, obj.unitY_), Vector3::dot(intern, obj.unitZ_));
      return intern;  // faster for orto-normal orientations
   }

   std::pair<Vector3, std::pair<double, double>> project_point(const Vector3& pos) const
   {
      Vector3 intern(to_internal(pos));
      Vector3 dxdydz(intern.abs() - half_extent_);
      const double min_dist(dxdydz.X() <= 0 && dxdydz.Y() <= 0 && dxdydz.Z() <= 0 ? std::max(dxdydz.X(), std::max(dxdydz.Y(), dxdydz.Z())) : 1.0);      // TODO make this is proper distance if we need it

      // TODO the projection of the point is not very well defined.
      // The projection of a point on a box.
      return std::make_pair(Vector3(), std::make_pair(double(), min_dist));
   }

   std::pair<Vector3, std::pair<double, double>> project_point_on_surface(const Vector3& pos) const
   {
      UNUSED(pos);
      // TODO. If we ever need it.
      // The projection of a point on a box.
      return std::make_pair(Vector3(), std::make_pair(double(), double()));
   }

   std::pair<Vector3, bool> deflect(const Vector3& r0, const Vector3& d) const
   {
      // Displacements are not deflected on cuboidal regions,
      // but this function has to be defined for every shape to be used in structure.
      // For now it just returns original pos. + displacement. The changeflage = false.
      return std::make_pair(r0 + d, false);
   }

   double distance(const Vector3& pos) const
   {
      Vector3 intern(to_internal(pos));
      Vector3 dxdydz(intern.abs() - half_extent_);
      Vector3 outside(std::max(0., dxdydz.X()), std::max(0., dxdydz.Y()), std::max(0., dxdydz.Z()));  // set negative values to zero, causes distance from box sides
      if (outside == Vector3::null) return std::max(std::max(dxdydz.X(), dxdydz.Y()), dxdydz.Z());    // when zero, its inside the box, return largest (negative) value, closed to a given side
      return outside.length();
   }

   Vector3 random_position(RandomNumberGenerator& rng) const
   {
      // -1 < rng() < 1. See for example CuboidalRegion.hpp.
      return position_ + Vector3(half_extent_.X() * rng.uniform(-1,1) , half_extent_.Y() * rng.uniform(-1, 1), half_extent_.Z() * rng.uniform(-1, 1));
   }


protected:

   Vector3 position_;                       // vector to middle of box.
   Vector3 half_extent_;                    // Extent: for a box of 2 by 2 by 2, half_extent is 1 by 1 by 1.
   Vector3 unitX_, unitY_, unitZ_;          // 3 unit vectors defining orientation in space, for now it always equals the ortho-normal base (x,y,z)
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Box& b)
{
   stream << "Box{P=" << b.position() << ", HE=" << b.half_extent() << ", U=x" << b.unit_x() << ",y" << b.unit_y() << ",z" << b.unit_z() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Box >
   {
      std::size_t operator()(const Box& b) const
      {
         return hash<Vector3>()(b.position()) ^
            hash<Vector3>()(b.unit_x()) ^ hash<Vector3>()(b.unit_y()) ^ hash<Vector3>()(b.unit_z()) ^
            hash<Vector3>()(b.half_extent());
      }
   };
}

// --------------------------------------------------------------------------------------------------------------------------------
