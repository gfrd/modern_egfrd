#ifndef SPHERE_HPP
#define SPHERE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "DefsEgfrd.hpp"
#include "Vector3.hpp"
#include <functional>
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Sphere
{
public:

   // constructor(s)

   Sphere() noexcept : Sphere(Vector3(), 0) {}

   Sphere(const Vector3& position, double radius) noexcept : position_(position), radius_(radius) { }

   // member functions

   const Vector3& position() const { return position_; }

   double radius() const { return radius_; }

   static int dof() { return 2; }        // degrees of freedom for particle movement

   bool operator==(const Sphere& rhs) const
   {
      return position_ == rhs.position_ && radius_ == rhs.radius_;
   }

   bool operator!=(const Sphere& rhs) const
   {
      return !operator==(rhs);
   }

   Sphere offset(Vector3 offset) const
   {
      return Sphere(position_ + offset, radius_);
   }

   // mutable methods

   Vector3& position() { return position_; }

   // some math utility functions

   double to_internal(const Vector3& pos) const
   {
      Vector3 intern(pos - position_);
      return intern.length();
   }

   std::pair<Vector3, std::pair<double, double>> project_point(const Vector3& pos) const
   {
      double r(to_internal(pos));
      // The projection of a point on a sphere is always the center point of the sphere.
      return std::make_pair(position_, std::make_pair(r, 0.0));
   }

   std::pair<Vector3, std::pair<double, double>> project_point_on_surface(const Vector3& pos) const
   {
      UNUSED(pos);
      // TODO. If we ever need it. The projection of a point on a sphere.
      return std::make_pair(Vector3(), std::make_pair(0.0, 0.0));
   }

   double distance(const Vector3& pos) const
   {
      return to_internal(pos) - radius_;
   }

   std::pair<Vector3, bool> deflect(const Vector3& r0, const Vector3& d) const
   {
      // Displacements are not deflected on spheres (yet),
      // but this function has to be defined for every shape to be used in structure.
      // For now it just returns the new position, no change
      return std::make_pair(r0 + d, false);
   }

   Vector3 random_position(RandomNumberGenerator& rng) const
   {
      // Random point in sphere!
      return position_ + rng.uniform(0, radius_) * Vector3::random(rng);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   Vector3 position_;
   double radius_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Sphere& s)
{
   stream << "Sphere{P=" << s.position() << ", R=" << s.radius() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Sphere >
   {
      std::size_t operator()(const Sphere& s) const
      {
         return hash<Vector3>()(s.position()) ^ hash<double>()(s.radius());
      }
   };
}

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* SPHERE_HPP */
