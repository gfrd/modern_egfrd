#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "exceptions.hpp"
#include "Vector3.hpp"
#include "DomainID.hpp"
#include "ShellID.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Shell
{
public:
   using shell_id_pair = std::pair<ShellID, Shell>;


   enum class Shape
   {
      SPHERE,
      CYLINDER,
   };

   enum class Code
   {
      INIT,                // to be initialized shell, always Spherical, with radius of particle
      NORMAL,              // constructed shell, of any kind, but no multi
      MULTI,               // shell is part of multi domain
   };

   explicit Shell() noexcept : did_(0), shape_(Shape::SPHERE), code_(Code::INIT), shell_(Sphere()) {}    // default constructor (use only for Persistence)

   explicit Shell(const DomainID did, const Sphere& sphere, Code code) noexcept : did_(did), shape_(Shape::SPHERE), code_(code), shell_(sphere) {}

   explicit Shell(const DomainID did, const Cylinder& cylinder, Code code) noexcept : did_(did), shape_(Shape::CYLINDER), code_(code), shell_(cylinder) {}

   DomainID did() const noexcept { return did_; }


   Shape shape() const noexcept { return shape_; }

   Code code() const noexcept { return code_; }

   const Vector3& position() const { return shape_ == Shape::SPHERE ? shell_.s.position() : shell_.c.position(); }

   const Sphere& get_sphere() const { THROW_UNLESS(illegal_state, shape_ == Shape::SPHERE); return shell_.s; }

   const Cylinder& get_cylinder() const { THROW_UNLESS(illegal_state, shape_ == Shape::CYLINDER); return shell_.c; }

   bool operator==(const Shell& rhs) const { return did_ == rhs.did_ && shape_ == rhs.shape_ && code_ == rhs.code_ && (shape_ == Shape::SPHERE ? shell_.s == rhs.get_sphere() : shape_ == Shape::CYLINDER ? shell_.c == rhs.get_cylinder() : true); }

   bool operator!=(const Shell& rhs) const { return !operator==(rhs); }


private:
   DomainID did_;
   Shape shape_;
   Code code_;
   union U
   {
      Sphere s;
      Cylinder c;

      explicit U(const Sphere& s) : s(s) {}
      explicit U(const Cylinder& c) : c(c) {}

   }  shell_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Shell& shell)
{
   stream << "Shell{" << shell.did() << ", Code=" << static_cast<std::underlying_type<Shell::Code>::type>(shell.code());
   if (shell.shape() == Shell::Shape::SPHERE) stream << shell.get_sphere();
   if (shell.shape() == Shell::Shape::CYLINDER) stream << shell.get_cylinder();
   stream << "}";
   return stream;
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Shell >
   {
      std::size_t operator()(const Shell& val)
      {
         std::size_t h = hash< DomainID >()(val.did()) ^ hash<std::underlying_type<Shell::Code>::type>()(static_cast<std::underlying_type<Shell::Code>::type>(val.code()));
         if (val.shape() == Shell::Shape::SPHERE) h ^= hash<Sphere>()(val.get_sphere());
         if (val.shape() == Shell::Shape::CYLINDER) h ^= hash<Cylinder>()(val.get_cylinder());
         return h;
      }
   };
}

// --------------------------------------------------------------------------------------------------------------------------------
