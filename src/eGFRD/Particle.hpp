#ifndef PARTICLE_HPP
#define PARTICLE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include "DefsEgfrd.hpp"
#include "Vector3.hpp"
#include "Sphere.hpp"
#include "StructureID.hpp"
#include "SpeciesTypeID.hpp"
#include "SpeciesType.hpp"
#include "ParticleID.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct Particle final
{
   using particle_id_pair = std::pair<ParticleID, Particle>;

   // constructors

   Particle() noexcept : shape_(), species_id_(), structure_id_(), D_(0.), v_(0.) {}

   explicit Particle(const SpeciesTypeID species_id, const Sphere& shape, const StructureID structure_id, double D, double v = 0.0) noexcept
      : shape_(shape), species_id_(species_id), structure_id_(structure_id), D_(D), v_(v) {}

   explicit Particle(const SpeciesType& species, const StructureID structure_id, Vector3 position) noexcept         // construct particle from species type
      : shape_(Sphere(position, species.radius())), species_id_(species.id()), structure_id_(structure_id), D_(species.D()), v_(species.v()) {}

   // Get basic properties

   const Vector3& position() const { return shape_.position(); }

   double radius() const { return shape_.radius(); }

   double D() const { return D_; }

   double v() const { return v_; }

   const Sphere& shape() const { return shape_; }

   SpeciesTypeID sid() const { return species_id_; }

   // Get the id of the structure that the particle lives on.
   StructureID structure_id() const { return structure_id_; }

   bool operator==(const Particle& rhs) const
   {
      return species_id_ == rhs.sid() && shape_ == rhs.shape() && structure_id_ == rhs.structure_id();      // No D_ and v_ ???
   }

   bool operator!=(const Particle& rhs) const
   {
      return !operator==(rhs);
   }

   // mutable methods

   Vector3& position() { return shape_.position(); }           // NOTE: remember to update the MatrixSpace, when position is modified!! 

private:
   Sphere          shape_;
   SpeciesTypeID   species_id_;
   StructureID     structure_id_;
   double          D_;           // diffusion constant
   double          v_;           // drift v
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Particle& p)
{
   stream << "Particle{" << p.shape() << ", D=" << p.D() << ", v=" << p.v() << ", " << p.sid() << ", " << p.structure_id() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Particle >
   {
      std::size_t operator()(const Particle& val)
      {
         return hash<Vector3>()(val.position()) ^ hash<double>()(val.radius()) ^
            hash<double>()(val.D()) ^ hash<double>()(val.v()) ^
            hash<SpeciesTypeID>()(val.sid()) ^ hash<StructureID>()(val.structure_id());
      }
   };
} // namespace std

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* PARTICLE_HPP */
