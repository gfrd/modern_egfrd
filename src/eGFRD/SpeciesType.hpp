#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "SpeciesTypeID.hpp"
#include "exceptions.hpp"
#include "makeString.hpp"
#include "StructureTypeID.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class SpeciesType
{
public:
   SpeciesType() = default;                                 // use default construct
   SpeciesType(const SpeciesType&) = default;               // copy operation
   SpeciesType& operator=(const SpeciesType&) = default;
   SpeciesType(SpeciesType&&) = default;                    // move operation
   SpeciesType& operator=(SpeciesType&&) = default;

   explicit SpeciesType(const std::string& name, const StructureTypeID structure, double D = 0., double r = 0., double v = 0.) noexcept
      : name_(name), diffusion_coef_(D), drift_velocity_(v), radius_(r), structure_type_id_(structure) {}

   SpeciesTypeID id() const
   {
      if (!id_) throw illegal_state("Not bound to Model");
      return id_;
   };

   // Get the name
   const std::string & name() const { return name_; }

   // Get the particle radius
   double radius() const { return radius_; }

   // Get the id of the structure type the species lives on
   StructureTypeID structure_type_id() const { return structure_type_id_; }

   // Get the diffusion constant
   double D() const { return diffusion_coef_; }

   // Get the drift
   double v() const { return drift_velocity_; }

   // Check equality/inequality
   bool operator==(const SpeciesType& rhs) const
   {
      return id_ == rhs.id() && diffusion_coef_ == rhs.D() && drift_velocity_ == rhs.v() &&
         radius_ == rhs.radius() && structure_type_id_ == rhs.structure_type_id();
   }

   bool operator!=(const SpeciesType& rhs) const { return !operator==(rhs); }

protected:
   friend class Model;  // class that can set id
   void set_id(const SpeciesTypeID id) { id_ = id; }

   friend class Persistence;

private:
   SpeciesTypeID  id_;
   std::string    name_;
   double         diffusion_coef_;
   double         drift_velocity_;
   double         radius_;
   StructureTypeID  structure_type_id_;      // The structure type that the particle lives on, default means world

};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const SpeciesType& s)
{
   stream << "SpeciesType{" << s.id() << " '" << s.name() << "', D=" << s.D() << ", v=" << s.v() << ", r=" << s.radius() << ", " << s.structure_type_id() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

