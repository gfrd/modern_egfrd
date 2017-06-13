#ifndef PARTICLE_SECTION_HPP
#define PARTICLE_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SpeciesType.hpp"
#include "SectionBase.hpp"

   // --------------------------------------------------------------------------------------------------------------------------------

struct SpeciesTypeSection final : SectionBase
{
   explicit SpeciesTypeSection() : name_(std::string()), D_(0), v_(0), r_(0) { }
   ~SpeciesTypeSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "SpeciesType"; }

   const std::string key_name = "Name";
   const std::string& name() const { return name_; }
   const std::string key_diffusion = "D";
   double D() const { return D_; }
   const std::string key_drift_velocity = "v";
   double v() const { return v_; }
   const std::string key_radius = "r";
   double r() const { return r_; }
   const std::string key_structure_id = "StructureType";

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_keypair(const std::string& key, const std::string& value) override
   {
      if (key == key_name) { name_ = value; return; }
      if (key == key_radius) { r_ = std::stod(value); return; }
      if (key == key_drift_velocity) { v_ = std::stod(value); return; }
      if (key == key_diffusion) { D_ = std::stod(value); return; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   SpeciesType create_species(const StructureTypeID& structure_type) const
   {
      return SpeciesType(name_, structure_type, D_, r_, v_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   std::string name_;
   double D_;
   double v_;
   double r_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const SpeciesTypeSection& sts)
{
   stream << "[" << sts.section_name() << "]" << std::endl;
   stream << sts.key_name << " = " << sts.name() << std::endl;
   stream << sts.key_radius << " = " << sts.r() << std::endl;
   stream << sts.key_diffusion << " = " << sts.D() << std::endl;
   if (sts.v()) stream << sts.key_drift_velocity << " = " << sts.v() << std::endl;
   if (false) stream << sts.key_structure_id << " = " << "world" << std::endl;
   stream << std::endl;
   return stream;
}

#endif