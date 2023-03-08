#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SpeciesType.hpp"
#include "Model.hpp"
#include "SectionBase.hpp"
#include <iomanip>

// --------------------------------------------------------------------------------------------------------------------------------

struct ME_EXPORT SpeciesTypeSection final : SectionBase
{
   explicit SpeciesTypeSection() : SectionBase()
   {
      init_auto_vars( { { key_diffusion, 0.0}, { key_drift_velocity, 0.0}, { key_radius, 0.0} } );
   }

   ~SpeciesTypeSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "SpeciesType"; }

   const std::string key_name = "Name";
   const std::string& name() const { return name_; }
   const std::string key_diffusion = "D";
   double D() const { return auto_var_value(key_diffusion); }
   const std::string key_drift_velocity = "v";
   double v() const { return auto_var_value(key_drift_velocity); }
   const std::string key_radius = "r";
   double r() const { return auto_var_value(key_radius); }
   const std::string key_structure_type = "Structure";
   const std::string& structureType() const { return structure_type_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionBase::set_keypair(key,value)) return true;
      if (key == key_name)
      {
          THROW_UNLESS_MSG(illegal_section_value, name_.empty(), "Name already set! Use a new SpeciesType section to define a new species type.");
          name_ = format_check(value);
          return true;
      }
      if (key == key_structure_type)
      {
          THROW_UNLESS_MSG(illegal_section_value, structure_type_.empty(), "Structure type already set!");
          structure_type_ = value;
          return true;
      }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void create_species(Model& model, const VariablesSection& vars) const
   {
      auto sid = model.get_def_structure_type_id();

      if(!structure_type_.empty()) {
          sid = model.get_structure_type_id_by_name(structure_type_);
          if(sid == StructureTypeID(0)) {
              THROW_EXCEPTION(illegal_section_value, "Structure type '" << structure_type_ << "' was not previously created, please check the name.")
          }
      }

      for (auto name : names_)
         model.add_species_type(SpeciesType(name, sid, D(), r(), v()));
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      std::cout << std::setw(14) << "species = " << "'" << name_ << "'" << ", D = " << D() << " [m^2*s^-1], r = " << r() << " [m]";
      if (v() != 0) std::cout << ", v = " << v();
      std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   std::string format_check(std::string input)
   {
      auto regex = std::regex("[^\\s,]+");
      auto begin =  std::sregex_iterator(input.begin(), input.end(), regex);
      auto end = std::sregex_iterator();
 
      auto size = std::distance(begin, end);
      THROW_UNLESS_MSG(illegal_section_value, size > 0, "SpeciesTypeName not valid");
      names_.reserve(size);
    
      std::stringstream ss;
      for (auto i = begin; i != end; ++i) 
      {
         std::string name = i->str();
         if (!is_valid_speciestype_name(name)) THROW_EXCEPTION(illegal_section_value, "SpeciesTypeName '" << name << "'not valid");
         names_.emplace_back(name);
         if (i != begin) ss << ", ";
         ss << name;
      }

      return ss.str();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name_, structure_type_;
   std::vector<std::string> names_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const SpeciesTypeSection& sts)
{
   stream << "[" << sts.section_name() << "]" << std::endl;
   stream << sts.key_name << " = " << sts.name() << std::endl;
   stream << sts.key_radius << " = " << sts.r() << std::endl;
   stream << sts.key_diffusion << " = " << sts.D() << std::endl;
   if (sts.v() != 0) stream << sts.key_drift_velocity << " = " << sts.v() << std::endl;
   stream << sts.key_structure_type << " = " << sts.structureType() << std::endl;
   stream << std::endl;
   return stream;
}
