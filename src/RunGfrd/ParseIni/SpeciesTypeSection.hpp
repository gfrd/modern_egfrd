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

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (key == key_name) { name_ = format_check(value); return true; }
      if (key == key_radius) { r_ = std::stod(value); return true; }
      if (key == key_drift_velocity) { v_ = std::stod(value); return true; }
      if (key == key_diffusion) { D_ = std::stod(value); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void create_species(Model& model) const
   {
      auto sid = model.get_def_structure_type_id();
      
      for (auto name : names_)
         model.add_species_type(SpeciesType(name, sid, D_, r_, v_));
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
         std::smatch match;
         auto regex2 = std::regex("[a-zA-Z][\\w\\-*_']*");
         if (!std::regex_match(name, match, regex2)) THROW_EXCEPTION(illegal_section_value, "SpeciesTypeName '" << name << "'not valid");
         names_.emplace_back(name);
         if (i != begin) ss << ", ";
         ss << name;
      }

      return ss.str();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name_;
   std::vector<std::string> names_;
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