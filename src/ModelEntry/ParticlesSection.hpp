#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "randomNumberGenerator.hpp"
#include "World.hpp"
#include "SectionBase.hpp"

   // --------------------------------------------------------------------------------------------------------------------------------

struct ME_EXPORT ParticlesSection : SectionModeBase
{
   explicit ParticlesSection() : SectionModeBase() { mode_ = modes::On; }
   ~ParticlesSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "Particles"; }
   const std::vector<std::pair<const std::string, const int>>& particles() const { return particles_; }
   const std::vector<std::pair<const std::string, const std::vector<Vector3>>>& placement() const { return placement_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key, value)) return true;

      if (!is_valid_speciestype_name(key)) THROW_EXCEPTION(illegal_section_value, "SpeciesTypeName '" << key << "'not valid");

      std::smatch match;
      if (is_placement_string(value, match))
      {
         int N = 1;
         if (match.size() == 5 && match[4].matched)
         {
            auto d = vars_->evaluate_value_expression(match[4].str(), key);
            THROW_UNLESS_MSG(illegal_section_value, d >= 0, "Number of particles must be positive for '" << key << "'");
            N = static_cast<uint>(d);
         }

         auto positions = vars_->evaluate_value_expressions_repeat(match[1].str(), match[2].str(), match[3].str(), N, key);
         placement_.emplace_back(std::make_pair(key, std::move(positions)));
      }
      else
      {
         auto d = vars_->evaluate_value_expression(value, key);
         THROW_UNLESS_MSG(illegal_section_value, d >= 0, "Number of particles must be positive for '" << key << "'");
         particles_.emplace_back(std::make_pair(key, static_cast<uint>(d)));
      }

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_particles_to_world(const Model& model, World& world, RandomNumberGenerator& rng_) const
   {
      for (auto& ps_pair : particles_)
      {
         auto sid = model.get_species_type_id_by_name(ps_pair.first);
         world.throwInParticles(sid, ps_pair.second, rng_);
      }
      for (auto& pp_pair : placement_)
      {
         auto sid = model.get_species_type_id_by_name(pp_pair.first);
         auto species = model.get_species_type_by_id(sid);
         auto structureTypeId = species.structure_type_id();
         auto structids = world.get_structure_ids(structureTypeId);

         auto stid = world.get_def_structure_id();
         if(structids.size() == 1) {
             stid = *structids.begin();
         }

         for (auto &v : pp_pair.second)
            world.add_particle(sid, stid, v);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      std::string mode;
      if (mode_ != modes::On) mode = mode_ == modes::Off ? ", Mode = Off" : ", Mode = Run";
      for (auto& particles : particles_)
         std::cout << std::setw(14) << "particle(s) = " << "'" << particles.first << "'" << ", N = " << particles.second << ", pos=random" << mode << "\n";
      for (auto& placement : placement_)
      {
         if (placement.second.size() == 1)
            std::cout << std::setw(14) << "particle(s) = " << "'" << placement.first << "'" << ", N = 1, pos=" << placement.second[0] << mode << "\n";
         else
         {
            std::cout << std::setw(14) << "particle(s) = " << "'" << placement.first << "'" << ", N = " << placement.second.size() << mode << "\n";
            for (auto& v : placement.second)
               std::cout << std::setw(14) << " " << " pos=" << v << "\n";
         }
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   bool is_placement_string(const std::string& value, std::smatch& match)
   {
      auto regex = std::regex(R"(\(([^,]+),([^,]+),([^,]+)\)\s*(?:,\s*(.+))*)");
      return std::regex_search(value, match, regex);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::vector<std::pair<const std::string, const int>> particles_;
   std::vector<std::pair<const std::string, const std::vector<Vector3>>> placement_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ParticlesSection& ps)
{
   stream << "[" << ps.section_name() << "]" << std::endl;
   if (ps.mode() != SectionModeBase::modes::On) stream << ps.key_mode << " = " << (ps.mode() == SectionModeBase::modes::Run ? "Run" : ps.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   for (auto& particle : ps.particles()) { stream << particle.first + " = " << particle.second << std::endl; }
   //for (auto& placement : ps.placement()) { stream << placement.first + " = " << placement.second << std::endl; }
   stream << std::endl;
   return stream;
}
