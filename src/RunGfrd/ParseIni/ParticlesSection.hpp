#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "randomNumberGenerator.hpp"
#include "World.hpp"
#include "SectionBase.hpp"

   // --------------------------------------------------------------------------------------------------------------------------------

struct ParticlesSection : SectionModeBase
{
   explicit ParticlesSection(): SectionModeBase() { mode_ = modes::On; }
   ~ParticlesSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "Particles"; }
   const std::vector<std::pair<const std::string, const int>>& particles() const { return particles_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key, value)) return true;

      if (!is_valid_speciestype_name(key)) THROW_EXCEPTION(illegal_section_value, "SpeciesTypeName '" << key << "'not valid");
      auto d = vars_->evaluate_value_expression(value, key);
      THROW_UNLESS_MSG(illegal_section_value, d >= 0, "Number of particles must be positive for '" << key << "'");
      particles_.emplace_back(std::make_pair(key, static_cast<uint>(d)));
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_particles_to_world(const Model& model, World& world, RandomNumberGenerator& rng_) const
   {
      for (auto ps_pair : particles_)
      {
         auto sid = model.get_species_type_id_by_name(ps_pair.first);
         world.throwInParticles(sid, ps_pair.second, rng_);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      std::string mode;
      if (mode_ != modes::On) mode = mode_ == modes::Off ? ", Mode = Off" : ", Mode = Run";
         for (auto& particles : particles_)
            std::cout << std::setw(14) << "particle = " << "'" << particles.first << "'" << ", N = " << particles.second << mode << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   std::vector<std::pair<const std::string, const int>> particles_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ParticlesSection& ps)
{
   stream << "[" << ps.section_name() << "]" << std::endl;
   if (ps.mode() != SectionModeBase::modes::On) stream << ps.key_mode << " = " << (ps.mode() == SectionModeBase::modes::Run ? "Run" : ps.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   for (auto particle : ps.particles()) { stream << particle.first + " = " << particle.second << std::endl; }
   stream << std::endl;
   return stream;
}
