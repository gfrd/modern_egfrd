#ifndef THROWINPARTICLES_SECTION_HPP
#define THROWINPARTICLES_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "randomNumberGenerator.hpp"
#include "World.hpp"
#include "SectionBase.hpp"

   // --------------------------------------------------------------------------------------------------------------------------------

struct ParticlesSection : SectionBase
{
   explicit ParticlesSection() { }
   ~ParticlesSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "Particles"; }
   const std::vector<std::pair<const std::string, const int>>& particles() const { return particles_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_keypair(const std::string& key, const std::string& value) override
   {
      particles_.emplace_back(std::make_pair(key, std::stoi(value)));
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

private:
   std::vector<std::pair<const std::string, const int>> particles_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ParticlesSection& ps)
{
   stream << "[" << ps.section_name() << "]" << std::endl;
   for (auto particle : ps.particles()) { stream << particle.first + " = " << particle.second << std::endl; }
   stream << std::endl;
   return stream;
}

#endif
