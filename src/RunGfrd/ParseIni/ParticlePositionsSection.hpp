#ifndef PARTICLEPOSITIONS_SECTION_HPP
#define PARTICLEPOSITIONS_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ParticlePositionSection final : SectionModeBase
{
   explicit ParticlePositionSection() : interval_(1E-3) { }

   ~ParticlePositionSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "ParticlePositions"; }

   const std::string key_mode = "Mode";
   const std::string key_interval = "Interval";
   const std::string key_file = "File";

   double interval() const { return interval_; }
   std::string file() const { return file_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key, value)) return true;
      if (key == key_interval) { interval_ = std::stod(value); return true; }
      if (key == key_file) { file_ = value; return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   double interval_;
   std::string file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ParticlePositionSection& cns)
{
   stream << "[" << cns.section_name() << "]" << std::endl;
   stream << cns.key_mode << " = " << (cns.mode() == SectionModeBase::modes::Run ? "Run" : cns.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   stream << cns.key_interval << " = " << cns.interval() << std::endl;
   stream << cns.key_file << " = " << cns.file() << std::endl;
   stream << std::endl;
   return stream;
}

#endif