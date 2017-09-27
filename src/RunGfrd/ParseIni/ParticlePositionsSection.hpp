#ifndef PARTICLEPOSITIONS_SECTION_HPP
#define PARTICLEPOSITIONS_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ParticlePositionSection final : SectionBase
{
   explicit ParticlePositionSection() : mode_(modes::Off), interval_(1E-3) { }

   ~ParticlePositionSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   enum class modes { Off = 0, On = 1, Run = 2, };

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "ParticlePositions"; }

   const std::string key_mode = "Mode";
   const std::string key_interval = "Interval";
   const std::string key_file = "File";

   modes mode() const { return mode_; }
   double interval() const { return interval_; }
   std::string file() const { return file_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_keypair(const std::string& key, const std::string& value) override
   {
      if (key == key_mode) { mode_ = get_mode(value); return; }
      if (key == key_interval) { interval_ = std::stod(value); return; }
      if (key == key_file) { file_ = value; return; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   static modes get_mode(const std::string& value)
   {
      std::string lcvalue;
      std::transform(value.begin(), value.end(), std::back_inserter(lcvalue), std::bind(std::tolower<char>, std::placeholders::_1, std::locale()));
      if (lcvalue == "off" || lcvalue == "0") return modes::Off;
      if (lcvalue == "on" || lcvalue == "yes" || lcvalue == "1") return modes::On;
      if (lcvalue == "run" || lcvalue == "2") return modes::Run;
      return modes::Off;
   }

private:
   modes mode_;
   double interval_;
   std::string file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ParticlePositionSection& cns)
{
   stream << "[" << cns.section_name() << "]" << std::endl;
   stream << cns.key_mode << " = " << (cns.mode() == ParticlePositionSection::modes::Run ? "Run" : cns.mode() == ParticlePositionSection::modes::On ? "On" : "Off") << std::endl;
   stream << cns.key_interval << " = " << cns.interval() << std::endl;
   stream << cns.key_file << " = " << cns.file() << std::endl;
   stream << std::endl;
   return stream;
}

#endif