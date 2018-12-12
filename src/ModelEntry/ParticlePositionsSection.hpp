#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ME_EXPORT ParticlePositionSection final : SectionModeBase
{
   explicit ParticlePositionSection() : SectionModeBase() { mode_ = modes::Run; init_auto_vars({ {key_interval, 1E-3 }}); }

   ~ParticlePositionSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "ParticlePositions"; }

   const std::string key_mode = "Mode";
   const std::string key_interval = "Interval";
   const std::string key_file = "File";

   double interval() const { return auto_var_value(key_interval); }
   std::string file() const { return file_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key, value)) return true;
      if (key == key_file) { file_ = vars_->evaluate_string_expression(value, key_file); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      std::string mode;
      if (mode_ != modes::Run) mode = mode_ == modes::Off ? ", Mode = Off" : ", Mode = On";
      std::cout  << std::setw(14) << "pos_record = " << "'" << (file_.empty() ? "stdout" : file_) << "'" << ", Interval = " << interval()  << " [s]"  << mode << "\n";
      std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
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
