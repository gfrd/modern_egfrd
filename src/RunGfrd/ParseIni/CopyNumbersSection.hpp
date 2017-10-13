#ifndef COPYNUMBERS_SECTION_HPP
#define COPYNUMBERS_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct CopyNumbersSection final : SectionModeBase
{
   explicit CopyNumbersSection() : type_(types::Instantaneous), interval_(1E-3) { }
   ~CopyNumbersSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   enum class types { Average = 0, Instantaneous = 1, };

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "CopyNumbers"; }

   const std::string key_interval = "Interval";
   const std::string key_file = "File";
   const std::string key_type = "Type";

   double interval() const { return interval_; }
   std::string file() const { return file_; }
   types type() const { return type_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key, value)) return true;
      if (key == key_interval) { interval_ = std::stod(value); return true; }
      if (key == key_file) { file_ = value; return true; }
      if (key == key_type) { type_ = get_type(value); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   static types get_type(const std::string& value)
   {
      std::string lcvalue;
      std::transform(value.begin(), value.end(), std::back_inserter(lcvalue), std::bind(std::tolower<char>, std::placeholders::_1, std::locale()));
      if (lcvalue == "average" || lcvalue == "0") return types::Average;
      if (lcvalue == "instantaneous" || lcvalue == "1") return types::Instantaneous;
      return types::Average;
   }

private:
   types type_;
   double interval_;
   std::string file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const CopyNumbersSection& cns)
{
   stream << "[" << cns.section_name() << "]" << std::endl;
   stream << cns.key_mode << " = " << (cns.mode() == SectionModeBase::modes::Run ? "Run" : cns.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   stream << cns.key_interval << " = " << cns.interval() << std::endl;
   stream << cns.key_file << " = " << cns.file() << std::endl;
   stream << cns.key_type << " = " << (cns.type() == CopyNumbersSection::types::Average ? "Average" : "Instantaneous") << std::endl;
   stream << std::endl;
   return stream;
}

#endif