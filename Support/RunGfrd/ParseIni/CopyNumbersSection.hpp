#ifndef COPYNUMBERS_SECTION_HPP
#define COPYNUMBERS_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct CopyNumbersSection final : SectionBase
{
   explicit CopyNumbersSection() : mode_(modes::Off), type_(types::Instantaneous), interval_(1E-3) { }

   ~CopyNumbersSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   enum class modes { Off = 0, On = 1, Run = 2, };
   enum class types { Average = 0, Instantaneous = 1, };

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "CopyNumbers"; }

   const std::string key_mode = "Mode";
   const std::string key_interval = "Interval";
   const std::string key_file = "File";
   const std::string key_type = "Type";

   modes mode() const { return mode_; }
   double interval() const { return interval_; }
   std::string file() const { return file_; }
   types type() const { return type_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_keypair(const std::string& key, const std::string& value) override
   {
      if (key == key_mode) { mode_ = get_mode(value); return; }
      if (key == key_interval) { interval_ = std::stod(value); return; }
      if (key == key_file) { file_ = value; return; }
      if (key == key_type) { type_ = get_type(value); return; }
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

   static types get_type(const std::string& value)
   {
      std::string lcvalue;
      std::transform(value.begin(), value.end(), std::back_inserter(lcvalue), std::bind(std::tolower<char>, std::placeholders::_1, std::locale()));
      if (lcvalue == "average" || lcvalue == "0") return types::Average;
      if (lcvalue == "instantaneous" || lcvalue == "1") return types::Instantaneous;
      return types::Average;
   }

private:
   modes mode_;
   types type_;
   double interval_;
   std::string file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const CopyNumbersSection& cns)
{
   stream << "[" << cns.section_name() << "]" << std::endl;
   stream << cns.key_mode << " = " << (cns.mode() == CopyNumbersSection::modes::Run ? "Run" : cns.mode() == CopyNumbersSection::modes::On ? "On" : "Off") << std::endl;
   stream << cns.key_interval << " = " << cns.interval() << std::endl;
   stream << cns.key_file << " = " << cns.file() << std::endl;
   stream << cns.key_type << " = " << (cns.type() == CopyNumbersSection::types::Average ? "Average" : "Instantaneous") << std::endl;
   stream << std::endl;
   return stream;
}

#endif