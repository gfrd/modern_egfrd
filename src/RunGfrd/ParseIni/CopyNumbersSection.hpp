#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"
#include <iomanip>

// --------------------------------------------------------------------------------------------------------------------------------

struct CopyNumbersSection final : SectionModeBase
{
   explicit CopyNumbersSection() : SectionModeBase(), type_(types::Instantaneous) { mode_ = modes::Run; init_auto_vars({ {key_interval, 1E-3 }}); }
   
   ~CopyNumbersSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   enum class types { Average = 0, Instantaneous = 1, };

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "CopyNumbers"; }

   const std::string key_interval = "Interval";
   const std::string key_file = "File";
   const std::string key_type = "Type";

   double interval() const { return auto_var_value(key_interval); }
   std::string file() const { return file_; }
   types type() const { return type_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key, value)) return true;
      if (key == key_file) { file_ = vars_->evaluate_string_expression(value, key_file); return true; }
      if (key == key_type) { type_ = get_type(value); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      std::string mode;
      if (mode_ != modes::Run) mode = mode_ == modes::Off ? ", Mode = Off" : ", Mode = On";
      std::cout  << std::setw(14) << "copynumbers = " << "'" << (file_.empty() ? "stdout" : file_) << "'" << ", Interval = " << interval()  << " [s], Type = " << (type_ == types::Average ? "Average" : "Instantaneous")  << mode << "\n";
      std::cout << "\n";
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
   std::string file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const CopyNumbersSection& cns)
{
   stream << "[" << cns.section_name() << "]" << std::endl;
   stream << cns.key_mode << " = " << (cns.mode() == SectionModeBase::modes::Run ? "Run" : cns.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   stream << cns.key_interval << " = " << cns.interval() << std::endl;
   if (!cns.file().empty()) stream << cns.key_file << " = " << cns.file() << std::endl;
   stream << cns.key_type << " = " << (cns.type() == CopyNumbersSection::types::Average ? "Average" : "Instantaneous") << std::endl;
   stream << std::endl;
   return stream;
}
