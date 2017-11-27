#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ReactionRecordSection final : SectionModeBase
{
   explicit ReactionRecordSection() : SectionModeBase() { mode_ = modes::Run; }
   ~ReactionRecordSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "ReactionRecord"; }
   const std::string key_file = "File";
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
      std::cout  << std::setw(14) << "react_rec = " << "'" << (file_.empty() ? "stdout" : file_) << "'" << mode << "\n";
      std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   std::string file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ReactionRecordSection& rrs)
{
   stream << "[" << rrs.section_name() << "]" << std::endl;
   stream << rrs.key_mode << " = " << (rrs.mode() == SectionModeBase::modes::Run ? "Run" : rrs.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   stream << rrs.key_file << " = " << rrs.file() << std::endl;
   stream << std::endl;
   return stream;
}
