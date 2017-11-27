#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"

   // --------------------------------------------------------------------------------------------------------------------------------

struct ProgressSection : SectionModeBase
{
   explicit ProgressSection() : SectionModeBase(), column_width_(80) { mode_ = modes::On; }
   ~ProgressSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "Progress"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   const std::string key_width = "Width";

   uint column_width() const { return column_width_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key, value)) return true;
      if (key == key_width) { column_width_ = static_cast<uint>(vars_->evaluate_value_expression(value, key_width)); return true; }
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      // don't print (progress will reveal itself when enables)
      //std::cout << "Progress = " << (mode_ == modes::Off ? "Off" : (mode_ == modes::Run ? "Run" : "On"));
      //if (column_width_ != 80) std::cout << std::setw(14) << "Width = " << column_width_ << "\n";
      //std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   uint column_width_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ProgressSection& ps)
{
   stream << "[" << ps.section_name() << "]" << std::endl;
   if (ps.mode() != SectionModeBase::modes::On) stream << ps.key_mode << " = " << (ps.mode() == SectionModeBase::modes::Run ? "Run" : ps.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   stream << "Width" " = " << ps.column_width() << std::endl;
   stream << std::endl;
   return stream;
}
