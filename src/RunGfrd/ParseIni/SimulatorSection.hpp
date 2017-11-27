#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct SimulatorSection final : SectionBase
{
   explicit SimulatorSection() : SectionBase(), seed_(0), maintenance_step_(0)
   {
      init_auto_vars( { {key_prepare_time, 0.0}, {key_end_time, 0.0},  } );
   }

   ~SimulatorSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "Simulator"; }

   const std::string key_seed = "Seed";
   const std::string key_prepare_time = "PrepareTime";
   const std::string key_end_time = "EndTime";
   const std::string key_maintenance_step = "MaintenanceStep";
   const std::string key_maintenance_file = "MaintenanceFile";

   uint seed() const { return seed_; }
   double prepare_time() const { return auto_var_value(key_prepare_time); }
   double end_time() const { return auto_var_value(key_end_time); }
   size_t maintenance_step() const { return maintenance_step_; }
   std::string maintenance_file() const { return simstate_file_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionBase::set_keypair(key,value)) return true;
      if (key == key_seed) { seed_ = parse_seed(value); return true; }
      if (key == key_maintenance_step) { maintenance_step_ = static_cast<uint>(vars_->evaluate_value_expression(value, key_maintenance_step)); return true; }
      if (key == key_maintenance_file) { simstate_file_ = vars_->evaluate_string_expression(value, key_maintenance_file); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override {}

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   
   uint parse_seed(const std::string& value) const
   {
      // if parameter is a string-expression, evaluate it and convert to integer (with stroul for hex or decimal format)
      if (*value.begin() == '$')
      {
         const auto eval = vars_->evaluate_string_expression(value, key_seed);
         try { return std::stoul(eval, nullptr, 0); }
         catch (std::runtime_error) { THROW_EXCEPTION(illegal_section_key, "Value for key '" << key_seed << "' not valid. Failed to parse " << eval << " to integer."); }
      }
      
      // or parameter is a value-expression, evaluate it and convert to integer (just floor the value)
      const auto eval = vars_->evaluate_value_expression(value, key_seed);
      return static_cast<uint>(std::floor(eval));
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   uint seed_;
   size_t maintenance_step_;
   std::string simstate_file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const SimulatorSection& ss)
{
   stream << "[" << ss.section_name() << "]" << std::endl;
   if (ss.seed()) stream << ss.key_seed << " = 0x" << std::uppercase << std::hex << ss.seed() << std::dec << std::endl;
   if (ss.prepare_time()) stream << ss.key_prepare_time << " = " << ss.prepare_time() << std::endl;
   if (ss.end_time()) stream << ss.key_end_time << " = " << ss.end_time() << std::endl;
   if (ss.maintenance_step()) stream << ss.key_maintenance_step << " = " << ss.maintenance_step() << std::endl;
   if (ss.maintenance_step() && !ss.maintenance_file().empty()) stream << ss.key_maintenance_file << " = " << ss.maintenance_file() << std::endl;
   stream << std::endl;
   return stream;
}
