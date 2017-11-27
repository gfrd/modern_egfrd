#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"

   // --------------------------------------------------------------------------------------------------------------------------------

struct WorldSection final : SectionBase
{
   explicit WorldSection() : SectionBase(), matrix_space_(MATRIXSIZE) { init_auto_vars({{key_world_size, 0.0}}); }

   ~WorldSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "World"; }

   const std::string key_matrix_space = "Matrix";
   uint matrix_space() const { return matrix_space_; }

   const std::string key_world_size = "Size";
   double world_size() const { return auto_var_value(key_world_size); }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionBase::set_keypair(key, value)) return true;
      if (key == key_matrix_space) { matrix_space_ = static_cast<uint>(vars_->evaluate_value_expression(value, key)); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override {}

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   uint matrix_space_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const WorldSection& ws)
{
   stream << "[" << ws.section_name() << "]" << std::endl;
   if (ws.matrix_space()) stream << ws.key_matrix_space << " = " << ws.matrix_space() << std::endl;
   stream << ws.key_world_size << " = " << ws.world_size() << std::endl;
   stream << std::endl;
   return stream;
}
