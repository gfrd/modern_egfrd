#ifndef WORLDSECTION_HPP
#define WORLDSECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"

   // --------------------------------------------------------------------------------------------------------------------------------

struct WorldSection final : SectionBase
{
   explicit WorldSection() : matrix_space_(0), world_size_(0) { }
   ~WorldSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "World"; }

   const std::string key_matrix_space = "Matrix";
   uint matrix_space() const { return matrix_space_; }

   const std::string key_world_size = "Size";
   double world_size() const { return world_size_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_keypair(const std::string& key, const std::string& value) override
   {
      if (key == key_matrix_space) { matrix_space_ = std::stoul(value); return; }
      if (key == key_world_size) { world_size_ = std::stod(value); return; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   uint matrix_space_;
   double world_size_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const WorldSection& ws)
{
   stream << "[" << ws.section_name() << "]" << std::endl;
   if (ws.matrix_space()) stream << ws.key_matrix_space << " = " << ws.matrix_space() << std::endl;
   if (ws.world_size()) stream << ws.key_world_size << " = " << ws.world_size() << std::endl;
   stream << std::endl;
   return stream;
}

#endif
