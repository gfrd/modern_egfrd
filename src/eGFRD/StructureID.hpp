#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct StructureID final : Identifier < StructureID >
{
   explicit StructureID(const idtype value = idtype(0)) noexcept : Identifier<StructureID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < StructureID >
   {
      size_t  operator()(const StructureID id) const
      {
         return hash<idtype>()(id());
      }
   };
} // namespace std

inline std::ostream& operator<<(std::ostream& stream, const StructureID id)
{
   stream << "StructureID(" << id() << ")";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------
