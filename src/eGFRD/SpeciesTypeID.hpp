#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct SpeciesTypeID final : Identifier < SpeciesTypeID >
{
   explicit SpeciesTypeID(const idtype value = idtype(0)) noexcept : Identifier<SpeciesTypeID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < SpeciesTypeID >
   {
      size_t  operator()(const SpeciesTypeID id) const
      {
         return hash<idtype>()(id());
      }
   };
} // namespace std

inline std::ostream& operator<<(std::ostream& stream, const SpeciesTypeID id)
{
   stream << "SpeciesTypeID(" << id() << ")";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------
