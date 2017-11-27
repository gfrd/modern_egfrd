#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ShellID final : Identifier < ShellID >
{
   explicit ShellID(const idtype value = idtype(0)) noexcept : Identifier<ShellID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < ShellID >
   {
      size_t  operator()(const ShellID id) const
      {
         return hash<idtype>()(id());
      }
   };
}

inline std::ostream& operator<<(std::ostream& stream, const ShellID id)
{
   stream << "ShellID(" << id() << ")";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------
