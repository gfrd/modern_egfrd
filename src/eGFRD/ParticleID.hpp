#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ParticleID final : Identifier < ParticleID >
{
   explicit ParticleID(const idtype value = idtype(0)) noexcept : Identifier<ParticleID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < ParticleID >
   {
      size_t  operator()(const ParticleID id) const
      {
         return hash<idtype>()(id());
      }
   };
} // namespace std

inline std::ostream& operator<<(std::ostream& stream, const ParticleID id)
{
   stream << "PID(" << id() << ")";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------
