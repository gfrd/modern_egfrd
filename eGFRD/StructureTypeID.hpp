#ifndef STRUCTURE_TYPE_ID_HPP
#define STRUCTURE_TYPE_ID_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct StructureTypeID final : Identifier < StructureTypeID >
{
   explicit StructureTypeID(const idtype value = idtype(0)) noexcept : Identifier<StructureTypeID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < StructureTypeID >
   {
      size_t  operator()(const StructureTypeID id) const
      {
         return hash<idtype>()(id());
      }
   };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const StructureTypeID id)
{
   stream << "StructureTypeID(" << id() << ")";
   return stream;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* STRUCTURE_TYPE_ID_HPP */
