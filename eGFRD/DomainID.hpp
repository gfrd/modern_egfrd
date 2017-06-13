#ifndef DOMAIN_ID_HPP
#define DOMAIN_ID_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct DomainID final : Identifier < DomainID >
{
   explicit DomainID(const idtype value = idtype(0)) noexcept : Identifier<DomainID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < DomainID >
   {
      size_t  operator()(const DomainID id) const
      {
         return hash<idtype>()(id());
      }
   };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, DomainID id)
{
   stream << "DID(" << id() << ")";
   return stream;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* DOMAIN_ID_HPP */
