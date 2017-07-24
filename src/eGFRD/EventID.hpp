#ifndef EVENT_ID_HPP
#define EVENT_ID_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct EventID final : Identifier < EventID >
{
   explicit EventID(const idtype value = idtype(0)) noexcept : Identifier<EventID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < EventID >
   {
      size_t  operator()(const EventID id) const
      {
         return hash<idtype>()(id());
      }
   };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const EventID id)
{
   stream << "EventID(" << id() << ")";
   return stream;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* EVENT_ID_HPP */
