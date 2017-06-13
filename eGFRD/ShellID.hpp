#ifndef SHELL_ID_HPP
#define SHELL_ID_HPP

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
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const ShellID id)
{
   stream << "ShellID(" << id() << ")";
   return stream;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* SHELL_ID_HPP */
