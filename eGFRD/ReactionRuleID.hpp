#ifndef NETWORK_RULE_ID_HPP
#define NETWORK_RULE_ID_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "DefsEgfrd.hpp"
#include "Identifier.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ReactionRuleID final : Identifier < ReactionRuleID >
{
   explicit ReactionRuleID(const idtype value = idtype(0)) noexcept : Identifier<ReactionRuleID>(value) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

namespace std
{
   template<>
   struct hash < ReactionRuleID >
   {
      size_t  operator()(const ReactionRuleID id) const
      {
         return hash<idtype>()(id());
      }
   };
}; // namespace std

inline std::ostream& operator<<(std::ostream& stream, const ReactionRuleID id)
{
   stream << "RRID(" << id() << ")";
   return stream;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* NETWORK_RULE_ID_HPP */
