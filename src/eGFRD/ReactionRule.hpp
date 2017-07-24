#ifndef REACTION_RULE_HPP
#define REACTION_RULE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <vector>
#include <ostream>
#include <algorithm>
#include "SpeciesTypeID.hpp"
#include "ReactionRuleID.hpp"
#include "Reactants.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class ReactionRule
{
public:
   using reactants = Reactants<SpeciesTypeID, SpeciesTypeID>;
   using products = std::vector<SpeciesTypeID>;

   explicit ReactionRule(const SpeciesTypeID s1, double k, const products& products) : ReactionRule(reactants(s1), k, products) {}

   explicit ReactionRule(const SpeciesTypeID s1, const SpeciesTypeID s2, double k, const products& products) : ReactionRule(reactants(s1, s2), k, products) {}

   explicit ReactionRule(const reactants& reactants, double k, const products& products) : id_(), reactants_(reactants), products_(products), k_(k)
   {
      THROW_UNLESS_MSG(illegal_argument, reactants_.size() > 0, "ReactionRule need at least one SpeciesTypeID");
      std::sort(products_.begin(), products_.end());        // sort for easier duplicate checks
   }

   ReactionRuleID id() const { return id_; }

   double getK() const { return k_; }

   const reactants& get_reactants() const { return reactants_; }

   const products& get_products() const { return products_; }

   bool operator<(const ReactionRule& rhs) const
   {
      // compare reactants and products (ignore id and k) will detect duplicates in set
      if (reactants_ < rhs.get_reactants()) return true;
      return reactants_ == rhs.get_reactants() ? products_ < rhs.get_products() : false;
   }

   bool operator==(const ReactionRule& rhs) const
   {
      return reactants_ == rhs.get_reactants() && products_ == rhs.get_products() && k_ == rhs.k_;
   }

   bool operator!=(const ReactionRule& rhs) const { return !(*this == rhs); }

   static ReactionRule empty() { return ReactionRule(SpeciesTypeID(1), 0.0, std::vector<SpeciesTypeID>()); }

protected:
   friend class ReactionRuleCollection;
   void set_id(const ReactionRuleID id) { ASSERT(id); id_ = id; }

   friend class Persistence;

private:
   ReactionRuleID id_;
   reactants reactants_;
   products products_;
   double k_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ReactionRule& rr)
{
   stream << "ReactionRule{" << rr.id() << ", reactants=[" << rr.get_reactants();
   stream << "], k=" << rr.getK() << ", products=[";
   bool first = true;
   for (auto& p : rr.get_products())
   {
      if (!first) stream << ", ";
      stream << p;
      first = false;
   }
   stream << "]}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* REACTION_RULE_HPP */
