#ifndef INTERACTION_RULE_HPP
#define INTERACTION_RULE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <vector>
#include <ostream>
#include <algorithm>
#include "SpeciesTypeID.hpp"
#include "ReactionRuleID.hpp"
#include "Reactants.hpp"
#include "StructureTypeID.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class InteractionRule
{
public:
   using reactants = Reactants<SpeciesTypeID, StructureTypeID>;
   using products = std::vector<SpeciesTypeID>;

   explicit InteractionRule(const SpeciesTypeID species, const StructureTypeID structure, double k, const products& products) : InteractionRule(reactants(species, structure), k, products) {}

   explicit InteractionRule(const reactants& reactants, double k, const products& products) : id_(), reactants_(reactants), products_(products), k_(k)
   {
      THROW_UNLESS_MSG(illegal_argument, reactants_.size() == 2,"InteractionRule needs SpeciesTypeID and StructrureTypeID");
      std::sort(products_.begin(), products_.end());        // sort for easier duplicate checks
   }

   ReactionRuleID id() const { return id_; }

   double getK() const { return k_; }

   const reactants& get_reactants() const { return reactants_; }

   const products& get_products() const { return products_; }

   bool operator<(const InteractionRule& rhs) const
   {
      if (reactants_ < rhs.get_reactants()) return true;
      return reactants_ == rhs.get_reactants() ? products_ < rhs.get_products() : false;
   }

   bool operator==(const InteractionRule& rhs) const
   {
      return reactants_ == rhs.get_reactants() && products_ == rhs.get_products() && k_ == rhs.k_;
   }

   bool operator!=(const InteractionRule& rhs) const { return !(*this == rhs); }


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

inline std::ostream& operator<<(std::ostream& stream, const InteractionRule& rr)
{
   stream << "InteractionRule{" << rr.id() << ", reactants=[" << rr.get_reactants();
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

#endif /* INTERACTION_RULE_HPP */
