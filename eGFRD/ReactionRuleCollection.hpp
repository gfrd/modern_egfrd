#ifndef REACTIONRULECOLLECTION_HPP
#define REACTIONRULECOLLECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <map>
#include <set>
#include "genericIterator.hpp"
#include "ReactionRule.hpp"
#include "exceptions.hpp"
#include "SerialIDGenerator.hpp"
#include "makeString.hpp"
#include "InteractionRule.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class ReactionRuleCollection
{
public:
   using reaction_rule_set = std::set<ReactionRule>;
   using reaction_rules_map = std::map<ReactionRule::reactants, reaction_rule_set>;
   using reaction_rules_range = gi::iteratorRange<reaction_rules_map>;
   using interaction_rule_set = std::set<InteractionRule>;
   using interaction_rules_map = std::map<InteractionRule::reactants, interaction_rule_set>;
   using interaction_rules_range = gi::iteratorRange<interaction_rules_map>;

   ReactionRuleCollection() noexcept : sgen_(), rr_empty_() {}

   template<typename uniRR, typename = std::enable_if<std::is_same<uniRR, ReactionRule>::value>>
   ReactionRuleID add_reaction_rule(uniRR&& rr)
   {
      auto i(reaction_rules_map_[rr.get_reactants()].insert(std::forward<ReactionRule>(rr)));
      if (!i.second) throw already_exists(make_string() << "ReactionRule already exists: " << (*i.first).get_reactants());
      const_cast<ReactionRule&>(*i.first).set_id(sgen_());
      return (*i.first).id();
   }

   template<typename uniIR, typename = std::enable_if<std::is_same<uniIR, InteractionRule>::value>>
   ReactionRuleID add_interaction_rule(uniIR&& ir)
   {
      auto i(interaction_rules_map_[ir.get_reactants()].insert(std::forward<InteractionRule>(ir)));
      if (!i.second) throw already_exists(make_string() << "InteractionRule already exists: " << (*i.first).get_reactants());
      const_cast<InteractionRule&>(*i.first).set_id(sgen_());
      return (*i.first).id();
   }

   void remove_reaction_rule(const ReactionRule& rr)
   {
      reaction_rules_map_.erase(rr.get_reactants());
   }
   void remove_interaction_rule(const InteractionRule& ir)
   {
      interaction_rules_map_.erase(ir.get_reactants());
   }

   const reaction_rule_set& query_reaction_rules(const SpeciesTypeID s1) const
   {
      auto i(reaction_rules_map_.find(ReactionRule::reactants(s1)));
      if (i == reaction_rules_map_.end()) return rr_empty_;
      return (*i).second;
   };

   const reaction_rule_set& query_reaction_rules(const SpeciesTypeID s1, const SpeciesTypeID s2) const
   {
      auto i(reaction_rules_map_.find(ReactionRule::reactants(s1, s2)));
      if (i == reaction_rules_map_.end()) return rr_empty_;
      return (*i).second;
   }

   const interaction_rule_set& query_interaction_rules(const SpeciesTypeID s1, const StructureTypeID s2) const
   {
      auto i(interaction_rules_map_.find(InteractionRule::reactants(s1, s2)));
      if (i == interaction_rules_map_.end()) return ir_empty_;
      return (*i).second;
   }

   reaction_rules_range get_reaction_rules() const { return reaction_rules_range(reaction_rules_map_); }

   interaction_rules_range get_interaction_rules() const { return interaction_rules_range(interaction_rules_map_); }

   friend class Persistence;

private:
   reaction_rules_map reaction_rules_map_;
   interaction_rules_map interaction_rules_map_;
   SerialIDGenerator<ReactionRuleID> sgen_;
   reaction_rule_set rr_empty_;
   interaction_rule_set ir_empty_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* REACTIONRULECOLLECTION_HPP */
