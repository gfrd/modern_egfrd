#ifndef REACTIONRULE_SECTION_HPP
#define REACTIONRULE_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <array>
#include "ReactionRule.hpp"
#include "ParserExceptions.hpp"
#include "Model.hpp"
#include "SectionBase.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ReactionRuleSection final : SectionBase
{
   explicit ReactionRuleSection() : rule_(std::string()), k_(0) { }
   ~ReactionRuleSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "ReactionRule"; }

   const std::string key_rule = "Rule";
   const std::string& rule() const { return rule_; }
   const std::string key_k = "k";
   double k() const { return k_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_keypair(const std::string& key, const std::string& value) override
   {
      if (key == key_rule) { rule_ = value; return; }
      if (key == key_k) { k_ = std::stod(value); return; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   ReactionRule create_reaction_rule(const Model& model) const
   {
      auto pattern = std::regex("(\\w|-?->)");
      std::smatch match;

      auto reactant = std::vector<SpeciesTypeID>();
      auto product = std::vector<SpeciesTypeID>();
      bool has_arrow = false;

      std::string line = rule_;
      while (std::regex_search(line, match, pattern))
      {
         if (match[0].str().find("->") != std::string::npos) has_arrow = true;
         else
         {
            auto sid = model.get_species_type_id_by_name(match[0]);
            (has_arrow ? product : reactant).emplace_back(sid);
         }
         line = match.suffix();
      }

      THROW_UNLESS_MSG(illegal_section_value, has_arrow, "ReactionRule '" << rule_ << "' is not valid.");
      THROW_UNLESS_MSG(illegal_section_value, product.size() < 3, "ReactionRule '" << rule_ << "' has to many products.");

      switch (reactant.size())
      {
      case 1:
         return ReactionRule(reactant[0], k_, product);
      case 2:
         return ReactionRule(reactant[0], reactant[1], k_, product);
      default:
         THROW_EXCEPTION(illegal_section_value, "ReactionRule '" << rule_ << "' has to invalid number of reactants.");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   std::string rule_;
   double k_;
};


// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ReactionRuleSection& rrs)
{
   stream << "[" << rrs.section_name() << "]" << std::endl;
   stream << rrs.key_rule << " = " << rrs.rule() << std::endl;
   stream << rrs.key_k << " = " << rrs.k() << std::endl;
   stream << std::endl;
   return stream;
}

#endif
