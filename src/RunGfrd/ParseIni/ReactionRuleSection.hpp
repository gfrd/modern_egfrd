#ifndef REACTIONRULE_SECTION_HPP
#define REACTIONRULE_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <array>
#include "ReactionRule.hpp"
#include "ParserExceptions.hpp"
#include "Model.hpp"
#include "SectionBase.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ReactionRuleSection final : SectionModeBase
{
   explicit ReactionRuleSection() : rule_(std::string()), k_(-1), ka_(-1), kd_(-1), is_bidirectional_(false) { mode_ = modes::On; }
   ~ReactionRuleSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "ReactionRule"; }

   const std::string key_rule = "Rule";
   const std::string& rule() const { return rule_; }
   const std::string key_k = "k";
   const std::string key_ka = "ka";
   const std::string key_kd = "kd";
   double k() const { return k_; }
   double ka() const { return ka_; }
   double kd() const { return kd_; }
   bool is_bidirectional() const { return is_bidirectional_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key,value)) return true;
      if (key == key_rule) { rule_ = format_check(value); return true; }
      if (key == key_k) { k_ = std::stod(value); return true; }
      if (key == key_ka) { ka_ = std::stod(value); return true; }
      if (key == key_kd) { kd_ = std::stod(value); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void create_reaction_rule(const Model& model, ReactionRuleCollection& rules) const
   {
      auto reactants = reactants_.size();
      auto products = products_.size();
      
      THROW_UNLESS_MSG(illegal_section_value, reactants > 0 && reactants <= 2, "Rule should have one or two reactant particles.");
      THROW_UNLESS_MSG(illegal_section_value, products >= 0 && products <= 2, "Rule should have no more than three product particles.");
      
      // find species ID's for given names (may throw exception when not found)
      auto reactant = std::vector<SpeciesTypeID>();
      for (int i=0; i<reactants; ++i) reactant.emplace_back(model.get_species_type_id_by_name(reactants_[i]));
      auto product = std::vector<SpeciesTypeID>();
      for (int i=0; i<products; ++i) product.emplace_back(model.get_species_type_id_by_name(products_[i]));

      if (is_bidirectional_)
      {
         THROW_UNLESS_MSG(illegal_section_value, ka_ >= 0 && kd_ >= 0, "Bidirectional rule should specify ka and kd values.");
         THROW_UNLESS_MSG(illegal_section_value, products == 1  , "Bidirectional rule should have one product.");

         // associate
         switch (reactants)
         {
         case 1: rules.add_reaction_rule( ReactionRule(reactant[0], ka_, product)); break;
         case 2: rules.add_reaction_rule( ReactionRule(reactant[0], reactant[1], ka_, product)); break;
         default: break;
         }

         // dissociate
         rules.add_reaction_rule( ReactionRule(product[0], kd_, reactant));
      }
      else
      {
         THROW_UNLESS_MSG(illegal_section_value, k_ >= 0, "Rule should specify a value for k.");
         switch (reactants)
         {
         case 1: rules.add_reaction_rule( ReactionRule(reactant[0], k_, product)); break;
         case 2: rules.add_reaction_rule( ReactionRule(reactant[0], reactant[1], k_, product)); break;
         default: break;
         }
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   std::string format_check(std::string input)
   {
      auto regex = std::regex("[^\\s]+");
      auto begin =  std::sregex_iterator(input.begin(), input.end(), regex);
      auto end = std::sregex_iterator();
      bool has_arrow = false;
 
      std::stringstream ss;
      for (auto i = begin; i != end; ++i) 
      {
         std::string name = i->str();
         if (name == "+") { ss << " + "; continue; }
         if (name.find("<") != std::string::npos) { is_bidirectional_ = has_arrow = true; ss << " <=> "; continue;}
         if (name.find(">") != std::string::npos) { has_arrow = true; ss << " -> "; continue;}
         
         std::smatch match;
         auto regex2 = std::regex("[a-zA-Z][\\w\\-*_']*");
         if (!std::regex_match(name, match, regex2)) THROW_EXCEPTION(illegal_section_value, "SpeciesTypeName '" << name << "'not valid");
         if (has_arrow) products_.emplace_back(name); else reactants_.emplace_back(name);
         ss << name;
      }
      return ss.str();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string rule_;
   double k_,ka_,kd_;
   bool is_bidirectional_;
   std::vector<std::string> reactants_;
   std::vector<std::string> products_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ReactionRuleSection& rrs)
{
   stream << "[" << rrs.section_name() << "]" << std::endl;
   stream << rrs.key_rule << " = " << rrs.rule() << std::endl;
   if (rrs.is_bidirectional())
   {
      stream << rrs.key_ka << " = " << rrs.ka() << std::endl;
      stream << rrs.key_kd << " = " << rrs.kd() << std::endl;
   }
   else
      stream << rrs.key_k << " = " << rrs.k() << std::endl;
   stream << std::endl;
   return stream;
}

#endif
