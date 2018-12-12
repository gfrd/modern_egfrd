#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "ReactionRule.hpp"
#include "ParserExceptions.hpp"
#include "Model.hpp"
#include "ReactionRuleCollection.hpp"
#include "SectionBase.hpp"
#include <iomanip>

// --------------------------------------------------------------------------------------------------------------------------------

struct ME_EXPORT ReactionRuleSection final : SectionModeBase
{
   explicit ReactionRuleSection() : SectionModeBase(), rule_(std::string()), is_bidirectional_(false), is_bimolecular_(false)
   {
      mode_ = modes::On;
      init_auto_vars( { { key_k, -1}, { key_k1, -1}, { key_k2, -1}, { key_ka, -1}, { key_kd, -1} } );
   }

   ~ReactionRuleSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "ReactionRule"; }

   const std::string key_rule = "Rule";
   const std::string& rule() const { return rule_; }
   const std::string key_k = "k";
   const std::string key_k1 = "k1";
   const std::string key_k2 = "k2";
   const std::string key_ka = "ka";
   const std::string key_kd = "kd";
   double k() const { return auto_var_value(key_k); }
   double k1() const { return auto_var_value(key_k1); }
   double k2() const { return auto_var_value(key_k2); }
   double ka() const { return auto_var_value(key_ka); }
   double kd() const { return auto_var_value(key_kd); }
   bool is_bidirectional() const { return is_bidirectional_; }
   bool is_bimolecular() const { return is_bimolecular_; }


   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionModeBase::set_keypair(key,value)) return true;
      if (key == key_rule) { rule_ = format_check(value); return true; }
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
         double rk1 = is_bimolecular() ? ka() : k1();
         double rk2 = is_bimolecular() ? kd() : k2();

         THROW_UNLESS_MSG(illegal_section_value, rk1 >= 0 && rk2 >= 0, "Bidirectional " << (is_bimolecular() ? "bi-" : "mono-") << "molecular rule should specify " << (is_bimolecular() ? "ka and kd" : "k1 and k2") << " values.");
         THROW_UNLESS_MSG(illegal_section_value, products == 1  , "Bidirectional rule should have one product.");

         // associate
         switch (reactants)
         {
         case 1: rules.add_reaction_rule( ReactionRule(reactant[0], rk1, product)); break;
         case 2: rules.add_reaction_rule( ReactionRule(reactant[0], reactant[1], rk1, product)); break;
         default: break;
         }

         // dissociate
         rules.add_reaction_rule( ReactionRule(product[0], rk2, reactant));
      }
      else
      {
         double rk = k();
         THROW_UNLESS_MSG(illegal_section_value, rk >= 0, "Rule should specify a value for k.");
         switch (reactants)
         {
         case 1: rules.add_reaction_rule( ReactionRule(reactant[0], rk, product)); break;
         case 2: rules.add_reaction_rule( ReactionRule(reactant[0], reactant[1], rk, product)); break;
         default: break;
         }
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      std::string mode;
      if (mode_ != modes::On) mode = mode_ == modes::Off ? ", Mode = Off" : ", Mode = Run";
      const std::string unit = is_bimolecular_ ? " [m^3*s^-1]" : " [s^-1]";

      if (!is_bidirectional_)
         std::cout << std::setw(14) << "rule = " << "'" << rule() << "'" << ", k = " << k() << unit << mode << "\n";
      else
         if (!is_bimolecular_)
            std::cout << std::setw(14) << "rule = " << "'" << rule() << "'" << ", k1 = " << k1() << " [s^-1], k2 = " << k2() << " [s^-1]" << mode << "\n";
         else 
            std::cout << std::setw(14) << "rule = " << "'" << rule() << "'" << ", ka = " << ka() << " [m^3*s^-1], kd = " << kd() << " [s^-1]" << mode << "\n";
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
         if (name == "+") { ss << " + "; is_bimolecular_ = true; continue; }
         if (name.find("<") != std::string::npos) { is_bidirectional_ = has_arrow = true; ss << " <=> "; continue;}
         if (name.find(">") != std::string::npos) { has_arrow = true; ss << " -> "; continue;}
         
         if (!is_valid_speciestype_name(name)) THROW_EXCEPTION(illegal_section_value, "SpeciesTypeName '" << name << "'not valid");
         if (has_arrow) products_.emplace_back(name); else reactants_.emplace_back(name);
         ss << name;
      }
      return ss.str();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string rule_;
   bool is_bidirectional_, is_bimolecular_;
   std::vector<std::string> reactants_;
   std::vector<std::string> products_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const ReactionRuleSection& rrs)
{
   stream << "[" << rrs.section_name() << "]" << std::endl;
   if (rrs.mode() != SectionModeBase::modes::On) stream << rrs.key_mode << " = " << (rrs.mode() == SectionModeBase::modes::Run ? "Run" : rrs.mode() == SectionModeBase::modes::On ? "On" : "Off") << std::endl;
   stream << rrs.key_rule << " = " << rrs.rule() << std::endl;
   if (rrs.is_bidirectional())
   {
      stream << rrs.key_k1 << " = " << rrs.k1() << std::endl;
      stream << rrs.key_k2 << " = " << rrs.k2() << std::endl;
   }
   else
      stream << rrs.key_k << " = " << rrs.k() << std::endl;
   stream << std::endl;
   return stream;
}
