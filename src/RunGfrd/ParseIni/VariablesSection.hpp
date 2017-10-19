#ifndef VARIABLES_SECTION_HPP
#define VARIABLES_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "../../eGFRD/DefsEgfrd.hpp"
#include <iomanip>
#include "ParserExceptions.hpp"
#include <iostream>

// --------------------------------------------------------------------------------------------------------------------------------

// forward declare sybol_table type
namespace exprtk
{
   template<typename T>
   class symbol_table;
}

// --------------------------------------------------------------------------------------------------------------------------------

struct VariablesSection final : SectionBase
{
   // --------------------------------------------------------------------------------------------------------------------------------

   const int variable_table_size = 100;
   
   // --------------------------------------------------------------------------------------------------------------------------------

   template<typename T>
   struct variable_t
   {
      std::string name;
      std::string expr;
      T value;
      bool force;
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   explicit VariablesSection();
   ~VariablesSection();

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "Variables"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (!is_valid_variable_name(key)) THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' is not a valid variable name.");
      
      add_variable(key, value, false);
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   const std::vector<variable_t<double>>& value_vars() const { return value_vars_; }
   const std::vector<variable_t<std::string>>& string_vars() const { return string_vars_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   double evaluate_value_expression(std::string expression, std::string name) const;

   std::string evaluate_string_expression(std::string expression, std::string name) const;

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_variable(const std::string name, const std::string expression, bool force)
   {
      if (*expression.begin() == '$')
         add_string_variable(name, expression, force);
      else
         add_value_variable(name, expression, force);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_value_variable(const std::string name, const std::string expression, bool force);

   void add_string_variable(const std::string name, const std::string expression, bool force);

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() const override
   {
      for (auto& var : value_vars_)
         std::cout << std::setw(11) << var.name << " = " << var.expr << " => " << var.value << "\n";
      for (auto& var : string_vars_)
         std::cout << std::setw(11) << var.name << " = " << var.expr << " => " << var.value << "\n";
      if (!value_vars_.empty() || !string_vars_.empty()) std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   std::vector<variable_t<double>> value_vars_;
   std::vector<variable_t<std::string>> string_vars_;
   std::unique_ptr<exprtk::symbol_table<double>> table_;
   std::string dummy_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const VariablesSection& vs)
{
   stream << "[" << vs.section_name() << "]" << std::endl;

      for (auto& var : vs.value_vars())
         stream << var.name << " = " << var.expr << std::endl;
      for (auto& var : vs.string_vars())
         stream << var.name << " = " << var.expr << std::endl;
   stream << std::endl;
   return stream;
}

#endif