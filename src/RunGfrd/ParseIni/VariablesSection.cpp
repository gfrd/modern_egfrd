#define exprtk_disable_enhanced_features
#define exprtk_disable_break_continue
#define exprtk_disable_return_statement
#define exprtk_disable_rtl_io_file
#define exprtk_disable_rtl_vecops
#define exprtk_disable_caseinsensitivity

#include "SectionBase.hpp"
#include "VariablesSection.hpp"
#include "ParserExceptions.hpp"
#include "../../exprtk/exprtk.hpp"           // <-- only include it here (it's huge and takes long to compile)
#include "../../eGFRD/DefsEgfrd.hpp"
#include "convert.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// functions callable from settings-file
namespace functions
{
   inline double concentration(double mol, double volume)
   {
      return convert::particles_in_volume_liter(mol, volume);
   }

   inline double per_nM_per_sec_to_m3_per_sec(double rate)
   {
      return convert::per_nM_per_sec_to_m3_per_sec(rate);
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

VariablesSection::VariablesSection() : SectionBase()
{
   table_ = std::make_unique<exprtk::symbol_table<double>>();
   table_->add_constants();
   table_->add_constant("NA" , N_A);

   value_vars_.reserve(variable_table_size);
   string_vars_.reserve(variable_table_size);

   table_->add_function("concentration", functions::concentration);
   table_->add_function("per_nM_per_sec_to_m3_per_sec", functions::per_nM_per_sec_to_m3_per_sec);
   table_->add_stringvar("dummy", dummy_, false);
}

// --------------------------------------------------------------------------------------------------------------------------------

VariablesSection::~VariablesSection()
{
}

// --------------------------------------------------------------------------------------------------------------------------------

double VariablesSection::evaluate_value_expression(std::string expression, std::string name) const
{
   exprtk::expression<double> expr;
   expr.register_symbol_table(*table_.get());
      
   exprtk::parser<double> parser;
   parser.settings().disable_all_control_structures();

   if (!parser.compile(expression, expr))
      THROW_EXCEPTION(illegal_section_value, "ParseError: " << parser.error() << " in '" << name << " = " << expression << "'\n");

   return expr.value();
}

// --------------------------------------------------------------------------------------------------------------------------------

std::string VariablesSection::evaluate_string_expression(std::string expression, std::string name) const
{
   if (*expression.begin() == '$')
   {
      exprtk::expression<double> expr;
      expr.register_symbol_table(*table_.get());
      
      exprtk::parser<double> parser;
      parser.settings().disable_all_control_structures();

      if (!parser.compile("dummy := " + expression.substr(1), expr))
         THROW_EXCEPTION(illegal_section_value, "ParseError: " << parser.error() << " in '" << name << " = " << expression << "'\n");

      THROW_UNLESS_MSG(illegal_section_value, std::isnan(expr.value()), "Failed to evaluate expression!")
      return dummy_;
   }
      
   // no $ prefix for string expression, return content as result
   return expression;
}

// --------------------------------------------------------------------------------------------------------------------------------

void VariablesSection::add_value_variable(const std::string name, const std::string expression, bool force)
{
   auto i = std::find_if(value_vars_.begin(), value_vars_.end(), [&name](const variable_t<double>& v) { return v.name == name; } );
   if (i != value_vars_.end())
   {
      if ((*i).force && !force) return;  // already set and forced by command line
      THROW_EXCEPTION(already_exists, "Variable already defined: " << name);
   }

   if (value_vars_.size() >= variable_table_size)   // resize of std::vector causes all table reference to be destroyed (no relocation update possible)
      THROW_EXCEPTION(illegal_size, "Variable table is full (max=" << variable_table_size);
      
   double val = evaluate_value_expression(expression, name);
   value_vars_.emplace_back(variable_t<double>{name, expression, val, force});
   table_->add_variable(name, (*value_vars_.rbegin()).value, true);
}

// --------------------------------------------------------------------------------------------------------------------------------

void VariablesSection::add_string_variable(const std::string name, const std::string expression, bool force)
{
   auto i = std::find_if(string_vars_.begin(), string_vars_.end(), [&name](const variable_t<std::string>& v) { return v.name == name; } );
   if (i != string_vars_.end())
   {
      if ((*i).force && !force) return;  // already set and forced by command line
      THROW_EXCEPTION(already_exists, "Variable already defined: " << name);
   }

   if (string_vars_.size() >= variable_table_size)   // resize of std::vector causes all table reference to be destroyed (no relocation update possible)
      THROW_EXCEPTION(illegal_size, "Variable table is full (max=" << variable_table_size);
      
   std::string val = evaluate_string_expression(expression, name);
   string_vars_.emplace_back(variable_t<std::string>{name, expression, val, force});
   table_->add_stringvar(name, (*string_vars_.rbegin()).value, false);
}

// --------------------------------------------------------------------------------------------------------------------------------
