#include <string>
#include <iostream>
#include "ParserExceptions.hpp"
#include "SectionBase.hpp"
#include "VariablesSection.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

bool SectionBase::set_keypair(const std::string& key, const std::string& value)
{
   if (auto_values_.find(key) == auto_values_.end()) return false;
   auto_values_[key] = vars_->evaluate_value_expression(value, key);
   return true;
}

// --------------------------------------------------------------------------------------------------------------------------------
