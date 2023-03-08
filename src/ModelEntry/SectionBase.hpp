#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <regex>
#include <string>
#include <functional>
#include <unordered_map>
#include <Vector3.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

#if defined(_MSC_VER)
#if defined(ModelEntry_EXPORTS)
#define ME_EXPORT __declspec(dllexport)
#else
#define ME_EXPORT __declspec(dllimport)
#endif
#else
#define ME_EXPORT
#endif

// --------------------------------------------------------------------------------------------------------------------------------

struct VariablesSection;

// --------------------------------------------------------------------------------------------------------------------------------

struct ME_EXPORT SectionBase
{
   explicit SectionBase() : vars_(nullptr) { }
   virtual ~SectionBase() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual bool set_keypair(const std::string& key, const std::string& value);

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_vars(VariablesSection *vars) { vars_ = vars; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void init_auto_vars(std::initializer_list<std::pair<std::string, double>> keys)
   {
      for (auto& kvp : keys)
         auto_values_[kvp.first] = kvp.second;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double auto_var_value(const std::string& name) const { return auto_values_.at(name); }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void PrintSettings() const = 0;

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   bool is_valid_speciestype_name(const std::string& name) const
   {
      const std::regex regex(regex_speciestype_name_);
      std::smatch match;
      return std::regex_match(name, match, regex);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool is_valid_variable_name(const std::string& name) const
   {
      const std::regex regex(regex_variable_name_);
      std::smatch match;
      return std::regex_match(name, match, regex);
   }

    // --------------------------------------------------------------------------------------------------------------------------------

    static bool get_bool(const std::string& value)
    {
        std::string lcvalue;
        std::transform(value.begin(), value.end(), std::back_inserter(lcvalue), std::bind(std::tolower<char>, std::placeholders::_1, std::locale()));
        return lcvalue == "true" || lcvalue == "yes" || lcvalue == "1";
    }

    // --------------------------------------------------------------------------------------------------------------------------------

   VariablesSection* vars_;

private:

   std::unordered_map<std::string, double> auto_values_;

   const std::string regex_speciestype_name_ = "[a-zA-Z][\\w\\-*_']*";
   const std::string regex_variable_name_ = "[a-zA-Z][a-zA-Z0-9_]*";

   // --------------------------------------------------------------------------------------------------------------------------------
};


struct ME_EXPORT SectionModeBase : SectionBase
{
   explicit SectionModeBase() : mode_(modes::Off) { }

   virtual ~SectionModeBase() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   enum class modes { Off = 0, On = 1, Run = 2, };
   const std::string key_mode = "Mode";

   // --------------------------------------------------------------------------------------------------------------------------------

   modes mode() const { return mode_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (SectionBase::set_keypair(key, value)) return true;
      if (key == key_mode) { mode_ = get_mode(value); return true; }
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   static modes get_mode(const std::string& value)
   {
      std::string lcvalue;
      std::transform(value.begin(), value.end(), std::back_inserter(lcvalue), std::bind(std::tolower<char>, std::placeholders::_1, std::locale()));
      if (lcvalue == "off" || lcvalue == "0") return modes::Off;
      if (lcvalue == "on" || lcvalue == "yes" || lcvalue == "1") return modes::On;
      if (lcvalue == "run" || lcvalue == "2") return modes::Run;
      return modes::Off;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:
   modes mode_;

   // --------------------------------------------------------------------------------------------------------------------------------
};
