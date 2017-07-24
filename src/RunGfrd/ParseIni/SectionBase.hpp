#ifndef SECTIONBASE_HPP
#define SECTIONBASE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <regex>
#include <string>

// --------------------------------------------------------------------------------------------------------------------------------

struct SectionBase
{
   explicit SectionBase() {}
   virtual ~SectionBase() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void set_keypair(const std::string& key, const std::string& value) = 0;

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   // --------------------------------------------------------------------------------------------------------------------------------

   static bool get_bool(const std::string& value)
   {
      std::string lcvalue;
      std::transform(value.begin(), value.end(), std::back_inserter(lcvalue), std::bind(std::tolower<char>, std::placeholders::_1, std::locale()));
      if (lcvalue == "true" || lcvalue == "yes" || lcvalue == "1") return true;
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------
};

#endif