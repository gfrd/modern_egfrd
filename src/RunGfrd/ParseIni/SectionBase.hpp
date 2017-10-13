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

   virtual bool set_keypair(const std::string& key, const std::string& value) = 0;

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


struct SectionModeBase : SectionBase
{
   explicit SectionModeBase() : mode_(modes::Off)  {}
   virtual ~SectionModeBase() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   enum class modes { Off = 0, On = 1, Run = 2, };
   const std::string key_mode = "Mode";

   // --------------------------------------------------------------------------------------------------------------------------------

   modes mode() const { return mode_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
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




#endif