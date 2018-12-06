#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <cstring>
#include <string>
#include <vector>
#include <regex>

// --------------------------------------------------------------------------------------------------------------------------------

class getoptions
{
public:
   explicit getoptions(int argc, char**argv)
   {
      parse(argc, argv);
   }

private:

   using OptionsVector = std::vector<std::pair<bool, std::string>>;

   void parse(int argc, char**argv)
   {
      options_.clear();
      for (int i = 1; i < argc; i++)      // skip executable path
      {
         // when start with - (but not a negative number)
         bool param = (argv[i][0] == '-') && !(argv[i][1] == 0 || argv[i][1] == '.' || (argv[i][1] >= '0' && argv[i][1] <= '1'));
         std::string option = param ? argv[i] + 1 : argv[i];

         if (option[0] == '\'' || option[0] == '"')
         {
            char quote = option[0];
            size_t len = option.length();

            if (option[len - 1] == quote)
               option = option.substr(1, len - 2);
            else
            {
               option = option.substr(1, len - 1);
               while (++i < argc)
               {
                  option = option + ' ' + argv[i];
                  if (argv[i][std::strlen(argv[i]) - 1] == quote) break;
               }

               len = option.length();
               if (option[len - 1] == quote)
                  option.resize(len - 1);
            }
         }

         //std::cout << "Handle cmd-arg: " << i << " = " << argv[i] << "  => Param: " << param << " Option: " << option << "\n";

         options_.emplace_back(std::make_pair(param, option));
      }
   }

public:
   size_t size() const { return options_.size(); }
   bool isparam(size_t i) const { return i >= 0 && i < size() ? options_[i].first : false; }
   bool isvalue(size_t i) const { return i >= 0 && i < size() ? !options_[i].first : false; }

   bool isvalue_NH(size_t i) const { return (i >= 0 && i < size() ? !options_[i].first : false) && ismatch(i, "(0[xX])?[0-9a-fA-F]+"); }     // may be hex
   bool isvalue_NP(size_t i) const { return (i >= 0 && i < size() ? !options_[i].first : false) && ismatch(i, "[+]?[\\d]+"); }     // positive
   bool isvalue_N(size_t i) const { return (i >= 0 && i < size() ? !options_[i].first : false) && ismatch(i, "[+\\-]?[\\d]+"); }    // integer number
   bool isvalue_F(size_t i) const { return (i >= 0 && i < size() ? !options_[i].first : false) && ismatch(i, "[^<>\"|?*%]+"); }         // file/path (no < > " | ? * %)
   bool isvalue_D(size_t i) const { return (i >= 0 && i < size() ? !options_[i].first : false) && ismatch(i, "[+\\-]?\\d*(?:\\.\\d*)?(?:[eE][+\\-]?\\d+)?"); }    // scientific number

   std::string option(size_t i) const { return i >= 0 && i < size() ? options_[i].second : ""; }
   double value(size_t i) const { return i >= 0 && i < size() ? std::stod(options_[i].second) : 0.0; }
   int number(size_t i) const { return i >= 0 && i < size() ? std::stoi(options_[i].second) : 0; }

private:

   bool ismatch(size_t i, std::string pattern) const
   {
      const std::regex regex(pattern);
      std::smatch match;
      return std::regex_match(options_[i].second, match, regex);
   }

   OptionsVector options_;
};

// --------------------------------------------------------------------------------------------------------------------------------
