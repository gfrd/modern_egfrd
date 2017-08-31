// --------------------------------------------------------------------------------------------------------------------------------

#include <cassert>
#include <unordered_map>
#include <memory>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <locale>
#include <functional>
#include "Logger.hpp"
#include "gfrd_compat.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

#if defined(_DEBUG) && defined(_MSC_VER)
#include <windows.h>
#endif

// --------------------------------------------------------------------------------------------------------------------------------

Logger& Logger::get_logger(const char* name)
{
   static std::unordered_map<std::string, std::unique_ptr<Logger>> logger_table;

   std::string lcname;
   std::transform(name, name + std::strlen(name), std::back_inserter(lcname), std::bind(std::tolower<char>, std::placeholders::_1, std::locale()));

   auto i = logger_table.insert(std::make_pair(lcname, nullptr));
   if (i.second) (*i.first).second = std::make_unique<Logger>(name);
   return *(*i.first).second;
}

// --------------------------------------------------------------------------------------------------------------------------------

const char* Logger::stringize_error_level(loglevel lv)
{
   static const char* names[] = { "OFF", "DEBUG", "INFO", "WARN", "ERROR", "FATAL" };
   assert(static_cast<size_t>(lv) < sizeof(names) / sizeof(*names));
   return names[static_cast<size_t>(lv)];
}

// --------------------------------------------------------------------------------------------------------------------------------

void Logger::header(loglevel lv, std::ostream& stream) const
{
   if (prefix_) stream << prefix_;

   if ((flags_ & logflags::Clock) != logflags::None)
   {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      auto start = std::chrono::system_clock::from_time_t(time); // convert back floors to whole seconds
      auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();  // get milli seconds
      tm* tm = std::localtime(&time);
      stream << std::setfill('0') << std::setw(2) << tm->tm_hour << ":" << std::setw(2) << tm->tm_min << ":" << std::setw(2) << tm->tm_sec << "." << std::setw(3) << msec << std::setfill(' ') << " ";
   }

   if ((flags_ & logflags::Name) != logflags::None) stream << name_ << " ";
   if ((flags_ & logflags::Level) != logflags::None) stream << "[" << stringize_error_level(lv) << "] ";
   if ((flags_ & logflags::Separator) != logflags::None) stream << ": ";
}

// --------------------------------------------------------------------------------------------------------------------------------

#if defined(_MSC_VER) && defined(_DEBUG)
LOGGER_API void Logger::iostub::ToDebugConsole(const std::stringstream& s) const
{
   if ((log_.flags_ & logflags::Debugger) != logflags::None) OutputDebugStringA(s.str().c_str());
}
#endif

// --------------------------------------------------------------------------------------------------------------------------------
