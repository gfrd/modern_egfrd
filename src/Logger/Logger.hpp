#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#if defined(_MSC_VER)
#if defined(Logger_EXPORTS)
#define LOGGER_API __declspec(dllexport)
#else
#define LOGGER_API __declspec(dllimport)
#endif
#else
#define LOGGER_API
#endif

// --------------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include "bitmask_operators.hpp"

#if defined(_MSC_VER) && defined(_DEBUG)
#include <sstream>
#endif

// --------------------------------------------------------------------------------------------------------------------------------

class Logger final
{
public:

   // levels
   enum class loglevel { Off = 0, Debug = 1, Info = 2, Warning = 3, Error = 4, Fatal = 5 };
   enum class logflags { None = 0, Clock = 1, Name = 2, Level = 4, Separator = 8, Debugger = 16, All = 0xffff };

   // static methods
   LOGGER_API static Logger& get_logger(const char* name);
   LOGGER_API static const char* stringize_error_level(loglevel lv);

   // constructors (no default, need name)
   Logger() = delete;
   explicit Logger(std::string name, loglevel level = loglevel::Info, std::ostream& stream = std::clog) noexcept :
      name_(name), level_(level), stream_(&stream), prefix_(nullptr), postfix_("\n"), flags_(logflags::All) {}

   // setters
   void set_level(loglevel level) { level_ = level; }
   void set_stream(std::ostream& stream) { stream_ = &stream; }
   void set_prefix(const char* prefix) { prefix_ = prefix; }
   void set_postfix(const char* postfix) { postfix_ = postfix; }
   void set_flags(logflags flags) { flags_ = flags; }

   // getters
   loglevel level() const { return level_; }
   std::ostream& stream() const { return *stream_; }
   const std::string& name() const { return name_; }
   const char*prefix() const { return prefix_; }
   const char*postfix() const { return postfix_; }
   logflags flags() const { return flags_; }

   // ------------------------------------------------------------
   // internal stub for C++ insertion (<<) operators
   friend class iostub;
   class iostub final
   {
   public:
      iostub() = delete;
      iostub(const iostub&) = default;
      iostub& operator=(const iostub&) = delete;
      iostub& operator=(const iostub&&) = delete;

      // constructor outputs prefix (time, name, level, etc)
      explicit iostub(loglevel lv, const Logger& log) : log_(log), trace_(lv >= log_.level_)
      {
         if (trace_) log_.header(lv, log_.stream());
#if defined(_MSC_VER) && defined(_DEBUG)
         std::stringstream debug;
         log_.header(lv, debug);
         ToDebugConsole(debug);
#endif
      }

      // insertion method
      template<typename T> iostub& operator<<(const T& v)
      {
         if (!trace_) return *this;
         log_.stream() << v;
#if defined(_MSC_VER) && defined(_DEBUG)
         std::stringstream debug;
         debug << v;
         ToDebugConsole(debug);
#endif
         return *this;
      }

      // destructor outputs postfix (usually newline)
      ~iostub() { if (trace_ && log_.postfix_) *this << log_.postfix_; }

   protected:
      const Logger& log_;
      bool trace_;
#if defined(_MSC_VER) && defined(_DEBUG)
      LOGGER_API void ToDebugConsole(const std::stringstream& s) const;
#endif
   };

   // ------------------------------------------------------------

   // loggers (will return with RVO!)
   iostub debug() const { return iostub(loglevel::Debug, *this); }
   iostub info() const { return iostub(loglevel::Info, *this); }
   iostub warn() const { return iostub(loglevel::Warning, *this); }
   iostub error() const { return iostub(loglevel::Error, *this); }
   iostub fatal() const { return iostub(loglevel::Fatal, *this); }


protected:

   LOGGER_API void header(loglevel lv, std::ostream& stream) const;

   const std::string name_;
   loglevel level_;
   std::ostream* stream_;
   const char* prefix_;
   const char* postfix_;
   logflags flags_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline static Logger& Log(const char* name) { return Logger::get_logger(name); }

// --------------------------------------------------------------------------------------------------------------------------------

template<> struct enable_bitmask_operators<Logger::logflags> { static const bool enable = true; };

// --------------------------------------------------------------------------------------------------------------------------------
