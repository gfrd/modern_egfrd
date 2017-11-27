#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <sstream>

// --------------------------------------------------------------------------------------------------------------------------------

class make_string
{
public:

   make_string() noexcept : stream_() {}
   
   // no copy/assign/move operations 
   make_string(const make_string&) = delete;
   make_string& operator=(const make_string&) = delete;
   make_string(const make_string&&) = delete;
   make_string& operator=(const make_string&&) = delete;

   operator std::string() const { return stream_.str(); }

   template<typename T>
   make_string& operator<<(const T& v) { stream_ << v; return *this; }

private:
   std::stringstream stream_;
};

// --------------------------------------------------------------------------------------------------------------------------------
