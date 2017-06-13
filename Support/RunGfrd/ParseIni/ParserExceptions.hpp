#ifndef PARSER_EXCEPTIONS_HPP
#define PARSER_EXCEPTIONS_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <stdexcept>

   // --------------------------------------------------------------------------------------------------------------------------------

class parser_exception : public std::runtime_error
{
public:
   explicit parser_exception(const std::string& msg) : runtime_error(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class illegal_section : public parser_exception
{
public:
   explicit illegal_section(const std::string& msg) : parser_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class illegal_section_key : public parser_exception
{
public:
   explicit illegal_section_key(const std::string& msg) : parser_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class illegal_section_value : public parser_exception
{
public:
   explicit illegal_section_value(const std::string& msg) : parser_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class illegal_matrix_size : public parser_exception
{
public:
   explicit illegal_matrix_size(const std::string& msg) : parser_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif
