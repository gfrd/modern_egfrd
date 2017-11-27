#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <exception>
#include <stdexcept>

// --------------------------------------------------------------------------------------------------------------------------------

class gfrd_exception : public std::runtime_error
{
public:
   gfrd_exception(const std::string& msg) : runtime_error(msg) { }
};

// --------------------------------------------------------------------------------------------------------------------------------

class illegal_state : public gfrd_exception
{
public:
   illegal_state(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class illegal_argument : public gfrd_exception
{
public:
   illegal_argument(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class not_found : public gfrd_exception
{
public:
   not_found(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class already_exists : public gfrd_exception
{
public:
   already_exists(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class unsupported : public gfrd_exception
{
public:
   unsupported(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class propagation_error : public gfrd_exception
{
public:
   propagation_error(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class illegal_propagation_attempt : public gfrd_exception
{
public:
   illegal_propagation_attempt(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class not_implemented : public gfrd_exception
{
public:
   not_implemented(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------

class no_space : public gfrd_exception
{
public:
   no_space(const std::string& msg) : gfrd_exception(msg) {}
};

// --------------------------------------------------------------------------------------------------------------------------------
