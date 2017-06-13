#ifndef IDENTIFIER_HPP
#define IDENTIFIER_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "DefsEgfrd.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TBase>
struct Identifier
{
   explicit Identifier(const idtype& value) : value_(value) { }

   idtype operator++() { return value_++; }

   const idtype& operator()() const { return value_; }

   bool operator!() const { return value_ == 0; }

   bool operator==(const TBase& rhs) const { return value_ == rhs.value_; }

   bool operator!=(const TBase& rhs) const { return value_ != rhs.value_; }

   bool operator<(const TBase& rhs) const { return value_ < rhs.value_; }

   bool operator>=(const TBase& rhs) const { return value_ >= rhs.value_; }

   bool operator>(const TBase& rhs) const { return value_ > rhs.value_; }

   bool operator<=(const TBase& rhs) const  { return value_ <= rhs.value_; }

   explicit operator idtype() const { return value_; }

   operator bool() const { return value_!= 0; }

private:
   idtype value_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* IDENTIFIER_HPP */
