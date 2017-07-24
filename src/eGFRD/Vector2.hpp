#ifndef VECTOR2_HPP
#define VECTOR2_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <cmath>
#include <ostream>
#include <iomanip>
#include "DefsEgfrd.hpp"
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct Vector2
{
   // constructor(s)

   Vector2() noexcept : x_(0), y_(0) { }

   Vector2(double x, double y) noexcept : x_(x), y_(y) { }


   // operator overloads with lhs is this

   Vector2 operator+(const Vector2& rhs) const
   {
      return add(*this, rhs);
   }

   Vector2 operator-(const Vector2& rhs) const
   {
      return substract(*this, rhs);
   }

   Vector2 operator*(double rhs) const
   {
      return multiply(*this, rhs);
   }

   Vector2 operator/(double rhs) const
   {
      return divide(*this, rhs);
   }

   bool operator==(const Vector2& rhs) const
   {
      return x_ == rhs.x_ && y_ == rhs.y_;
   }

   bool operator!=(const Vector2& rhs) const
   {
      return !operator==(rhs);
   }

   Vector2 operator-() const
   {
      return Vector2(-x_, -y_);
   }

   // operator overloads with lhs and assign to this

   Vector2& operator+=(const Vector2& rhs)
   {
      *this = add(*this, rhs);
      return *this;
   }
   Vector2& operator-=(const Vector2& rhs)
   {
      *this = substract(*this, rhs);
      return *this;
   }
   Vector2& operator*=(double rhs)
   {
      *this = multiply(*this, rhs);
      return *this;
   }
   Vector2& operator/=(double rhs)
   {
      *this = divide(*this, rhs);
      return *this;
   }

   // utility

   double length() const
   {
      return length(*this);
   }
   Vector2 normal() const
   {
      return normal(*this);
   }
   Vector2 abs() const
   {
      return abs(*this);
   }

   // static methods

   static Vector2 add(const Vector2& u, const Vector2& v)
   {
      return Vector2(u.x_ + v.x_, u.y_ + v.y_);
   }
   static Vector2 substract(const Vector2& u, const Vector2& v)
   {
      return Vector2(u.x_ - v.x_, u.y_ - v.y_);
   }
   static Vector2 multiply(const Vector2& u, double s)
   {
      return Vector2(u.x_ * s, u.y_ *s);
   }
   static Vector2 divide(const Vector2& u, double s)
   {
      return Vector2(u.x_ / s, u.y_ / s);
   }
   static double length(const Vector2& u)
   {
      return std::sqrt(dot(u, u));
   }
   static Vector2 normal(const Vector2& u)
   {
      return divide(u, u.length());
   }
   static Vector2 abs(const Vector2& u)
   {
      return Vector2(std::abs(u.x_), std::abs(u.y_));
   }
   static double dot(const Vector2& u, const Vector2& v)
   {
      return u.x_ * v.x_ + u.y_ * v.y_;
   }

   // static variables

   GFRD_EXPORT static const Vector2 null;        // the null vector
   GFRD_EXPORT static const Vector2 ux;          // the unit vector in x
   GFRD_EXPORT static const Vector2 uy;          // the unit vector in y


   // utility functions

   static Vector2 random(RandomNumberGenerator& rng)                // random vector with unit length! phi (0..2pi), draws 1 random number
   {
      double phi = rng.uniform(0, 2 * M_PI);
      double sin_phi = std::sin(phi);
      double cos_phi = std::cos(phi);
      return Vector2(sin_phi, cos_phi);
   }

   // member variable access

   double X() const { return x_; }
   double Y() const { return y_; }

private:
   double x_, y_;
};

// --------------------------------------------------------------------------------------------------------------------------------

// operator overloads with lhs double and rhs vector 

inline Vector2 operator*(double lhs, const Vector2& rhs)
{
   return Vector2::multiply(rhs, lhs);
}

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Vector2& v)
{
   stream << std::setprecision(12) << "(" << v.X() << ", " << v.Y() << ")";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Vector2 >
   {
      size_t operator()(const Vector2& v) const
      {
         return hash<double>()(v.X()) + 3 * hash<double>()(v.Y());
      }
   };
}

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* Vector2_HPP */
