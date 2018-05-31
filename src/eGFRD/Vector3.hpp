#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <cmath>
#include <ostream>
#include <iomanip>
#include "DefsEgfrd.hpp"
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct Matrix4;

// --------------------------------------------------------------------------------------------------------------------------------

struct Vector3
{
   // constructor(s)

   Vector3() noexcept : x_(0), y_(0), z_(0) { }

   Vector3(double x, double y, double z) noexcept : x_(x), y_(y), z_(z) { }

   // operator overloads with lhs is this

   Vector3 operator+(const Vector3& rhs) const
   {
      return add(*this, rhs);
   }
   Vector3 operator-(const Vector3& rhs) const
   {
      return substract(*this, rhs);
   }
   Vector3 operator*(double rhs) const
   {
      return multiply(*this, rhs);
   }
   Vector3 operator/(double rhs) const
   {
      return divide(*this, rhs);
   }

   bool operator==(const Vector3& rhs) const
   {
      return x_ == rhs.x_ && y_ == rhs.y_ && z_ == rhs.z_;
   }

   bool operator!=(const Vector3& rhs) const
   {
      return !operator==(rhs);
   }

   Vector3 operator-() const
   {
      return Vector3(-x_, -y_, -z_);
   }

   // operator overloads with lhs and assign to this

   Vector3& operator+=(const Vector3& rhs)
   {
      *this = add(*this, rhs);
      return *this;
   }
   Vector3& operator-=(const Vector3& rhs)
   {
      *this = substract(*this, rhs);
      return *this;
   }
   Vector3& operator*=(double rhs)
   {
      *this = multiply(*this, rhs);
      return *this;
   }
   Vector3& operator/=(double rhs)
   {
      *this = divide(*this, rhs);
      return *this;
   }

   // utility

   double length() const
   {
      return length(*this);
   }
   Vector3 normal() const
   {
      return normal(*this);
   }
   Vector3 abs() const
   {
      return abs(*this);
   }

   Vector3 modulo(const double s) const
   {
      double x = std::fmod(x_, s);
      if (x != 0 && (x > 0) == (s < 0)) x += s;
      double y = std::fmod(y_, s);
      if (y != 0 && (y > 0) == (s < 0)) y += s;
      double z = std::fmod(z_, s);
      if (z != 0 && (z > 0) == (s < 0)) z += s;
      return Vector3(x, y, z);
   }

   Vector3 modulo(const Vector3& s) const
   {
      double x = std::fmod(x_, s.X());
      if (x != 0 && (x > 0) == (s.X() < 0)) x += s.X();
      double y = std::fmod(y_, s.Y());
      if (y != 0 && (y > 0) == (s.Y() < 0)) y += s.Y();
      double z = std::fmod(z_, s.Z());
      if (z != 0 && (z > 0) == (s.Z() < 0)) z += s.Z();
      return Vector3(x, y, z);
   }

   // static methods

   static Vector3 add(const Vector3& u, const Vector3& v)
   {
      return Vector3(u.x_ + v.x_, u.y_ + v.y_, u.z_ + v.z_);
   }
   static Vector3 substract(const Vector3& u, const Vector3& v)
   {
      return Vector3(u.x_ - v.x_, u.y_ - v.y_, u.z_ - v.z_);
   }
   static Vector3 multiply(const Vector3& u, double s)
   {
      return Vector3(u.x_ * s, u.y_ *s, u.z_ * s);
   }
   static Vector3 divide(const Vector3& u, double s)
   {
      return Vector3(u.x_ / s, u.y_ / s, u.z_ / s);
   }
   static double length(const Vector3& u)
   {
      return std::sqrt(dot(u, u));
   }
   static Vector3 normal(const Vector3& u)
   {
      return divide(u, u.length());
   }
   static Vector3 abs(const Vector3& u)
   {
      return Vector3(std::abs(u.x_), std::abs(u.y_), std::abs(u.z_));
   }
   static double dot(const Vector3& u, const Vector3& v)
   {
      return u.x_ * v.x_ + u.y_ * v.y_ + u.z_ * v.z_;
   }
   static Vector3 cross(const Vector3& u, const Vector3& v)
   {
      return Vector3(u.y_*v.z_ - u.z_*v.y_, u.z_*v.x_ - u.x_*v.z_, u.x_*v.y_ - u.y_ *v.x_);
   }

   // matrix struff
   GFRD_EXPORT static Vector3 transformVector(const Vector3& v, const Matrix4& mat);

   // static variables

   GFRD_EXPORT static const Vector3 null;        // the null vector
   GFRD_EXPORT static const Vector3 one;         // the one vector (1,1,1)
   GFRD_EXPORT static const Vector3 ux;          // the unit vector in x
   GFRD_EXPORT static const Vector3 uy;          // the unit vector in y
   GFRD_EXPORT static const Vector3 uz;          // the unit vector in z

   // utility functions
   
   static Vector3 random(RandomNumberGenerator& rng)         // random vector with unit length! theta (-1.0 .. + 1.0), phi (0..2pi), draws 2 random numbers
   {
      double cos_theta = rng.uniform(-1, 1);
      double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
      double phi = rng.uniform(0, 2 * M_PI);
      double sin_phi = std::sin(phi);
      double cos_phi = std::cos(phi);
      return Vector3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
   }
   
   // member variable access

   double X() const { return x_; }
   double Y() const { return y_; }
   double Z() const { return z_; }
   
   friend struct Matrix4; 
   friend class Persistence;

private:
   double x_, y_, z_;
};

// --------------------------------------------------------------------------------------------------------------------------------

// operator overloads with lhs double and rhs vector 

inline Vector3 operator*(double lhs, const Vector3& rhs)
{
   return Vector3::multiply(rhs, lhs);
}

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Vector3& v)
{
   stream << std::defaultfloat << std::setprecision(12) << "(" << v.X() << ", " << v.Y() << ", " << v.Z() << ")";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Vector3 >
   {
      size_t operator()(const Vector3& v) const
      {
         return hash<double>()(v.X()) + 3 * hash<double>()(v.Y()) + 11 * hash<double>()(v.Z());
      }
   };
}

// --------------------------------------------------------------------------------------------------------------------------------
