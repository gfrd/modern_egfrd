#ifndef MATRIX4_HPP
#define MATRIX4_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include <iomanip>
#include "DefsEgfrd.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct Vector3;

// --------------------------------------------------------------------------------------------------------------------------------

struct Matrix4
{
   // constructor(s)

   Matrix4() noexcept : m00(0), m01(0), m02(0), m03(0), m10(0), m11(0), m12(0), m13(0), m20(0), m21(0), m22(0), m23(0), m30(0), m31(0), m32(0), m33(0) { }
   Matrix4(const double r00, const double r01, const double r02, const double r03, const double r10, const double r11, const double r12, const double r13,
      const double r20, const double r21, const double r22, const double r23, const double r30, const double r31, const double r32, const double r33) noexcept :
      m00(r00), m01(r01), m02(r02), m03(r03), m10(r10), m11(r11), m12(r12), m13(r13),
      m20(r20), m21(r21), m22(r22), m23(r23), m30(r30), m31(r31), m32(r32), m33(r33) { }


   // operator overloads with lhs is this

   Matrix4 operator+(const Matrix4& rhs) const
   {
      return add(*this, rhs);
   }
   Matrix4 operator-(const Matrix4& rhs) const
   {
      return substract(*this, rhs);
   }
   Matrix4 operator*(double rhs) const
   {
      return multiply(*this, rhs);
   }
   Matrix4 operator/(double rhs) const
   {
      return divide(*this, rhs);
   }

   //bool operator==(const Matrix4& rhs) const
   //{
   //   return x_ == rhs.x_ && y_ == rhs.y_ && z_ == rhs.z_;
   //}

   //bool operator!=(const Matrix4& rhs) const
   //{
   //   return !operator==(rhs);
   //}

   Matrix4 operator-() const
   {
      return Matrix4(-m00, -m01, -m02, -m03,
         -m10, -m11, -m12, -m13,
         -m20, -m21, -m22, -m23,
         -m30, -m31, -m32, -m33);
   }

   //// operator overloads with lhs and assign to this

   Matrix4& operator+=(const Matrix4& rhs)
   {
      *this = add(*this, rhs);
      return *this;
   }
   Matrix4& operator-=(const Matrix4& rhs)
   {
      *this = substract(*this, rhs);
      return *this;
   }
   Matrix4& operator*=(double rhs)
   {
      *this = multiply(*this, rhs);
      return *this;
   }
   Matrix4& operator/=(double rhs)
   {
      *this = divide(*this, rhs);
      return *this;
   }

   //// utility

   //inline double length() const
   //{
   //   return length(*this);
   //}
   //inline Matrix4 normal() const
   //{
   //   return normal(*this);
   //}
   //inline Matrix4 abs() const
   //{
   //   return abs(*this);
   //}

   GFRD_EXPORT Vector3 multiply(const Vector3& v) const;

   Matrix4 squared() const
   {
      return multiply(*this, *this);
   }

   Matrix4 transpose() const
   {
      return transpose(*this);
   }

   double* getGLMatrix(double *mm) const
   {
      mm[0] = m00; mm[1] = m10; mm[2] = m20; mm[3] = m30;
      mm[4] = m01; mm[5] = m11; mm[6] = m21; mm[7] = m31;
      mm[8] = m02; mm[9] = m12; mm[10] = m22; mm[11] = m32;
      mm[12] = m03; mm[13] = m13; mm[14] = m23; mm[15] = m33;
      return mm;
   }


   //Matrix4 modulo(const double s) const
   //{
   //   double x = std::fmod(x_, s);
   //   if (x != 0 && (x > 0) == (s < 0)) x += s;
   //   double y = std::fmod(y_, s);
   //   if (y != 0 && (y > 0) == (s < 0)) y += s;
   //   double z = std::fmod(z_, s);
   //   if (z != 0 && (z > 0) == (s < 0)) z += s;
   //   return Matrix4(x, y, z);
   //}

   // static methods

   static Matrix4 add(const Matrix4& u, const Matrix4& v)
   {
      return Matrix4(
         u.m00 + v.m00, u.m01 + v.m01, u.m02 + v.m02, u.m03 + v.m03,
         u.m10 + v.m10, u.m11 + v.m11, u.m12 + v.m12, u.m13 + v.m13,
         u.m20 + v.m20, u.m21 + v.m21, u.m22 + v.m22, u.m23 + v.m23,
         u.m30 + v.m30, u.m31 + v.m31, u.m32 + v.m32, u.m33 + v.m33);
   }
   static Matrix4 substract(const Matrix4& u, const Matrix4& v)
   {
      return Matrix4(
         u.m00 - v.m00, u.m01 - v.m01, u.m02 - v.m02, u.m03 - v.m03,
         u.m10 - v.m10, u.m11 - v.m11, u.m12 - v.m12, u.m13 - v.m13,
         u.m20 - v.m20, u.m21 - v.m21, u.m22 - v.m22, u.m23 - v.m23,
         u.m30 - v.m30, u.m31 - v.m31, u.m32 - v.m32, u.m33 - v.m33);
   }
   static Matrix4 multiply(const Matrix4& u, double s)
   {
      return Matrix4(
         u.m00 * s, u.m01 * s, u.m02 * s, u.m03 * s,
         u.m10 * s, u.m11 * s, u.m12 * s, u.m13 * s,
         u.m20 * s, u.m21 * s, u.m22 * s, u.m23 * s,
         u.m30 * s, u.m31 * s, u.m32 * s, u.m33 * s);
   }
   static Matrix4 divide(const Matrix4& u, double s)
   {
      return Matrix4(
         u.m00 / s, u.m01 / s, u.m02 / s, u.m03 / s,
         u.m10 / s, u.m11 / s, u.m12 / s, u.m13 / s,
         u.m20 / s, u.m21 / s, u.m22 / s, u.m23 / s,
         u.m30 / s, u.m31 / s, u.m32 / s, u.m33 / s);
   }
   static Matrix4 transpose(const Matrix4& u)
   {
      return Matrix4(
         u.m00 , u.m10 , u.m20 , u.m30 ,
         u.m01 , u.m11 , u.m21 , u.m31 ,
         u.m02 , u.m12 , u.m22 , u.m32 ,
         u.m03 , u.m13 , u.m23 , u.m33 );
   }

   //static double length(const Matrix4& u)
   //{
   //   return std::sqrt(dot(u, u));
   //}
   //static Matrix4 normal(const Matrix4& u)
   //{
   //   return divide(u, u.length());
   //}
   //static Matrix4 abs(const Matrix4& u)
   //{
   //   return Matrix4(std::abs(u.x_), std::abs(u.y_), std::abs(u.z_));
   //}
   //static double dot(const Matrix4& u, const Matrix4& v)
   //{
   //   return u.x_ * v.x_ + u.y_ * v.y_ + u.z_ * v.z_;
   //}
   //static Matrix4 cross(const Matrix4& u, const Matrix4& v)
   //{
   //   return Matrix4(u.y_*v.z_ - u.z_*v.y_, u.z_*v.x_ - u.x_*v.z_, u.x_*v.y_ - u.y_ *v.x_);
   //}

   // matrix stuff

   static Matrix4 createRotationX(double theta)
   {
      return Matrix4(1, 0, 0, 0,
         0, std::cos(theta), -std::sin(theta), 0,
         0, std::sin(theta), std::cos(theta), 0,
         0, 0, 0, 1);
   }

   static Matrix4 createRotationY(double theta)
   {
      return Matrix4(std::cos(theta), 0, std::sin(theta), 0,
         0, 1, 0, 0,
         -std::sin(theta), 0, std::cos(theta), 0,
         0, 0, 0, 1);
   }

   static Matrix4 createRotationZ(double theta)
   {
      return Matrix4(std::cos(theta), -std::sin(theta), 0, 0,
         std::sin(theta), std::cos(theta), 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1);
   }

   GFRD_EXPORT static Matrix4 createRotationA(double theta, const Vector3& axis);

   GFRD_EXPORT static Matrix4 createRotationAB(const Vector3& a, const Vector3& b);

   GFRD_EXPORT static Matrix4 createTranslate(const Vector3& t);
   
   GFRD_EXPORT static Matrix4 createScale(const Vector3& t);

   GFRD_EXPORT static Matrix4 createCrossProduct(const Vector3& a);
   
   GFRD_EXPORT static Matrix4 createTensorProduct(const Vector3& a, const Vector3& b);

   GFRD_EXPORT static Vector3 multiply(const Matrix4& mat, const Vector3& v);

   GFRD_EXPORT static Matrix4 multiply(const Matrix4& a, const Matrix4& b);

   // static variables

   GFRD_EXPORT static const Matrix4 null;          // the null matrix
   GFRD_EXPORT static const Matrix4 identity;      // the identity matrix

   friend struct Vector3;

private:
   double m00, m01, m02, m03,
      m10, m11, m12, m13,
      m20, m21, m22, m23,
      m30, m31, m32, m33;
};

// --------------------------------------------------------------------------------------------------------------------------------

//inline std::ostream& operator<<(std::ostream& stream, const Matrix4& v)
//{
//   stream << std::setprecision(12) << "(" << v.X() << ", " << v.Y() << ", " << v.Z() << ")";
//   return stream;
//}
//
//// --------------------------------------------------------------------------------------------------------------------------------
//
//namespace std {
//   template<>
//   struct hash < Matrix4 >
//   {
//      size_t operator()(const Matrix4& v) const
//      {
//         return hash<double>()(v.X()) + 3 * hash<double>()(v.Y()) + 11 * hash<double>()(v.Z());
//      }
//   };
//}

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* MATRIX4_HPP */
