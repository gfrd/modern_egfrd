// --------------------------------------------------------------------------------------------------------------------------------
// genBesselTables.cpp
// --------------------------------------------------------------------------------------------------------------------------------

#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <algorithm>

// --------------------------------------------------------------------------------------------------------------------------------

using uint = unsigned int;

// --------------------------------------------------------------------------------------------------------------------------------

double minz_j(uint) { return 0; }

double minz_y(uint n) { return std::max(1u, n); }     // start after asymptotic begin of curve, used to be 0.5 for spherical and 5.0 for cylindrical

double maxz_j(uint n)
{
   double z = (n * n + n + 1) / 1.221e-4;
   return z >= 1000 ? std::max(1000u, n * n) : z;
}

double maxz_y(uint n)
{
   double z = (n * n + n + 1) / 6.06e-6;
   return z >= 2000 ? std::max(2000u, n * n) : z;
}

// --------------------------------------------------------------------------------------------------------------------------------

// Derivative of bessel functions
double gsl_sf_bessel_jdotl(uint n, double z)
{
   if (z == 0) return 0.0;
   if (n == 0) return (z * std::cos(z) - std::sin(z)) / (z*z);
   return 0.5 * (gsl_sf_bessel_jl(n - 1, z) - (gsl_sf_bessel_jl(n, z) + z * gsl_sf_bessel_jl(n + 1, z)) / z);
}
double gsl_sf_bessel_ydotl(uint n, double z)
{
   if (z == 0) return 0.0;
   if (n == 0) return (z * std::sin(z) + std::cos(z)) / (z*z);
   return 0.5 * (gsl_sf_bessel_yl(n - 1, z) - (gsl_sf_bessel_yl(n, z) + z * gsl_sf_bessel_yl(n + 1, z)) / z);
}

double gsl_sf_bessel_Jdotn(uint n, double z)
{
   double p = 1.0;
   double s = gsl_sf_bessel_Jn(n - 1, z);
   for (int i = 1; i <= 2; ++i)
   {
      p = -1.0 * (p * (1 - i + 1)) / i;
      s += p * gsl_sf_bessel_Jn(n - 1 + i * 2, z);
   }
   return s / 2;
}

double gsl_sf_bessel_Ydotn(uint n, double z)
{
   double p = 1.0;
   double s = gsl_sf_bessel_Yn(n - 1, z);
   for (int i = 1; i <= 2; ++i)
   {
      p = -1.0 * (p * (1 - i + 1)) / i;
      s += p*gsl_sf_bessel_Yn(n - 1 + i * 2, z);
   }
   return s / 2;
}

// --------------------------------------------------------------------------------------------------------------------------------

void write_table(std::ofstream& file, uint N, double start, double delta)
{
   file.write(reinterpret_cast<const char*>(&N), sizeof(uint));
   file.write(reinterpret_cast<const char*>(&start), sizeof(double));
   file.write(reinterpret_cast<const char*>(&delta), sizeof(double));
}

// --------------------------------------------------------------------------------------------------------------------------------

void BuildSphericalBesselTable(const char* filename)
{
   std::cout << "Spherical: " << filename;

   const uint minn_j = 4;
   const uint maxn_j = 51;
   const uint minn_y = 3;
   const uint maxn_y = 40; // GSL always uses Olver asymptotic form for n > 40 
   const double delta = M_PI / 35;

   // Python/SciPy      <->   GSL
   // special.sph_jn(n, z)    gsl_sf_bessel_jl(n, z)           spherical!
   // special.sph_yn(n, z)    gsl_sf_bessel_yl(n, z)

   std::ofstream file(filename, std::ofstream::binary | std::ios::out);
   file.write(reinterpret_cast<const char*>(&minn_j), sizeof(uint));
   file.write(reinterpret_cast<const char*>(&maxn_j), sizeof(uint));
   file.write(reinterpret_cast<const char*>(&minn_y), sizeof(uint));
   file.write(reinterpret_cast<const char*>(&maxn_y), sizeof(uint));

   // Spherical
   //for (double z = 1; z < 2; z += 0.1)
   //   for (int n = 0; n < 10; ++n)
   //      printf("n= %d  z= %.1lf  j= %.12lf  y= %.12lf  jdot= %.12lf  ydot= %.12lf\r\n", n, z, gsl_sf_bessel_jl(n, z), gsl_sf_bessel_yl(n, z), gsl_sf_bessel_jdotl(n, z), gsl_sf_bessel_ydotl(n, z));

   // j
   uint size = 0;
   for (uint n = minn_j; n <= maxn_j; ++n) size = std::max(size, static_cast<uint>(std::ceil(maxz_j(n) / delta)));

   double **table = new double*[size];         // pre-calculate j table for z=0..maxz and n=0..maxn
   for (uint i = 0; i < size; ++i)
   {
      double z = i*delta;
      table[i] = new double[maxn_j + 2];
      gsl_sf_bessel_jl_array(maxn_j + 1, z, table[i]);
   }

   for (uint n = minn_j; n <= maxn_j; ++n)
   {
      uint from = static_cast<uint>(std::ceil(minz_j(n) / delta));
      uint to = static_cast<uint>(std::ceil(maxz_j(n) / delta));
      write_table(file, to - from, from * delta, delta);

      for (uint i = from; i < to; ++i)
      {
         double j = table[i][n];                                                         //   = gsl_sf_bessel_jl(n, z);
         double z = i*delta;
         double jdot = i == 0 ? 0.0 : 0.5 * (table[i][n - 1] - (j + z * table[i][n + 1]) / z);
#if defined(CHECK_RESULTS)
         //double jdot2 = gsl_sf_bessel_jdotl(n, z);       // precise
         //double aerr1 = fabs(jdot - jdot2);
         //if (aerr1 > 2 * DBL_EPSILON) printf("x");

         //double jdot0 = i == 0 ? 0.0 : table[i][n - 1] - (n + 1) * j / z;      // old version (one table lookup)
         //double aerr2 = fabs(jdot - jdot0);
         //if (aerr2 > 0) printf("x");
#endif

         double jdd = delta * jdot;
         file.write(reinterpret_cast<const char*>(&j), sizeof(double));
         file.write(reinterpret_cast<const char*>(&jdd), sizeof(double));
      }
   }
   for (uint i = 0; i < size; ++i) delete[] table[i];
   delete[] table;
   size = 0;


   // y
   for (uint n = minn_y; n <= maxn_y; ++n) size = std::max(size, static_cast<uint>(std::ceil(maxz_y(n) / delta)));
   table = new double*[size];         // pre-calculate y table for z=0..maxz and n=0..maxn

   for (uint i = 1; i < size; ++i)     // skip zero its infinity
   {
      double z = i*delta;
      table[i] = new double[maxn_y + 2];
      gsl_sf_bessel_yl_array(maxn_y + 1, z, table[i]);
   }

   for (uint n = minn_y; n <= maxn_y; ++n)
   {
      uint from = static_cast<uint>(std::ceil(minz_y(n) / delta));
      uint to = static_cast<uint>(std::ceil(maxz_y(n) / delta));
      write_table(file, to - from, from * delta, delta);

      for (uint i = from; i < to; ++i)
      {
         double y = table[i][n];                                                         //   = gsl_sf_bessel_yl(n, z);
         double z = i*delta;
         double ydot = i == 0 ? 0.0 : 0.5 * (table[i][n - 1] - (y + z * table[i][n + 1]) / z);

#if defined(CHECK_RESULTS)
         //double ydot2 = gsl_sf_bessel_ydotl(n, z);             // presice (althoug not always equal to wolfram alpha)
         //double aerr1 = fabs(ydot - ydot2);
         //if (aerr1 > 5 * DBL_EPSILON) printf("x");

         //double ydot0 = i == 0 ? 0.0 : table[i][n - 1] - (n + 1) * y / z;         // old version (one table lookup)
         //double aerr2 = fabs(ydot - ydot0);
         //if (aerr2 > 0) printf("x");
#endif
         double ydd = delta * ydot;

         file.write(reinterpret_cast<const char*>(&y), sizeof(double));
         file.write(reinterpret_cast<const char*>(&ydd), sizeof(double));
      }
   }

   for (uint i = 1; i < size; ++i) delete[] table[i];
   delete[] table;

   file.close();
   std::cout << std::endl;
}

// --------------------------------------------------------------------------------------------------------------------------------

void BuildCylindricalBesselTable(const char* filename)
{
   std::cout << "Cylindrical: " << filename;

   const uint minn_j = 0;
   const uint maxn_j = 50;
   const uint minn_y = 0;
   const uint maxn_y = 50;
   const double delta = M_PI / 35;

   // Python/SciPy      <->   GSL
   // special.jv(n, z)        gsl_sf_bessel_Jn(n, z)           cylindrical!
   // special.yv(n, z)        gsl_sf_bessel_Yn(n, z)

   std::ofstream file(filename, std::ofstream::binary | std::ios::out);
   file.write(reinterpret_cast<const char*>(&minn_j), sizeof(uint));
   file.write(reinterpret_cast<const char*>(&maxn_j), sizeof(uint));
   file.write(reinterpret_cast<const char*>(&minn_y), sizeof(uint));
   file.write(reinterpret_cast<const char*>(&maxn_y), sizeof(uint));

   // Cylindrical
   //for (double z = 1; z < 2; z += 0.1)
   //   for (int n = 0; n < 10; ++n)
   //      printf("n= %d  z= %.1lf  J= %.12lf  Y= %.12lf  Jdot= %.12lf  Ydot= %.12lf\r\n", n, z, gsl_sf_bessel_Jn(n, z), gsl_sf_bessel_Yn(n, z), gsl_sf_bessel_Jdotn(n, z), gsl_sf_bessel_Ydotn(n, z));

   // J
   uint size = 0;
   for (uint n = minn_j; n <= maxn_j; ++n) size = std::max(size, static_cast<uint>(std::ceil(maxz_j(n) / delta)));

   double **table = new double*[size];         // pre-calculate J table for z=0..maxz and n=0..maxn
   for (uint i = 0; i < size; ++i)
   {
      double z = i*delta;
      table[i] = new double[maxn_j + 1];
      gsl_sf_bessel_Jn_array(minn_j, maxn_j, z, table[i]);
   }

   for (uint n = minn_j; n <= maxn_j; ++n)
   {
      uint from = static_cast<uint>(std::ceil(minz_j(n) / delta));
      uint to = static_cast<uint>(std::ceil(maxz_j(n) / delta));
      write_table(file, to - from, from * delta, delta);
      //printf("n=%d, size=%d, from=%d\r\n", n, to - from, from);

      for (uint i = from; i < to; ++i)
      {
         double j = table[i][n];                                                                      //   = gsl_sf_bessel_Jn(n, z);
         double jdot = n < 40 ? ((n > 0 ? table[i][n - 1] : -table[i][1]) - table[i][n + 1]) / 2 :    // for large n the table values are inaccurate to calculate the derivative
            gsl_sf_bessel_Jdotn(n, i*delta);
         double jdd = delta * jdot;

         file.write(reinterpret_cast<const char*>(&j), sizeof(double));
         file.write(reinterpret_cast<const char*>(&jdd), sizeof(double));
      }
   }
   for (uint i = 0; i < size; ++i) delete[] table[i];
   delete[] table;
   size = 0;


   // Y
   for (uint n = minn_y; n <= maxn_y; ++n) size = std::max(size, static_cast<uint>(std::ceil(maxz_y(n) / delta)));
   table = new double*[size];         // pre-calculate y table for z=0..maxz and n=0..maxn

   for (uint i = 1; i < size; ++i)     // skip zero its infinity
   {
      double z = i*delta;
      table[i] = new double[maxn_y + 1];
      gsl_sf_bessel_Yn_array(minn_y, maxn_y, z, table[i]);
   }

   for (uint n = minn_y; n <= maxn_y; ++n)
   {
      uint from = static_cast<uint>(std::ceil(minz_y(n) / delta));
      uint to = static_cast<uint>(std::ceil(maxz_y(n) / delta));
      write_table(file, to - from, from * delta, delta);
      //printf("n=%d, size=%d, from=%d\r\n", n, to - from, from);

      for (uint i = from; i < to; ++i)
      {
         double y = table[i][n];                                                                      //   = gsl_sf_bessel_Yn(n, z);
         double ydot = n < 40 ? ((n > 0 ? table[i][n - 1] : -table[i][1]) - table[i][n + 1]) / 2 :    // for large n the table values are inaccurate to calculate the derivative
            gsl_sf_bessel_Ydotn(n, i*delta);
         double ydd = delta * ydot;

         file.write(reinterpret_cast<const char*>(&y), sizeof(double));
         file.write(reinterpret_cast<const char*>(&ydd), sizeof(double));
      }
   }

   for (uint i = 1; i < size; ++i) delete[] table[i];
   delete[] table;

   file.close();
   std::cout << std::endl;
}

// --------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   std::cout << "Generate Bessel function LookUp-table" << std::endl;

   bool print = true;
   for (int i = 0; i < argc; ++i)
   {
      if (strcmp(argv[i], "-s") == 0 && (i + 1 < argc))
      {
         BuildSphericalBesselTable(argv[++i]);
         print = false;
      }
      if (strcmp(argv[i], "-c") == 0 && (i + 1 < argc))
      {
         BuildCylindricalBesselTable(argv[++i]);
         print = false;
      }
   }

   if (print)
   {
      std::cout << " -c <filename>     for Cylindrical Bessel Table." << std::endl;
      std::cout << " -s <filename>     for Spherical Bessel Table." << std::endl;
   }

   return 0;
}

// --------------------------------------------------------------------------------------------------------------------------------
