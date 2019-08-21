#include <stdexcept>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include "helperFunctionsGf.hpp"
#include "SphericalBesselGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

static double hermite_interpolate(double x, double x0, double dx, double const* y_array)
{
   const double hinv = 1.0 / dx;
   const uint i = static_cast<uint>((x - x0) * hinv);
   const uint index = i * 2;

   const double x_lo = (x - x0) * hinv - i;
   const double x_hi = 1.0 - x_lo;
   const double y_lo = y_array[index];
   const double ydot_lo = y_array[index + 1];  // *dx;
   const double y_hi = y_array[index + 2];
   const double ydot_hi = y_array[index + 3];  // *dx; move to table

   return x_hi * x_hi * (y_lo + x_lo * (2 * y_lo + ydot_lo)) + x_lo * x_lo * (y_hi + x_hi * (2 * y_hi - ydot_hi));
}

// --------------------------------------------------------------------------------------------------------------------------------

const SphericalBesselGenerator& SphericalBesselGenerator::instance()
{
   static const SphericalBesselGenerator sbg_;
   return sbg_;
}

// --------------------------------------------------------------------------------------------------------------------------------

void LoadTable(std::ifstream &file, sbTable& table)
{
   uint N;
   file.read(reinterpret_cast<char *>(&N), sizeof(uint));
   table.N = N;
   file.read(reinterpret_cast<char *>(&table.start), sizeof(double));
   file.read(reinterpret_cast<char *>(&table.delta), sizeof(double));
   table.values = std::make_unique<double[]>(2 * N);
   for (uint i = 0; i < 2 * N; i++)
      file.read(reinterpret_cast<char *>(&table.values[i]), sizeof(double));
}

// --------------------------------------------------------------------------------------------------------------------------------

SphericalBesselGenerator::SphericalBesselGenerator()
{
   std::string filename = GetExecutablePath() + "SphericalBesselTable";
   if (!file_exists(filename)) filename = GetExecutablePath() + "../SphericalBesselTable";
   
   std::ifstream file(filename, std::ios::in | std::ios::binary);
   if (!file.is_open()) throw std::runtime_error("SphericalBesselTable not found!");

   file.read(reinterpret_cast<char *>(&minNj_), sizeof(uint));
   file.read(reinterpret_cast<char *>(&maxNj_), sizeof(uint));
   file.read(reinterpret_cast<char *>(&minNy_), sizeof(uint));
   file.read(reinterpret_cast<char *>(&maxNy_), sizeof(uint));

   tableJ_ = std::make_unique<sbTable[]>(maxNj_ + 1);
   for (uint n = 0; n <= maxNj_; n++)
   {
      if (n < minNj_) continue;
      LoadTable(file, tableJ_[n]);
   }

   tableY_ = std::make_unique<sbTable[]>(maxNy_ + 1);
   for (uint n = 0; n <= maxNy_; n++)
   {
      if (n < minNy_) continue;
      LoadTable(file, tableY_[n]);
   }

   file.close();
}

// --------------------------------------------------------------------------------------------------------------------------------

inline double _j_smalln(uint n, double z)
{
   ASSERT(n <= 3 && n >= 0);

   if (n == 0) return (z != 0.0) ? sin(z) / z : 1.0;
   if (z == 0.0) return 0.0;

   const double z_r(1. / z);
   switch (n)
   {
   case 1:
      return (sin(z) * z_r - cos(z)) * z_r;
   case 2:
   {
      const double _3_zsq(3. * z_r * z_r);
      return (_3_zsq - 1) * sin(z) * z_r - _3_zsq * cos(z);
   }
   case 3:
   {
      const double _15_zsq(15. * z_r * z_r);
      return ((_15_zsq - 6.) * sin(z) * z_r - (_15_zsq - 1) * cos(z)) * z_r;
   }
   default:
      throw std::invalid_argument("0 <= n <= 3");
   }
}

inline double _y_smalln(uint n, double z)
{
   ASSERT(n <= 3 && n >= 0);
   if (z == 0) throw  std::invalid_argument("z=0 is not defined");

   const double z_r(1. / z);
   switch (n)
   {
   case 0:
      return -cos(z) * z_r;
   case 1:
      return -(cos(z) * z_r + sin(z)) * z_r;
   case 2:
   {
      const double _3_zsq(3. * z_r * z_r);
      return (1 - _3_zsq) * cos(z) * z_r - _3_zsq * sin(z);
   }
   case 3:
   {
      const double _15_zsq(15. * z_r * z_r);
      return ((-_15_zsq + 6.) * cos(z) * z_r - (_15_zsq - 1) * sin(z)) * z_r;
   }
   default:
      throw std::invalid_argument("0 <= n <= 2");
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

double SphericalBesselGenerator::j(uint n, double z) const
{
   if (n <= 3) return _j_smalln(n, z);
   if (n > MaxNj() || GfCfg.NO_BESSEL_TABLE) return gsl_sf_bessel_jl(n, z);

   const sbTable& table = tableJ_[n];

   const double minz(table.start + table.delta * 2);
   const double maxz(table.start + table.delta * (table.N - 2));        // could work with -1 but due to rounding its safer to use -2

   if (z >= minz && z < maxz)
   {
      const uint i = static_cast<uint>((z - table.start) / table.delta);
      ASSERT(i + 1 < table.N);
      return hermite_interpolate(z, table.start, table.delta, table.values.get());
   }

   return gsl_sf_bessel_jl(n, z);
}

double SphericalBesselGenerator::y(const uint n, const double z) const
{
   if (n <= 3) return _y_smalln(n, z);
   if (n > MaxNy() || GfCfg.NO_BESSEL_TABLE) return gsl_sf_bessel_yl(n, z);

   const sbTable& table = tableY_[n];
   const double minz(table.start + table.delta * 2);
   const double maxz(table.start + table.delta * (table.N - 2));        // could work with -1 but due to rounding its safer to use -2

   if (z >= minz && z < maxz)
   {
      const uint i = static_cast<uint>((z - table.start) / table.delta);
      ASSERT(i + 1 < table.N);
      return hermite_interpolate(z, table.start, table.delta, table.values.get());
   }

   return gsl_sf_bessel_yl(n, z);
}

// --------------------------------------------------------------------------------------------------------------------------------
