#include <stdexcept>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>

#include "helperFunctionsGf.hpp"
#include "CylindricalBesselGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

static double hermite_interpolate(double x, double x0, double dx, double const* y_array)
{
   const double hinv = 1.0 / dx;
   const uint i = static_cast<uint>((x - x0) * hinv);
   const uint index = i * 2;

   const double x_lo = (x - x0) * hinv - i;
   const double x_hi = 1.0 - x_lo;
   const double y_lo = y_array[index];
   const double ydot_lo = y_array[index + 1];  // * dx;
   const double y_hi = y_array[index + 2];
   const double ydot_hi = y_array[index + 3];  // *dx; move to table

   return x_hi * x_hi * (y_lo + x_lo * (2 * y_lo + ydot_lo)) + x_lo * x_lo * (y_hi + x_hi * (2 * y_hi - ydot_hi));
}

// --------------------------------------------------------------------------------------------------------------------------------

const CylindricalBesselGenerator& CylindricalBesselGenerator::instance()
{
   static const CylindricalBesselGenerator cbg_;
   return cbg_;
}

// --------------------------------------------------------------------------------------------------------------------------------

void LoadTable(std::ifstream &file, cbTable& table)
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

CylindricalBesselGenerator::CylindricalBesselGenerator()
{
   std::string fileName(GetExecutablePath() + "CylindricalBesselTable.bin");

   std::ifstream file(fileName, std::ios::in | std::ios::binary);
   if (!file.is_open()) throw std::runtime_error("CylindricalBesselTable.bin not found!");

   file.read(reinterpret_cast<char *>(&minNJ_), sizeof(uint));
   file.read(reinterpret_cast<char *>(&maxNJ_), sizeof(uint));
   file.read(reinterpret_cast<char *>(&minNY_), sizeof(uint));
   file.read(reinterpret_cast<char *>(&maxNY_), sizeof(uint));

   tableJ_ = std::make_unique<cbTable[]>(maxNJ_ + 1);
   for (uint n = 0; n <= maxNJ_; n++)
   {
      if (n < minNJ_) continue;
      LoadTable(file, tableJ_[n]);
   }

   tableY_ = std::make_unique<cbTable[]>(maxNY_ + 1);
   for (uint n = 0; n <= maxNY_; n++)
   {
      if (n < minNY_) continue;
      LoadTable(file, tableY_[n]);
   }

   file.close();
}
// --------------------------------------------------------------------------------------------------------------------------------

double CylindricalBesselGenerator::J(uint n, double z) const
{
   if (n > MaxNJ() || GfCfg.NO_BESSEL_TABLE) return gsl_sf_bessel_Jn(n, z);

   const cbTable& table = tableJ_[n];
   const double minz(table.start);
   const double maxz(table.start + table.delta * (table.N - 2));

   if (z >= minz && z < maxz)
   {
      const uint i = static_cast<uint>((z - table.start) / table.delta);
      ASSERT(i + 1 < table.N);

      return hermite_interpolate(z, table.start, table.delta, table.values.get());
   }

   return gsl_sf_bessel_Jn(n, z);
}

double CylindricalBesselGenerator::Y(const uint n, const double z) const
{
   if (n > MaxNY() || GfCfg.NO_BESSEL_TABLE) return gsl_sf_bessel_Yn(n, z);

   const cbTable& table = tableY_[n];
   const double minz(table.start);
   const double maxz(table.start + table.delta * (table.N - 2));

   if (z >= minz && z < maxz)
   {
      const uint i = static_cast<uint>((z - table.start) / table.delta);
      ASSERT(i + 1 < table.N);

      return hermite_interpolate(z, table.start, table.delta, table.values.get());
   }

   return gsl_sf_bessel_Yn(n, z);
}

// --------------------------------------------------------------------------------------------------------------------------------
