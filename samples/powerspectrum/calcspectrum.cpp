// --------------------------------------------------------------------------------------------------------------------------------
//
// calcspectrum <power_rec.dat>
//
// where power_rec.dat is the output file of the eGFRD reaction-recorder
// that should contain 2 reactions (1=association, 2=dissociation)
//
// See: Kaizu K, de Ronde W, Paijmans J, Takahashi K, Tostevin F, ten Wolde PR(2014) The Berg - Purcell Limit Revisited Biophys J, 106 : 976 - 985.
//      doi(https://dx.doi.org/10.1016/j.bpj.2013.12.030 )
//
//
// compile:
// unx> g++ -std=c++11 calcspectrum.cpp -o calcspectrum
// win> cl calcspectrum.cpp
// --------------------------------------------------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <complex>

// --------------------------------------------------------------------------------------------------------------------------------

const int fourier_size = 1 << 15;
const int output_size = 1 << 8;
const int average_size = fourier_size / output_size;
const double PI = 3.1415926535897932384626433832795;
const double NA = 6.022140857E23;

// --------------------------------------------------------------------------------------------------------------------------------

void calc_spectrum(std::string filename)
{
   std::cout << "Fourier Size: " << fourier_size << ", Output Size: " << output_size << std::endl;

   // Take time
   auto start = std::chrono::high_resolution_clock::now();

   // pre-process file
   int   count = 0;
   int   offset = 0;
   double final_time;
   {
      std::ifstream infile(filename);
      std::string line, lastline;
      while (std::getline(infile, line))
      {
         if (line.find("Time") != std::string::npos) { offset = count + 1; count = 0; continue; }
         count++;
         lastline = line;
      }
      std::istringstream iss(lastline);
      iss >> final_time;
   }
   std::cout << "Events: " << count << ", Last: " << final_time << std::endl;

   // Calculate omega table
   auto omega_max = 10.0;
   auto omega_min = std::log10(2.0 * PI * omega_max / final_time);
   auto weight = (omega_max - omega_min) / fourier_size;
   std::cout << "Omega Min: " << std::fixed << std::setprecision(3) << omega_min << ", Max: " << omega_max << std::endl;
   std::vector<double> omega(fourier_size);
   for (int i = 0; i < fourier_size; ++i)
      omega[i] = std::pow(10, (i + 1.0) * weight + omega_min);

   // Calculate fourier
   std::vector<std::complex<double>> fourier(fourier_size);
   {
      std::ifstream infile(filename);
      std::string line;
      for (int i = 0; i < offset; ++i) std::getline(infile, line);
      for (int i = 0; i < count; ++i)
      {
         std::getline(infile, line);
         std::istringstream iss(line);
         double time; int rid;
         iss >> time >> rid;

         int flip_state = rid == 1 ? 1 : -1;
         for (int j = 0; j < fourier_size; j++)
            fourier[j] += std::polar(flip_state / omega[j], omega[j] * time);
      }
   }

   // Calculate power-spectrum
   std::vector<double> power(fourier_size);
   for (int i = 0; i < fourier_size; i++)
      power[i] = std::pow(std::abs(fourier[i]), 2) / final_time;

   // Reduce number of points by averaging
   std::vector<double> omega_avg(output_size);
   std::vector<double> power_avg(output_size);
   for (int i = 0; i < output_size; i++)
   {
      double sum = 0;
      for (int j = 0; j < average_size; j++)
         sum += power[average_size * i + j];
      power_avg[i] = sum / average_size;
      omega_avg[i] = omega[average_size * i + average_size / 2];
   }

   // Calculate intrinsic and effective
   double D = 1e-12;
   double r = 5e-9;
   double ka = 9.16639E-19;
   double kd = 220.8;
   double kD = 4 * PI * r * 2 * D;
   double kon = 1 / (1 / ka + 1 / kD);
   double koff = 1 / (1 / kd + 1 / kD);
   double c = 0.4E-3 * NA;
   double mu_int = ka * c + kd;
   double mu_eff = kon * c + koff;
   std::vector<double> intrinsic(output_size);
   std::vector<double> effective(output_size);
   for (int i = 0; i < output_size; i++)
   {
      intrinsic[i] = 0.5 * mu_int / (std::pow(mu_int, 2) + std::pow(omega_avg[i], 2));
      effective[i] = 0.5 * mu_eff / (std::pow(mu_eff, 2) + std::pow(omega_avg[i], 2));
   }

   // Print elapsed calculation time
   auto now = std::chrono::high_resolution_clock::now();
   auto elapsed = std::chrono::duration_cast<std::chrono::duration<float>>(now - start).count();
   std::cout << "Calculation Time: " << elapsed << " sec" << std::endl;

   // Print output table
   std::cout << std::endl << std::setw(20) << "omega" << std::setw(20) << "power" << std::setw(20) << "intrinsic" << std::setw(20) << "effective" << std::endl;
   std::cout << std::scientific << std::setprecision(12);
   for (int i = 0; i < output_size; ++i)
      std::cout << omega_avg[i] << std::setw(20) << power_avg[i] << std::setw(20) << intrinsic[i] << std::setw(20) << effective[i] << std::endl;
}

// --------------------------------------------------------------------------------------------------------------------------------

int main(int argc, const char** argv)
{
   calc_spectrum(std::string(argv[1]));
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------
