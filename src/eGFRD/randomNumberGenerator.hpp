#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <random>

// --------------------------------------------------------------------------------------------------------------------------------

class RandomNumberGenerator
{
public:

   RandomNumberGenerator() : rng_() { }

   double normal(double mean, double sigma)
   {
      return std::normal_distribution<>(mean, sigma)(rng_);
   }

   double uniform(double min, double max)
   {
      return std::uniform_real_distribution<>(min, max)(rng_);
   }

   int uniform_int(int min, int max)
   {
      return std::uniform_int_distribution<>(min, max)(rng_);
   }

   double operator()()
   {
      return std::uniform_real_distribution<>()(rng_);
   }

   void seed(unsigned long int seed)
   {
      rng_.seed(seed);
   }

private:
   friend class Persistence;
   std::mt19937 rng_;
};

// --------------------------------------------------------------------------------------------------------------------------------
