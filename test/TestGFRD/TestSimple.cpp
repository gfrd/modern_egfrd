#include <stdio.h>
#include <string>
#include "Single_test.hpp"
#include <EGFRDSimulator.hpp>
#include <CopyNumbers.hpp>
#include <gsl/gsl_errno.h>

// --------------------------------------------------------------------------------------------------------------------------------

int testSimple()
{
   gsl_set_error_handler(nullptr);  // TestGFRD need at least one reference to GSL (or GNU will optimize it away, while it still is required for libGreensFunctions)
   
   RandomNumberGenerator rng;
   ReactionRuleCollection rules;
   World w;

   rng.seed(1234);

   Model m;
   auto s1 = m.add_species_type(SpeciesType("P", m.get_def_structure_type_id(), 1e-12, 3e-9));

   w.initialize(1E-6, m);

   w.throwInParticles(s1, 60, rng, false);

   EGFRDSimulator s(w, rules, rng);
   while (s.time() < 1E-3)
   {
      if ((s.num_steps() % 1000) == 0)
      {
         std::cout << "*";
         //s.dump("D:\\dump2.log");
      }
      s.step();
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testCompare()
{
   RandomNumberGenerator rng;
   ReactionRuleCollection rules;
   World w;
   EGFRDSimulator s(w, rules, rng);

   rng.seed(918273645);

   Model m;
   auto s1 = m.add_species_type(SpeciesType("P", m.get_def_structure_type_id(), 1e-12, 3e-9));

   w.initialize(1E-6, m);

   const int particles = 9;
   for (int i = 0; i < particles; ++i)
      w.add_particle(s1, w.get_def_structure_id(), w.world_size() / 2 + Vector3((-particles / 2 + i) * w.world_size().X() / (particles + 2), 0, 0));

   TINYTEST_TIME_START();
   const uint steps = 15000;
   for (uint i = 0; i < steps; ++i)
   {
      if ((i % 1000) == 0)
      {
         std::cout << "*";
         //std::string name = make_string() << "dumpSingle_" << std::setfill('0') << std::setw(5) << s.num_steps() << ".log";
         //s.dump(name);
      }
      s.step();
   }
   TINYTEST_TIME_STOP("Sim15000 steps");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testPairs()
{
   RandomNumberGenerator rng;
   ReactionRuleCollection rules;
   World w;
   EGFRDSimulator s(w, rules, rng);

   rng.seed(918273645);

   Model m;
   auto s1 = m.add_species_type(SpeciesType("P", m.get_def_structure_type_id(), 1e-12, 1e-9));

   w.initialize(1E-6, m);

   const int particles = 9;
   for (int i = 0; i < particles; ++i)
   {
      w.add_particle(s1, w.get_def_structure_id(), w.world_size() / 2 + Vector3((-particles / 2 + i) * w.world_size().X() / (particles + 2), 0, 0));
      w.add_particle(s1, w.get_def_structure_id(), w.world_size() / 2 + Vector3((-particles / 2 + i) * w.world_size().X() / (particles + 2), 4e-9, 0));
   }

   TINYTEST_TIME_START();
   const uint steps = 15000;
   for (uint i = 0; i < steps; ++i)
   {
      if ((i % 1000) == 0)
      {
         std::cout << "*";
         //std::string name = make_string() << "dumpPair_" << std::setfill('0') << std::setw(5) << s.num_steps() << ".log";
         //s.dump(name);
      }
      s.step();
   }
   TINYTEST_TIME_STOP("Pair15000 steps");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testCopyNumbers()
{
   RandomNumberGenerator rng;
   ReactionRuleCollection rules;
   World w;

   Model m;
   auto s1 = m.add_species_type(SpeciesType("P", m.get_def_structure_type_id(), 1e-12, 1e-9));
   auto s2 = m.add_species_type(SpeciesType("Q", m.get_def_structure_type_id(), 1e-12, 1e-9));
   auto s3 = m.add_species_type(SpeciesType("R", m.get_def_structure_type_id(), 1e-12, 1e-9));

   w.initialize(1E-6, m);
   w.throwInParticles(s1, 20, rng);
   w.throwInParticles(s2, 20, rng);

   rules.add_reaction_rule(ReactionRule(s1, 1.0E6, std::vector < SpeciesTypeID > {s2}));
   rules.add_reaction_rule(ReactionRule(s2, 2.0E6, std::vector < SpeciesTypeID > {s3}));
   rules.add_reaction_rule(ReactionRule(s3, 1.0E6, std::vector < SpeciesTypeID > {s1}));

   CopyNumbersInst cn(w, 1E-6);
   EGFRDSimulator s(w, rules, rng);
   s.add_external_event(0, &cn);
   TINYTEST_TIME_START();
   while (s.time() < 1E-3)
   {
      if ((s.num_steps() % 1000) == 0) std::cout << "*";
      s.step();
   }
   TINYTEST_TIME_STOP("run CopyNumbers");

   return 1;
}
// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(SampleSimple);
TINYTEST_ADD_TEST(testSimple);
TINYTEST_ADD_TEST(testCompare);
TINYTEST_ADD_TEST(testPairs);
TINYTEST_ADD_TEST(testCopyNumbers);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
