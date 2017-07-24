#include <iostream>
#include <ParticleSimulator.hpp>
#include "ParticleSimulator_test.hpp"
#include <EGFRDSimulator.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

int testParticleSimulator()
{
   Model m;
   auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), .3, 0.05));
   auto s2 = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), .3, 0.025));

   World w;
   w.initialize(10,m);

   TINYTEST_ALMOST_EQUALS(10.0, w.world_size().X());
   TINYTEST_ALMOST_EQUALS(10.0 / CompileConfigSimulator::MatrixCellsX, w.cell_size());
   TINYTEST_EQUAL(CompileConfigSimulator::MatrixCellsX, w.matrix_size()[0]);
   TINYTEST_EQUAL(CompileConfigSimulator::MatrixCellsY, w.matrix_size()[1]);
   TINYTEST_EQUAL(CompileConfigSimulator::MatrixCellsZ, w.matrix_size()[2]);

   RandomNumberGenerator rng;
   w.throwInParticles(s1, 200, rng, false, Vector3(1, 1, 1), Vector3(9, 9, 9));


   ReactionRuleCollection rules;
   rules.add_reaction_rule(ReactionRule(s1, 0.2, std::vector < SpeciesTypeID > {s2}));
   rules.add_reaction_rule(ReactionRule(s2, 0.2, std::vector < SpeciesTypeID > {s2}));

   ParticleSimulator ps(w, rules, rng);
   ps.set_repulsive();

   //std::cout << std::endl;
   //auto rrr = rules.get_rules();
   //for (auto& rs : rrr) 
   //{
   //   std::cout << " REACTION: " << rs.first << " RULES: " << std::endl;
   //   for (auto& rr : rs.second) std::cout << "\t" << rr << std::endl;
   //   std::cout << std::endl;
   //}

   auto r1 = rules.get_reaction_rules();
   TINYTEST_EQUAL(5, std::distance(r1.begin(), r1.end()));

   // again, but now no changes
   ps.set_repulsive();
   auto r2 = rules.get_reaction_rules();
   TINYTEST_EQUAL(5, std::distance(r2.begin(), r2.end()));

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testEGFRDSimulator()
{
   Model m;
   auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), .3, 0.05));
   auto s2 = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), .3, 0.025));

   World w;
   w.initialize(10,m);

   RandomNumberGenerator rng;
   w.throwInParticles(s1, 20, rng, false, Vector3(1, 1, 1), Vector3(9, 9, 9));


   ReactionRuleCollection rules;
   rules.add_reaction_rule(ReactionRule(s1, 0.2, std::vector < SpeciesTypeID > {s2}));
   rules.add_reaction_rule(ReactionRule(s2, 0.2, std::vector < SpeciesTypeID > {s2}));

   EGFRDSimulator s(w, rules, rng);
   while (s.time() == 0.0)
      s.step();

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(ParticleSimulator);
TINYTEST_ADD_TEST(testParticleSimulator);
TINYTEST_ADD_TEST(testEGFRDSimulator);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------