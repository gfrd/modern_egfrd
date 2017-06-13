#include <iostream>
#include <BDSimulator.hpp>
#include <ReactionRule.hpp>
#include "BDSimulator_test.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int testBDSimulator()
{
   RandomNumberGenerator rng;
   const int n_steps = 100;

   Model m;
   auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), 0.3, 0.5));
   auto s2 = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 0.3, 0.6));

   ReactionRuleCollection rules;
   rules.add_reaction_rule(ReactionRule(s1, 0.2, std::vector <SpeciesTypeID> {s2}));
   rules.add_reaction_rule(ReactionRule(s2, 0.2, std::vector <SpeciesTypeID> {s1}));

   World world;
   world.initialize(10.0, m);
   world.throwInParticles(s1, 5, rng, false);
   world.throwInParticles(s2, 5, rng, false);

   BDSimulator bd_sim(world, rules, rng);
   bd_sim.set_reaction_length_factor(0.05, 0.01);
   for (auto i = n_steps; --i >= 0;) bd_sim.step();
   TINYTEST_ALMOST_EQUAL(n_steps * bd_sim.dt(), bd_sim.time(), 10E-8);
   
   bd_sim.step(2 * bd_sim.time());
   TINYTEST_ALMOST_EQUAL(2 * n_steps * bd_sim.dt(), bd_sim.time(), 10E-8);


   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(BDSimulator);
TINYTEST_ADD_TEST(testBDSimulator);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------