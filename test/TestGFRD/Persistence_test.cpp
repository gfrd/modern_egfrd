#include "DefsEgfrd.hpp"
#include "Persistence_test.hpp"
#include "World.hpp"
#include "Persistence.hpp"
#include "randomNumberGenerator.hpp"
#include "ReactionRuleCollection.hpp"
#include <EGFRDSimulator.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

int testPersistence_StoreWorld()
{
   Model m;
   auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id()));
   auto sPlane = m.add_structure_type(StructureType("plane"));

   World w;
   w.initialize(1.0, m);
   w.add_particle(s1, w.get_def_structure_id(), Vector3(1, 2, 3));

   w.add_structure(std::make_shared<PlanarSurface>(PlanarSurface("plane", sPlane, w.get_def_structure_id(), Plane(Vector3(), Vector3::ux, Vector3::uy, 1.0, 1.0, false))));

   Persistence p;
   p.store("WorldState.bin");
   p.store_world(w);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testPersistence_RestoreWorld()
{
   World w;
   Persistence p;
   p.retreive("WorldState.bin");
   p.retreive_world(w);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testPersistence_RandomNumberGen()
{
   RandomNumberGenerator rng;
   double d = rng.uniform(0, 1);
   TINYTEST_ASSERT(d >= 0 && d <= 1);

   {
      Persistence p;
      p.store("RngState.bin");
      p.store_rng(rng);
   }

   std::vector<double> random1(10);
   for (int i = 0; i < 10; ++i)
      random1[i] = rng.uniform(0, 1);

   {
      Persistence p;
      p.retreive("RngState.bin");
      p.retreive_rng(rng);
   }

   std::vector<double> random2(10);
   for (int i = 0; i < 10; ++i)
      random2[i] = rng.uniform(0, 1);

   for (int i = 0; i < 10; ++i)
      TINYTEST_EQUAL(random1[i], random2[i]);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testPersistence_ReactionRules()
{
   Model m;
   auto sPlane = m.add_structure_type(StructureType("plane"));

   auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), 1e-12, 1e-9));
   auto s2 = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-12, 0.5e-9));
   auto s3 = m.add_species_type(SpeciesType("C", m.get_def_structure_type_id(), 1e-12, 0.25e-9));

   auto sA = m.add_species_type(SpeciesType("A", sPlane, 2e-12, 1e-9));
   auto sB = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-12, 3e-9));

   ReactionRuleCollection rrc1;
   rrc1.add_reaction_rule(ReactionRule(s1, 1.0E-6, std::vector < SpeciesTypeID > {s2}));
   rrc1.add_reaction_rule(ReactionRule(s2, 1.0E-6, std::vector < SpeciesTypeID > {s3}));
   rrc1.add_reaction_rule(ReactionRule(s3, 1.0E-6, std::vector < SpeciesTypeID > {s1}));
   rrc1.add_interaction_rule(InteractionRule(sB, sPlane, 0.2, std::vector < SpeciesTypeID > {sA}));

   {
      Persistence p;
      p.store("RulesState.bin");
      p.store_reactionrules(rrc1);
   }

   ReactionRuleCollection rrc2;
   {
      Persistence p;
      p.retreive("RulesState.bin");
      p.retreive_reactionrules(rrc2);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testPersistence_Simulator()
{
   RandomNumberGenerator rng;
   ReactionRuleCollection rules;
   World w;

   Model m;
   auto s1 = m.add_species_type(SpeciesType("P", m.get_def_structure_type_id(), 1e-12, 3e-9));
   auto s2 = m.add_species_type(SpeciesType("Q", m.get_def_structure_type_id(), 1e-12, 3e-9));

   rules.add_reaction_rule(ReactionRule(s1, 1.0, std::vector<SpeciesTypeID>{s2}));
   rules.add_reaction_rule(ReactionRule(s2, 1.0, std::vector<SpeciesTypeID>{s1}));

   w.initialize(1E-6, m);
   w.throwInParticles(s1, 60, rng, false);

   EGFRDSimulator s(w, rules, rng);
   while (s.time() < 1E-5)
      s.step();

   {
      TINYTEST_TIME_START();
      Persistence p;
      p.store("SimState.bin");
      p.store_egfrd(s);
      TINYTEST_TIME_STOP("Store Simulator State to file");
   }

   const int steps = 120;
   size_t saveatstep = s.num_steps();
   std::vector<double> simtimes;
   simtimes.reserve(steps);
   std::vector<std::pair<bool, Vector3>> simParPos;
   simParPos.reserve(steps);

   for (int i = 0; i < steps; i++)
   {
      s.step();
      simtimes.emplace_back(s.time());

      bool found;
      auto pos = w.get_particle(ParticleID(i % 60), found);
      simParPos.emplace_back(std::make_pair(found, pos.second.position()));
   }

   {
      TINYTEST_TIME_START();
      Persistence p;
      p.retreive("SimState.bin");
      p.retreive_egfrd(s);
      TINYTEST_TIME_STOP("Retreive Simulator State from file");
   }

   TINYTEST_EQUAL_MSG(saveatstep, s.num_steps(), "Restored sim has different step.");
   for (int i = 0; i < steps; i++)
   {
      s.step();
      TINYTEST_EQUAL_MSG(simtimes[i], s.time(), "Continue sim has different time.");

      bool found;
      auto pos = w.get_particle(ParticleID(i % 60), found);
      TINYTEST_EQUAL_MSG(simParPos[i].first, found, "Continue sim has different particles.");
      if (found) TINYTEST_ASSERT_MSG(pos.second.position() == simParPos[i].second, "Continue sim has different particles.");
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(SimPersistance);
TINYTEST_ADD_TEST(testPersistence_StoreWorld);
TINYTEST_ADD_TEST(testPersistence_RestoreWorld);
TINYTEST_ADD_TEST(testPersistence_RandomNumberGen);
TINYTEST_ADD_TEST(testPersistence_ReactionRules);
TINYTEST_ADD_TEST(testPersistence_Simulator);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
