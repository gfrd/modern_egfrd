#include <sstream>
#include "DefsEgfrd.hpp"
#include "Model.hpp"
#include "Model_test.hpp"
#include "ReactionRuleCollection.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int testModelBasics()
{
   Model m;
   auto worldID = m.get_def_structure_type_id();

   SpeciesType s1("S", worldID, 1.5e-12, 5e-9);
   SpeciesType s2("P", worldID, 1e-12, 7e-9);

   m.add_species_type(std::move(s1));
   m.add_species_type(std::move(s2));

   ReactionRuleCollection rules;

   ReactionRule&& rr = ReactionRule(s1.id(), 0.2, std::vector<SpeciesTypeID>());
   rules.add_reaction_rule(rr);

   ReactionRule acopy = rr;

   TINYTEST_CHECK_THROW(rules.add_reaction_rule(acopy), already_exists);

   rules.add_reaction_rule(ReactionRule(s1.id(), 0.2, std::vector < SpeciesTypeID > {s2.id()}));

   TINYTEST_CHECK_THROW(rules.add_reaction_rule(ReactionRule(s1.id(), 0.2, std::vector < SpeciesTypeID > {s2.id()})), already_exists);

   rules.add_reaction_rule(ReactionRule(s1.id(), 0.2, std::vector < SpeciesTypeID > {s1.id(), s2.id()}));

   TINYTEST_CHECK_THROW(rules.add_reaction_rule(ReactionRule(s1.id(), 0.2, std::vector < SpeciesTypeID > {s1.id(), s2.id()})), already_exists);

   TINYTEST_CHECK_THROW(rules.add_reaction_rule(ReactionRule(s1.id(), 0.2, std::vector < SpeciesTypeID > {s2.id(), s1.id()})), already_exists);

   rules.add_reaction_rule(ReactionRule(s2.id(), 0.2, std::vector < SpeciesTypeID > {s2.id(), s1.id()}));

   TINYTEST_CHECK_THROW(rules.add_reaction_rule(ReactionRule(s2.id(), 0.2, std::vector < SpeciesTypeID > {s1.id(), s2.id()})), already_exists);

   rules.add_reaction_rule(ReactionRule(s1.id(), s2.id(), 0.2, std::vector < SpeciesTypeID > {s2.id(), s1.id() }));

   TINYTEST_CHECK_THROW(rules.add_reaction_rule(ReactionRule(ReactionRule::reactants(s1.id(), s2.id()), 0.2, std::vector < SpeciesTypeID > {s2.id(), s1.id()})), already_exists);

   TINYTEST_CHECK_THROW(rules.add_reaction_rule(ReactionRule(ReactionRule::reactants(s2.id(), s1.id()), 0.2, std::vector < SpeciesTypeID > {s1.id(), s2.id()})), already_exists);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testModelReactionRules()
{
   Model m;
   auto worldID = m.get_def_structure_type_id();

   SpeciesType s1("S", worldID, 1.5e-12, 5e-9);
   m.add_species_type(s1);
   SpeciesType s2("P", worldID, 1e-12, 7e-9);
   m.add_species_type(s2);

   ReactionRuleCollection rules;
   rules.add_reaction_rule(ReactionRule(s1.id(), 0.2, { s2.id() }));
   rules.add_reaction_rule(ReactionRule(s2.id(), 0.2, { s1.id() }));

   {
      auto gen(rules.query_reaction_rules(s1.id()));
      TINYTEST_ASSERT(gen.size() == 1);
      for (auto rr : gen)
      {
         TINYTEST_ASSERT(rr.get_reactants() == ReactionRule::reactants(s1.id()));
         TINYTEST_ASSERT(rr.get_products() == std::vector < SpeciesTypeID > { s2.id() });
      }
   }

   {
      auto gen(rules.query_reaction_rules(s1.id(), s2.id()));
      TINYTEST_ASSERT(gen.size() == 0);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testModelInteractionRules()
{
   Model m;
   auto worldID = m.get_def_structure_type_id();

   StructureType sA("membrane");
   m.add_structure_type(sA);

   SpeciesType sW("W", worldID, 1.5e-12, 5e-9);
   m.add_species_type(sW);
   SpeciesType sM("M", sA.id(), 1e-12, 7e-9);
   m.add_species_type(sM);

   ReactionRuleCollection rules;
   rules.add_reaction_rule(ReactionRule(sM.id(), 0.2, { sW.id() }));                   // membrane -> world
   rules.add_interaction_rule(InteractionRule(sW.id(), sA.id(), 0.2, { sM.id() }));    // world -> membrane


   {
      auto gen(rules.query_reaction_rules(sM.id()));
      TINYTEST_ASSERT(gen.size() == 1);
      for (auto rr : gen)
      {
         TINYTEST_ASSERT(rr.get_reactants() == ReactionRule::reactants(sM.id()));
         TINYTEST_ASSERT(rr.get_products() == std::vector < SpeciesTypeID > { sW.id() });
      }
   }

   {
      auto gen(rules.query_interaction_rules(sW.id(), sA.id()));
      TINYTEST_ASSERT(gen.size() == 1);
      for (auto ir : gen)
      {
         TINYTEST_ASSERT(ir.get_reactants() == InteractionRule::reactants(sW.id(), sA.id()));
         TINYTEST_ASSERT(ir.get_products() == std::vector < SpeciesTypeID > { sM.id() });
      }
   }

   {
      auto gen(rules.query_reaction_rules(sW.id(), sM.id()));
      TINYTEST_ASSERT(gen.size() == 0);
   }

   {
      auto gen(rules.query_reaction_rules(sM.id(), sW.id()));
      TINYTEST_ASSERT(gen.size() == 0);
   }

   {
      auto gen(rules.query_interaction_rules(sM.id(), sA.id()));
      TINYTEST_ASSERT(gen.size() == 0);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Model);
TINYTEST_ADD_TEST(testModelBasics);
TINYTEST_ADD_TEST(testModelReactionRules);
TINYTEST_ADD_TEST(testModelInteractionRules);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
