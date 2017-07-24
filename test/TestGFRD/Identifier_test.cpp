#include "ParticleID.hpp"
#include "SerialIDGenerator.hpp"
#include "Identifier_test.hpp"
#include <DomainID.hpp>
#include <Particle.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

int testIdentifiersConstruct()
{
   ParticleID id1(10);
   TINYTEST_EQUAL(id1(), 10);

   ParticleID id2(20);
   TINYTEST_EQUAL(id2(), 20);

   ParticleID id3 = id1;               // copy constructable
   TINYTEST_EQUAL(id3(), 10);
   TINYTEST_ASSERT(id3() < id2());
   TINYTEST_ASSERT(id2() > id1());

   ParticleID id4;
   TINYTEST_EQUAL(id4(), 0);

   id4 = id2;                          // assignable
   TINYTEST_EQUAL(id2(), id4());

   // id4++;                  // no math!
   // id4 = id1 + id2;
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testIdentifierGenerator()
{
   SerialIDGenerator<DomainID> dgen;

   auto did1 = dgen();
   TINYTEST_ASSERT(did1());            // first nonzero
   TINYTEST_EQUAL(did1(), 1);
   auto did2 = dgen();
   TINYTEST_ASSERT(did2());
   TINYTEST_ASSERT(did2() > did1());      // incremental

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testIdentifierPairs()
{
   using particle_id_pair = std::pair<ParticleID, Particle>;
   using cparticle_id_pair = std::pair<const ParticleID, const Particle&>;
   SerialIDGenerator<ParticleID> pgen;

   particle_id_pair pip1 = std::make_pair<ParticleID, Particle>(pgen(), Particle());
   cparticle_id_pair pip2 = pip1;         // pair to const pair assignment

   TINYTEST_EQUAL(pip1.first, pip2.first);
   TINYTEST_EQUAL(pip1.second, pip2.second);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Identifier);
TINYTEST_ADD_TEST(testIdentifiersConstruct);
TINYTEST_ADD_TEST(testIdentifierGenerator);
TINYTEST_ADD_TEST(testIdentifierPairs);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
