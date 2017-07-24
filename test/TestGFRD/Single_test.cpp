#include "Sphere.hpp"
#include "Single.hpp"
#include "Single_test.hpp"


// --------------------------------------------------------------------------------------------------------------------------------

int testSingleCreate()
{
   DomainID did;
   auto sid_pair = std::make_pair(ShellID(), Shell(did, Sphere(), Shell::Code::INIT));
   auto pid_pair = std::make_pair(ParticleID(), Particle());

   auto single = SingleSpherical(did, pid_pair, sid_pair, std::set<ReactionRule>());
   TINYTEST_STR_EQUAL("SingleSpherical", single.type_name());

   std::string info = make_string() << single;
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testSingleRandomPosition()
{
   RandomNumberGenerator rng;
   DomainID did;
   auto sid_pair = std::make_pair(ShellID(), Shell(did, Sphere(), Shell::Code::INIT));
   auto pid_pair = std::make_pair(ParticleID(), Particle(SpeciesTypeID(1), Sphere(Vector3(), 1E-9),StructureID(1), 1.0, 0.0));

   auto single = SingleSpherical(did, pid_pair, sid_pair, std::set<ReactionRule>());
   TINYTEST_STR_EQUAL("SingleSpherical", single.type_name());

   double r = 5E-9;     // 5 radii;
   auto displacement = single.create_position_vector(r, rng);     // new position at distance r from current position
   TINYTEST_ALMOST_EQUAL(displacement.length(), std::abs(r), 1E-8*pid_pair.second.radius());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testSingleMembers()
{
   //auto domain_id = DomainID(256);
   //auto pair_particle(std::make_pair(ParticleID(1024), Particle(SpeciesTypeID(2056), Sphere(Vector3(10, 9, 8), 5), StructureID(4112), 0.01)));
   //auto pair_shell = std::make_pair(ShellID(512), Shell(domain_id, Sphere(Vector3(1, 2, 3), 10)));
   //auto single = Single(domain_id, pair_particle, pair_shell);

   //TINYTEST_EQUAL(domain_id, single.domain().id());

   //TINYTEST_EQUAL_MSG(pair_particle.first, single.particle_id(), "particle id");
   //TINYTEST_EQUAL_MSG(pair_particle.second, single.particle(), "particle");

   //TINYTEST_EQUAL_MSG(pair_shell.first, single.shell_id(), "shell id");
   //TINYTEST_EQUAL_MSG(pair_shell.second, single.shell(), "shell");

   //TINYTEST_ALMOST_REL_EQUAL_MSG(5., single.mobility_radius(), 1E-15, "mobility radius");
   //TINYTEST_EQUAL_MSG(Vector3(1, 2, 3), single.position(), "position");

   //TINYTEST_ALMOST_REL_EQUAL_MSG(10., single.size(), 1E-15, "size");
   //TINYTEST_ASSERT_MSG(single.num_shells() == 1, "number of shells");
   //TINYTEST_ASSERT_MSG(single.multiplicity() == 1, "multiplicity");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testSinglePrint()
{
   //auto domain_id = DomainID(256);
   //auto pair_shell = std::make_pair(ShellID(512), Shell(domain_id, Sphere(Vector3(1, 2, 3), 10)));
   //auto pair_particle(std::make_pair(ParticleID(1024), Particle(SpeciesTypeID(2056), Sphere(Vector3(10, 9, 8), 10), StructureID(4112), 0.01)));
   //auto analytical_single = Single(domain_id, pair_particle, pair_shell);

   //std::stringstream sstream;
   //sstream << analytical_single;
   //TINYTEST_STR_EQUAL_MSG("Single{(DID(256),PID(1024)),dt=0}", sstream.str().c_str(), "print");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Single);
TINYTEST_ADD_TEST(testSingleCreate);
TINYTEST_ADD_TEST(testSingleRandomPosition);
TINYTEST_ADD_TEST(testSingleMembers);
TINYTEST_ADD_TEST(testSinglePrint);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------

