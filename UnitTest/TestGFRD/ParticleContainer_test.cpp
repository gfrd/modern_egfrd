#include <iostream>
#include <ParticleContainerImpl.hpp>
#include "ParticleContainer_test.hpp"


// --------------------------------------------------------------------------------------------------------------------------------

int testParticleContainer()
{
   ParticleContainerImpl pc;
   TINYTEST_EQUAL(pc.world_size().X(), 1.0);
   TINYTEST_EQUAL(pc.world_size().Y(), 1.0);
   TINYTEST_EQUAL(pc.world_size().Z(), 1.0);
   TINYTEST_EQUAL(pc.cell_size(), 1.0 / CompileConfigSimulator::MatrixCellsX);
   TINYTEST_EQUAL(pc.matrix_size()[0], CompileConfigSimulator::MatrixCellsX);
   TINYTEST_EQUAL(pc.matrix_size()[1], CompileConfigSimulator::MatrixCellsY);
   TINYTEST_EQUAL(pc.matrix_size()[2], CompileConfigSimulator::MatrixCellsZ);
   TINYTEST_EQUAL(pc.num_particles(), 0);

   ParticleID pid(10);
   ParticleID pidn(11);

   {
      auto ir(pc.update_particle(std::make_pair(pid, Particle(SpeciesTypeID(1), Sphere(Vector3(0.2, 0.6, 0.4), 0.05), StructureID(0), 1E-5))));
      TINYTEST_EQUAL(true, ir);
   }
   {
      auto ir(pc.update_particle(std::make_pair(pid, Particle(SpeciesTypeID(1), Sphere(Vector3(0.3, 0.7, 0.5), 0.05), StructureID(0), 1E-5, 123.4))));
      TINYTEST_EQUAL(false, ir);
   }

   TINYTEST_ASSERT(pc.has_particle(pid));
   TINYTEST_ASSERT(!pc.has_particle(pidn));

   auto pip = pc.get_particle(pid);
   TINYTEST_EQUAL(pip.first, pid);
   TINYTEST_ALMOST_EQUAL(pip.second.D(), 1E-5, 1E-15);
   TINYTEST_ALMOST_EQUAL(pip.second.radius(), 0.05, 1E-15);
   TINYTEST_ALMOST_EQUAL(pip.second.v(), 123.4, 1E-15);
   TINYTEST_EQUAL(pip.second.position(), Vector3(0.3, 0.7, 0.5));

   bool found;
   auto pip2 = pc.get_particle(pid, found);           // returns a copy! not a const ref to the real thing
   TINYTEST_ASSERT(found);
   TINYTEST_EQUAL(pip2.first, pid);
   TINYTEST_ALMOST_EQUAL(pip2.second.D(), 1E-5, 1E-15);
   TINYTEST_EQUAL(pip2.second.position(), Vector3(0.3, 0.7, 0.5));

   auto pip3 = pc.get_particle(pidn, found);
   TINYTEST_ASSERT(!found);
   TINYTEST_EQUAL(pip3.first(), 0);    // no id

   pip2.second.position() = Vector3(13.0, 14.0, 15.0);                    // modify local copy
   auto pip4 = pc.get_particle(pid);      // no change in particle container
   TINYTEST_ASSERT(pip2.second.position() != pip4.second.position());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testParticleContainerIterator()
{
   ParticleContainerImpl pc;
   TINYTEST_EQUAL(pc.world_size().X(), 1.0);
   TINYTEST_EQUAL(pc.world_size().Y(), 1.0);
   TINYTEST_EQUAL(pc.world_size().Z(), 1.0);
   TINYTEST_EQUAL(pc.cell_size(), 1.0 / CompileConfigSimulator::MatrixCellsX);
   TINYTEST_EQUAL(pc.matrix_size()[0], CompileConfigSimulator::MatrixCellsX);
   TINYTEST_EQUAL(pc.matrix_size()[1], CompileConfigSimulator::MatrixCellsY);
   TINYTEST_EQUAL(pc.matrix_size()[2], CompileConfigSimulator::MatrixCellsZ);
   TINYTEST_EQUAL(pc.num_particles(), 0);

   SerialIDGenerator<ParticleID> pgen;

   const int size = 10;
   for (int i = 0; i < size; i++)
   {
      double d = i / (size - 1.0);
      Vector3 pos(0.5 + 0.5*d*std::sin(2.0*M_PI*i / 100), 0.5 + 0.5*d*std::cos(2.0*M_PI*i / 100), d);    // spiral
      auto ir(pc.update_particle(std::make_pair(pgen(), Particle(SpeciesTypeID(1), Sphere(pos, 0.05), StructureID(0), 1E-5, 99.99))));
      TINYTEST_EQUAL(true, ir);
   }

   for (auto i : pc.get_particles())
   {
      const auto& particle = i.second;
      Particle* pp = const_cast<Particle*>(&particle);            // trick const_iterator into non_const, so we can modify the positions
      pp->position() = -pp->position();                           // not common practice, but just a check
      std::cout << particle << std::endl;
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testParticleContainerCopyConstruct()
{
   ParticleContainerImpl pc1;
   TINYTEST_EQUAL(pc1.world_size().X(), 1.0);
   TINYTEST_EQUAL(pc1.world_size().Y(), 1.0);
   TINYTEST_EQUAL(pc1.world_size().Z(), 1.0);
   TINYTEST_EQUAL(pc1.cell_size(), 1.0 / CompileConfigSimulator::MatrixCellsX);
   TINYTEST_EQUAL(pc1.matrix_size()[0], CompileConfigSimulator::MatrixCellsX);
   TINYTEST_EQUAL(pc1.matrix_size()[1], CompileConfigSimulator::MatrixCellsY);
   TINYTEST_EQUAL(pc1.matrix_size()[2], CompileConfigSimulator::MatrixCellsZ);
   TINYTEST_EQUAL(pc1.num_particles(), 0);

   ParticleID pid1(10);
   ParticleID pid2(11);
   ParticleID pid3(12);

   auto pc2 = pc1;                      // assignment
   ParticleContainerImpl pc3(pc1);      // copy construct

   {
      auto ir(pc1.update_particle(std::make_pair(pid1, Particle(SpeciesTypeID(1), Sphere(Vector3(0.2, 0.6, 0.4), 0.05), StructureID(0), 1E-5))));
      TINYTEST_EQUAL(true, ir);
   }
   {
      auto ir(pc2.update_particle(std::make_pair(pid2, Particle(SpeciesTypeID(1), Sphere(Vector3(0.2, 0.6, 0.4), 0.05), StructureID(0), 1E-5))));
      TINYTEST_EQUAL(true, ir);
   }
   {
      auto ir(pc3.update_particle(std::make_pair(pid3, Particle(SpeciesTypeID(1), Sphere(Vector3(0.2, 0.6, 0.4), 0.05), StructureID(0), 1E-5))));
      TINYTEST_EQUAL(true, ir);
   }

   TINYTEST_ASSERT(pc1.has_particle(pid1));
   TINYTEST_ASSERT(!pc1.has_particle(pid2));
   TINYTEST_ASSERT(!pc1.has_particle(pid3));

   TINYTEST_ASSERT(!pc2.has_particle(pid1));
   TINYTEST_ASSERT(pc2.has_particle(pid2));
   TINYTEST_ASSERT(!pc2.has_particle(pid3));

   TINYTEST_ASSERT(!pc3.has_particle(pid1));
   TINYTEST_ASSERT(!pc3.has_particle(pid2));
   TINYTEST_ASSERT(pc3.has_particle(pid3));

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(ParticleContainer);
TINYTEST_ADD_TEST(testParticleContainer);
TINYTEST_ADD_TEST(testParticleContainerIterator);
TINYTEST_ADD_TEST(testParticleContainerCopyConstruct);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------