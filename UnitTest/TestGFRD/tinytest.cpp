// --------------------------------------------------------------------------------------------------------------------------------

#include "../common/tinytest.h"
#include "../common/tinytest_specific.cpp"

// --------------------------------------------------------------------------------------------------------------------------------

#include "TestSimple.hpp"
#include "EventScheduler_test.hpp"
#include "Vector2_test.hpp"
#include "Vector3_test.hpp"
#include "Box_test.hpp"
#include "Cylinder_test.hpp"
#include "Disk_test.hpp"
#include "Plane_test.hpp"
#include "Sphere_test.hpp"
#include "BDSimulator_test.hpp"
#include "Single_test.hpp"
#include "MatrixSpace_test.hpp"
#include "World_test.hpp"
#include "Model_test.hpp"
#include "Identifier_test.hpp"
#include "ParticleContainer_test.hpp"
#include "ParticleSimulator_test.hpp"
#include "GenericIterator_test.hpp"
#include "Persistence_test.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_MAIN();

/* Samples */
TINYTEST_RUN_SUITE(SampleSimple);

/* Scheduler */
TINYTEST_RUN_SUITE(EventScheduler);

/* Math */
TINYTEST_RUN_SUITE(Vector2);
TINYTEST_RUN_SUITE(Vector3);

/* Model */
TINYTEST_RUN_SUITE(Model);
TINYTEST_RUN_SUITE(Identifier);

/* Shapes */
TINYTEST_RUN_SUITE(Box);
TINYTEST_RUN_SUITE(Cylinder);
TINYTEST_RUN_SUITE(Disk);
TINYTEST_RUN_SUITE(Sphere);
TINYTEST_RUN_SUITE(Plane);

/* Simulator */
TINYTEST_RUN_SUITE(BDSimulator);
TINYTEST_RUN_SUITE(MatrixSpace);
TINYTEST_RUN_SUITE(ParticleContainer);
TINYTEST_RUN_SUITE(ParticleSimulator);
TINYTEST_RUN_SUITE(Single);
TINYTEST_RUN_SUITE(World);

/* Utils */
TINYTEST_RUN_SUITE(GenericIterator);
TINYTEST_RUN_SUITE(SimPersistance);


TINYTEST_END_MAIN();

// --------------------------------------------------------------------------------------------------------------------------------
