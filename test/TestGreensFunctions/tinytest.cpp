// --------------------------------------------------------------------------------------------------------------------------------

#include "../common/tinytest.h"

// --------------------------------------------------------------------------------------------------------------------------------

#include "SphericalBesselGenerator_test.hpp"
#include "CylindricalBesselGenerator_test.hpp"
#include "GreensFunction3D_test.hpp"
#include "GreensFunction3DAbs_test.hpp"
#include "GreensFunction3DAbsSym_test.hpp"
#include "GreensFunction3DRadAbs_test.hpp"
#include "GreensFunction3DRadInf_test.hpp"
#include "GslSumLevinU_test.hpp"
#include "Factorial_test.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_MAIN();

/* GreensFunctions */
TINYTEST_RUN_SUITE(GreensFunction3D);
TINYTEST_RUN_SUITE(GreensFunction3DAbs);
TINYTEST_RUN_SUITE(GreensFunction3DAbsSym);
TINYTEST_RUN_SUITE(GreensFunction3DRadAbs);
TINYTEST_RUN_SUITE(GreensFunction3DRadInf);

/* Math */
TINYTEST_RUN_SUITE(SphericalBesselGenerator);
TINYTEST_RUN_SUITE(CylindricalBesselGenerator);

/* Gsl */
TINYTEST_RUN_SUITE(Factorial);
TINYTEST_RUN_SUITE(GslSumLevinU);

TINYTEST_END_MAIN();

