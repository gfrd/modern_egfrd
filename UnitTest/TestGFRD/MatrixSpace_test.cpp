#include "MatrixSpace.hpp"
#include "MatrixSpace_test.hpp"
#include <iostream>
#include <set>
#include <Sphere.hpp>
#include <Cylinder.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

#if !defined(M_PI)
#define M_PI                  3.1415926535897932384626433832795
#endif
#if defined(_DEBUG)
#define NUMBEROFTESTS             100
#define STEPSIZE                  0.01
#else
#define NUMBEROFTESTS             500
#define STEPSIZE                  0.001
#endif
#include "randomNumberGenerator.hpp"
#include <Particle.hpp>
#include <ParticleID.hpp>


// --------------------------------------------------------------------------------------------------------------------------------

template<typename TRng, typename TContainer>
inline void shuffle(TRng& rng, TContainer& cntnr)
{
   for (size_t i = cntnr.size(); i > 0;)
   {
      --i;
      const int j(rng.uniform_int(0, static_cast<int>(i)));
      std::swap(cntnr[i], cntnr[j]);
   }
}

template<typename TCollector>
struct collector_set
{
   collector_set() : result() {}

   void operator()(typename TCollector::const_iterator i, const Vector3&)
   {
      result.insert((*i).first);
   }

   std::set<typename TCollector::key_type> result;
};

template<typename TCollector>
struct collector_count
{
   collector_count() : count(0) {}

   void operator()(typename TCollector::const_iterator i, const Vector3&)
   {
      UNUSED(i);
      ++count;
   }

   int count;
};

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceSized()
{
   using ms_type = MatrixSpace<Sphere, int, 10, 10, 10>;

   ms_type ms;
   TINYTEST_EQUAL(ms.size(), 0);
   ms.update(std::make_pair(0, Sphere(Vector3(0.2, 0.6, 0.4), 0.05)));
   TINYTEST_EQUAL(ms.size(), 1);
   ms.update(std::make_pair(1, Sphere(Vector3(0.2, 0.6, 0.4), 0.05)));
   TINYTEST_EQUAL(ms.size(), 2);
   ms.update(std::make_pair(1, Sphere(Vector3(0.1, 0.2, 0.3), 0.05)));
   TINYTEST_EQUAL(ms.size(), 2);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceUpdate()
{
   using ms_type = MatrixSpace<Sphere, int, 10, 10, 10>;

   ms_type ms;
   TINYTEST_ALMOST_EQUAL(0.1, ms.cell_size(), 0.001);

   {
      std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(0, Sphere(Vector3(0.2, 0.6, 0.4), 0.05))));
      TINYTEST_EQUAL(true, ir.second);
      TINYTEST_ASSERT(ms.end() != ms.find(0));
      TINYTEST_ASSERT(ms.end() == ms.find(1));
   }
   {
      std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(0, Sphere(Vector3(0.2, 0.65, 0.4), 0.05))));
      TINYTEST_EQUAL(false, ir.second);
      TINYTEST_EQUAL(Sphere(Vector3(0.2, 0.65, 0.4), 0.05), (*ir.first).second);
      TINYTEST_ASSERT(ms.end() != ms.find(0));
      TINYTEST_ASSERT(ms.end() == ms.find(1));
   }
   {
      std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(0, Sphere(Vector3(0.2, 0.2, 0.4), 0.05))));
      TINYTEST_EQUAL(false, ir.second);
      TINYTEST_EQUAL(Sphere(Vector3(0.2, 0.2, 0.4), 0.05), (*ir.first).second);
      TINYTEST_ASSERT(ms.end() != ms.find(0));
      TINYTEST_ASSERT(ms.end() == ms.find(1));
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceErase()
{
   using ms_type = MatrixSpace<Sphere, int, 10, 10, 10>;

   TINYTEST_TIME_START()

      for (int i = 0; i < NUMBEROFTESTS; ++i)
      {
         ms_type ms;

         if (i % 10 == 0) std::cout << "*";

         for (int j = 0; j < i; ++j)
         {
            TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j));
            std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(j, Sphere(Vector3(0.2 + 0.0001 * i, 0.6 + 0.0002 * i, 0.4), 0.001))));
            TINYTEST_EQUAL(true, ir.second);
            for (int k = 0; k <= j; ++k) TINYTEST_ASSERT(ms.end() != ms.find(k));
            TINYTEST_ASSERT(ms.end() == ms.find(j + 1));
            TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j + 1));
         }
         for (int j = i; --j >= 0;)
         {
            for (int k = 0; k <= j; ++k)  TINYTEST_ASSERT(ms.end() != ms.find(k));
            TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j + 1));
            TINYTEST_ASSERT(ms.erase(j));
            TINYTEST_ASSERT(ms.end() == ms.find(j));
            TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j));
         }
      }

   std::cout << std::endl;

   for (int i = 0; i < NUMBEROFTESTS; ++i)
   {
      ms_type ms;

      if (i % 10 == 0) std::cout << "*";

      for (int j = 0; j < i; ++j)
      {
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j));
         std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(j, Sphere(Vector3(0.2 + 0.0001 * i, 0.6 + 0.0002 * i, 0.4), 0.001))));
         TINYTEST_EQUAL(true, ir.second);
         for (int k = 0; k <= j; ++k) TINYTEST_ASSERT(ms.end() != ms.find(k));
         TINYTEST_ASSERT(ms.end() == ms.find(j + 1));
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j + 1));
      }
      for (int j = 0; j < i; ++j)
      {
         for (int k = j; k < i; ++k) TINYTEST_ASSERT(ms.end() != ms.find(k));
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(i - j));
         TINYTEST_ASSERT(ms.erase(j));
         TINYTEST_ASSERT(ms.end() == ms.find(j));
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(i - j - 1));
      }
   }

   std::cout << std::endl;

   for (int i = 0; i < NUMBEROFTESTS; ++i)
   {
      ms_type ms;
      std::vector<int> id_list;

      if (i % 10 == 0) std::cout << "*";

      for (int j = 0; j < i; ++j)
      {
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j));
         std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(j, Sphere(Vector3(0.2 + 0.0001 * i, 0.6 + 0.0002 * i, 0.4), 0.001))));
         TINYTEST_EQUAL(true, ir.second);
         for (int k = 0; k <= j; ++k) TINYTEST_ASSERT(ms.end() != ms.find(k));
         TINYTEST_ASSERT(ms.end() == ms.find(j + 1));
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(j + 1));
         id_list.emplace_back(j);
      }

      RandomNumberGenerator rng;
      shuffle(rng, id_list);

      for (int j = 0; j < i; ++j)
      {
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(i - j));
         for (int k = 0; k < j; ++k) TINYTEST_ASSERT(ms.end() == ms.find(id_list[k]));
         for (int k = j; k < i; ++k)  TINYTEST_ASSERT(ms.end() != ms.find(id_list[k]));
         TINYTEST_ASSERT(ms.erase(id_list[j]));
         TINYTEST_ASSERT(ms.size() == static_cast<size_t>(i - j - 1));
      }
   }

   TINYTEST_TIME_STOP("Erase");
   std::cout << std::endl;
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceEachNeighbor()
{
   using ms_type = MatrixSpace<Sphere, int, 10, 10, 10>;

   ms_type ms;
   TINYTEST_ALMOST_EQUAL(0.1, ms.cell_size(), 0.001);

   ms.update(std::make_pair(0, Sphere(Vector3(0.2, 0.6, 0.4), 0.05)));
   TINYTEST_ASSERT(ms.end() != ms.find(0));
   TINYTEST_ASSERT(ms.end() == ms.find(1));
   ms.update(std::make_pair(1, Sphere(Vector3(0.2, 0.7, 0.5), 0.05)));
   TINYTEST_ASSERT(ms.end() != ms.find(0));
   TINYTEST_ASSERT(ms.end() != ms.find(1));
   TINYTEST_ASSERT(ms.end() == ms.find(2));
   ms.update(std::make_pair(2, Sphere(Vector3(0.9, 0.1, 0.4), 0.05)));
   TINYTEST_ASSERT(ms.end() != ms.find(0));
   TINYTEST_ASSERT(ms.end() != ms.find(1));
   TINYTEST_ASSERT(ms.end() != ms.find(2));
   TINYTEST_ASSERT(ms.end() == ms.find(3));
   ms.update(std::make_pair(3, Sphere(Vector3(0.9, 0.95, 0.4), 0.05)));
   TINYTEST_ASSERT(ms.end() != ms.find(0));
   TINYTEST_ASSERT(ms.end() != ms.find(1));
   TINYTEST_ASSERT(ms.end() != ms.find(2));
   TINYTEST_ASSERT(ms.end() != ms.find(3));
   TINYTEST_ASSERT(ms.end() == ms.find(4));

   {
      collector_set<ms_type> col;
      ms.each_neighbor(ms.index(Vector3(0.2, 0.6, 0.4)), col);
      TINYTEST_EQUAL(col.result.size(), 2);
      TINYTEST_ASSERT(col.result.find(0) != col.result.end());
      TINYTEST_ASSERT(col.result.find(1) != col.result.end());
   }
   {
      collector_set<ms_type> col;
      ms.each_neighbor(ms.index(Vector3(0.0, 0.1, 0.4)), col);
      TINYTEST_EQUAL(col.result.size(), 0);
   }
   {
      collector_set<ms_type> col;
      ms.each_neighbor_cyclic(ms.index(Vector3(0.0, 0.1, 0.4)), col);
      TINYTEST_EQUAL(col.result.size(), 1);
      TINYTEST_ASSERT(col.result.find(2) != col.result.end());
   }
   {
      collector_set<ms_type> col;
      ms.each_neighbor_cyclic(ms.index(Vector3(0.9, 0.0, 0.4)), col);
      TINYTEST_EQUAL(col.result.size(), 2);
      TINYTEST_ASSERT(col.result.find(2) != col.result.end());
      TINYTEST_ASSERT(col.result.find(3) != col.result.end());
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceReupdate()
{
   using ms_type = MatrixSpace<Sphere, int, 10, 10, 10>;

   ms_type ms;
   TINYTEST_ALMOST_EQUAL(0.1, ms.cell_size(), 0.001);

   TINYTEST_ASSERT(ms.update(std::make_pair(0, Sphere(Vector3(0.2, 0.6, 0.4), 0.05))).second);
   TINYTEST_ASSERT(ms.end() != ms.find(0));
   TINYTEST_ASSERT(ms.end() == ms.find(1));

   {
      collector_set<ms_type> col;
      ms.each_neighbor(ms.index(Vector3(0.2, 0.6, 0.4)), col);
      TINYTEST_EQUAL(col.result.size(), 1);
      TINYTEST_ASSERT(col.result.find(0) != col.result.end());
      TINYTEST_ASSERT(col.result.find(1) == col.result.end());
   }

   TINYTEST_ASSERT(!ms.update(std::make_pair(0, Sphere(Vector3(0.13, 0.83, 0.43), 0.05))).second);
   TINYTEST_ASSERT(ms.end() != ms.find(0));
   TINYTEST_ASSERT(ms.end() == ms.find(1));

   {
      collector_set<ms_type> col;
      ms.each_neighbor(ms.index(Vector3(0.25, 0.62, 0.43)), col);
      TINYTEST_EQUAL(col.result.size(), 0);
      TINYTEST_ASSERT(col.result.find(0) == col.result.end());
      TINYTEST_ASSERT(col.result.find(1) == col.result.end());
   }
   {
      collector_set<ms_type> col;
      ms.each_neighbor(ms.index(Vector3(0.23, 0.84, 0.45)), col);
      TINYTEST_EQUAL(col.result.size(), 1);
      TINYTEST_ASSERT(col.result.find(0) != col.result.end());
      TINYTEST_ASSERT(col.result.find(1) == col.result.end());
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceEachNeighbor2()
{
   using ms_type = MatrixSpace<Sphere, int, 10, 10, 10>;

   TINYTEST_TIME_START()

      for (double r = 0.01; r < 0.1; r += 0.01)
      {
         std::cout << "*";

         for (double o = 0.0; o < 0.9; o += STEPSIZE)
         {
            ms_type ms;
            TINYTEST_ALMOST_EQUAL(0.1, ms.cell_size(), 0.001);

            Vector3 centre(o, o, o);
            for (int i = 0; i < NUMBEROFTESTS / 5; ++i)
            {
               double t1(M_PI * 2 * (i % 10) / 10.), t2(M_PI * 2 * (i / 10) / 10.);
               const double _x = cos(t1) * r;
               const Vector3 p(centre.X() + _x * cos(t2), centre.Y() + sin(t1) * r, centre.Z() + _x * sin(t2));
               ms.update(std::make_pair(i, Sphere(p, r)));

               collector_count<ms_type> col;
               ms.each_neighbor(ms.index(centre), col);
               TINYTEST_EQUAL(col.count, i + 1);
            }
         }
      }
   std::cout << std::endl;

   for (double r = 0.01; r < 0.1; r += 0.01)
   {
      std::cout << "*";

      for (double o = 0.0; o < 0.9; o += STEPSIZE)
      {
         ms_type ms;
         TINYTEST_ALMOST_EQUAL(0.1, ms.cell_size(), 0.001);

         Vector3 centre(o, o, o);

         for (int i = 0; i < NUMBEROFTESTS / 5; ++i)
         {
            double t1(M_PI * 2 * (i % 10) / 10.), t2(M_PI * 2 * (i / 10) / 10.);
            const double _x = cos(t1) * r;
            const Vector3 p(centre.X() + _x * cos(t2), centre.Y() + sin(t1) * r, centre.Z() + _x * sin(t2));
            ms.update(std::make_pair(i, Sphere(p, r)));

            collector_count<ms_type> col;
            ms.each_neighbor_cyclic(ms.index(centre), col);
            TINYTEST_EQUAL(col.count, i + 1);
         }
      }
   }

   TINYTEST_TIME_STOP("Neighbor2");
   std::cout << std::endl;
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceTiming()
{
   // one liner time macro is a bit cryptic
   // TINYTEST_TIME_ELAPSE("Test", *"*/  Disk d();    /*"*);

   using ms_type = MatrixSpace<Sphere, int, 10, 10, 10>;
   std::cout << std::endl;

   TINYTEST_TIME_START()
      ms_type ms;
   TINYTEST_TIME_STOP("Create");

   ms_type ms;
   TINYTEST_TIME_START()
      ms.update(std::make_pair(0, Sphere(Vector3(0.2, 0.6, 0.4), 0.05)));
   TINYTEST_TIME_STOP("Update1");

   TINYTEST_TIME_START()
      ms.update(std::make_pair(1, Sphere(Vector3(0.3, 0.9, 0.2), 0.05)));
   TINYTEST_TIME_STOP("Update2");

   TINYTEST_TIME_START()
      for (int i = 0; i < 100; i++) ms.update(std::make_pair(i, Sphere(Vector3(0.3, 0.9, 0.2), 0.05)));
   TINYTEST_TIME_STOP("Update100");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceCylinderInsert()
{
   using ms_type = MatrixSpace<Cylinder, int, 10, 10, 10>;
   ms_type ms;
   ms.initialize(100.0);
   TINYTEST_ALMOST_EQUAL(100., ms.cell_size(), 0.001);

   {
      std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(0, Cylinder(Vector3(500, 500, 500), 25, Vector3(0, 0, 1), 1000))));
      TINYTEST_EQUAL(true, ir.second);  // Normal insert.
      TINYTEST_ASSERT(ms.end() != ms.find(0)); // Key 0 exists.
      TINYTEST_ASSERT(ms.end() == ms.find(1)); // Key 1 doesn't exist.
   }
   {
      // Update.
      std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(0, Cylinder(Vector3(500, 500, 500), 25, Vector3(0, 0, 1), 1000))));
      TINYTEST_EQUAL(false, ir.second); // False: this was an update.
      // ir.first is an iterator to the value you inserted. So accessing 
      // it's second element should return you the object.
      TINYTEST_EQUAL(Cylinder(Vector3(500, 500, 500), 25, Vector3(0, 0, 1), 1000), (*ir.first).second);
      TINYTEST_ASSERT(ms.end() != ms.find(0));
      TINYTEST_ASSERT(ms.end() == ms.find(1));
   }
   {
      // Another update.
      std::pair<ms_type::iterator, bool> ir(ms.update(std::make_pair(0, Cylinder(Vector3(500, 500, 500), 25, Vector3(0, 0, 1), 1000))));
      TINYTEST_EQUAL(false, ir.second);
      TINYTEST_EQUAL(Cylinder(Vector3(500, 500, 500), 25, Vector3(0, 0, 1), 1000), (*ir.first).second);
      TINYTEST_ASSERT(ms.end() != ms.find(0));
      TINYTEST_ASSERT(ms.end() == ms.find(1));
   }
   return 1;
}

int testMatrixSpaceCylinderEachNeighbor()
{
   using ms_type = MatrixSpace<Cylinder, int, 10, 10, 10>;
   ms_type ms;
   ms.initialize(100.0);
   TINYTEST_ALMOST_EQUAL(100., ms.cell_size(), 0.001);

   // Insert value 0.
   ms.update(std::make_pair(0, Cylinder(Vector3(500, 500, 0), 25, Vector3(0, 0, 1), 50)));
   TINYTEST_ASSERT(ms.end() != ms.find(0));
   TINYTEST_ASSERT(ms.end() == ms.find(1));

   {
      collector_set<ms_type> col;
      // Should return value 0.
      ms.each_neighbor(ms.index(Vector3(500, 500, 100)), col);
      TINYTEST_EQUAL(col.result.size(), 1);
      TINYTEST_ASSERT(col.result.find(0) != col.result.end());
   }

   {
      collector_set<ms_type> col;
      // No periodic boundary condition. Should return no values.
      // Behavior is unspecified for values at the boundary or out of the MatrixSpace (x,y,z >= 1000).
      ms.each_neighbor(ms.index(Vector3(500, 500, 900)), col);
      TINYTEST_EQUAL(col.result.size(), 0);
   }

   {
      collector_set<ms_type> col2;
      // Periodic boundary condition. Should return element 0 after applying 
      // periodic boundary condition in z (add 1000 to z coordinate of the 
      // origin of the cylinder to be in the same neighborhood as reference point), so: (0,0,1000).
      ms.each_neighbor_cyclic(ms.index(Vector3(500, 500, 900)), col2);
      TINYTEST_EQUAL(col2.result.size(), 1);
      TINYTEST_ASSERT(col2.result.find(0) != col2.result.end());
   }

   // Insert value 1.
   ms.update(std::make_pair(1, Cylinder(Vector3(500, 500, 900), 25, Vector3(0, 0, 1), 50)));
   {
      TINYTEST_ASSERT(ms.end() != ms.find(0));
      TINYTEST_ASSERT(ms.end() != ms.find(1));
      TINYTEST_ASSERT(ms.end() == ms.find(2));
   }

   {
      collector_set<ms_type> col2;
      // Periodic boundary condition. Should return element 0 (0, 0, 0) and element 1 (0,0,-1000).
      ms.each_neighbor_cyclic(ms.index(Vector3(500, 500, 0)), col2);
      TINYTEST_EQUAL(col2.result.size(), 2);
      TINYTEST_ASSERT(col2.result.find(0) != col2.result.end());
      TINYTEST_ASSERT(col2.result.find(1) != col2.result.end());
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSpaceParticles()
{
   using particle_matrix_type = MatrixSpace<Particle, ParticleID, 16, 16, 16>;
   particle_matrix_type pmat;
   pmat.initialize(100.0);
   TINYTEST_ALMOST_EQUAL(100., pmat.cell_size(), 1E-15);
   TINYTEST_ALMOST_EQUAL(1600., pmat.world_size().X(), 1E-15);
   TINYTEST_ALMOST_EQUAL(1600., pmat.world_size().Y(), 1E-15);
   TINYTEST_ALMOST_EQUAL(1600., pmat.world_size().Z(), 1E-15);

   auto pi_pair = std::make_pair<ParticleID, Particle>(ParticleID(100), Particle());
   pmat.update(pi_pair);

   // find particle
   particle_matrix_type::iterator i(pmat.find(pi_pair.first));
   TINYTEST_ASSERT(pmat.end() != i);

   pmat.erase(ParticleID(100));
   TINYTEST_ASSERT(pmat.size() == 0);

   i = pmat.find(pi_pair.first);
   TINYTEST_ASSERT(pmat.end() == i);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testMatrixSmall()
{

   using ms_type = MatrixSpace<Sphere, int, 1, 1, 1>;

   ms_type ms;
   Vector3 centre(0.5, 0.5, 0.5);
   ms.update(std::make_pair(1, Sphere(centre, 0.1)));

   collector_count<ms_type> col1;
   ms.each_neighbor(ms.index(centre), col1);
   TINYTEST_EQUAL(1, col1.count);

   collector_count<ms_type> col2;
   ms.each_neighbor_cyclic(ms.index(centre), col2);
   TINYTEST_EQUAL(3*3*3, col2.count);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(MatrixSpace);
TINYTEST_ADD_TEST(testMatrixSpaceSized);
TINYTEST_ADD_TEST(testMatrixSpaceUpdate);
TINYTEST_ADD_TEST(testMatrixSpaceErase);
TINYTEST_ADD_TEST(testMatrixSpaceEachNeighbor);
TINYTEST_ADD_TEST(testMatrixSpaceReupdate);
TINYTEST_ADD_TEST(testMatrixSpaceEachNeighbor2);
TINYTEST_ADD_TEST(testMatrixSpaceTiming);
TINYTEST_ADD_TEST(testMatrixSpaceCylinderInsert);
TINYTEST_ADD_TEST(testMatrixSpaceCylinderEachNeighbor);
TINYTEST_ADD_TEST(testMatrixSpaceParticles);
TINYTEST_ADD_TEST(testMatrixSmall);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
