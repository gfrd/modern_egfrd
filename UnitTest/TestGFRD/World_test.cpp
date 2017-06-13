#include "World.hpp"
#include "World_test.hpp"
#include "Model.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldCreate()
{
   World w;

   Model m;                                                    // we need a model, because it generates SpeciesTypeID numbers and the default (world) StructureTypeID
   StructureTypeID sid = m.get_def_structure_type_id();
   w.add_structure_type(m.get_structure_type_by_id(sid));
   w.set_def_structure_type_id(sid);

   SpeciesType s1("s1", sid, .3, 0.05);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .2, 0.03);
   m.add_species_type(s2);
   SpeciesType s3("s3", sid, .1, 0.02);
   m.add_species_type(s3);

   TINYTEST_EQUAL(w.num_species(), 0);
   w.add_species_type(s1);
   TINYTEST_EQUAL(w.num_species(), 1);
   w.add_species_type(s2);
   TINYTEST_EQUAL(w.num_species(), 2);
   w.add_species_type(s3);
   TINYTEST_EQUAL(w.num_species(), 3);

   // static type checking
    //w.has_particle(sid);       // pass speciesTypeId as ParticleID

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

#if defined(_MSC_VER)
__declspec(noinline)
#endif
int TestWorldParameter(const World& world, uint count)
{
   TINYTEST_EQUAL(world.num_species(), count);
   return 1;
}

#if defined(_MSC_VER)
__declspec(noinline)
#endif
int TestWorldMoveParameter(World&& world, uint count)
{
   TINYTEST_EQUAL(world.num_species(), count);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldContainerCopyConstruct()
{
   World w;

   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();
   w.add_structure_type(m.get_structure_type_by_id(sid));
   w.set_def_structure_type_id(sid);

   SpeciesType s1("s1", sid, .3, 0.05);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .2, 0.03);
   m.add_species_type(s2);
   SpeciesType s3("s3", sid, .1, 0.02);
   m.add_species_type(s3);


   //World w2 = w;         // cannot value copy a World (are you nuts?)

   TestWorldParameter(w, 0);        // pass world as ref, not as value-type!
   w.add_species_type(s1);
   TestWorldParameter(w, 1);
   w.add_species_type(s2);
   TestWorldParameter(w, 2);
   w.add_species_type(s3);
   TestWorldParameter(w, 3);

   TestWorldMoveParameter(std::move(w), 3);        // NOTE you can move a world around, but beware that its implementation may be still be a value copy (probably due to const members in some structures)!

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldNewParticles()
{
   World w;

   Model m;                                                    // we need a model, because it generates SpeciesTypeID numbers and the default (world) StructureTypeID
   StructureTypeID sid = m.get_def_structure_type_id();
   w.add_structure_type(m.get_structure_type_by_id(sid));
   w.set_def_structure_type_id(sid);

   Vector3 x = w.world_size() / 2;
   auto region = std::make_shared<CuboidalRegion>(CuboidalRegion("world", sid, StructureID(0), Box(x, x)));
   w.set_def_structure(region);
   StructureID wid = w.get_def_structure_id();

   SpeciesType s1("s1", sid, .3, 0.05);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .2, 0.06);
   m.add_species_type(s2);

   for (auto s : m.get_species())
      w.add_species_type(s.second);

   auto p1(w.add_particle(s1.id(), wid, Vector3(.2, .2, .2)));
   auto p2(w.add_particle(s2.id(), wid, Vector3(.29, .27, .28)));
   TINYTEST_ASSERT(p2.first != p1.first);

   auto check1 = w.check_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check1.size() == 0);
   auto check2 = w.check_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(check2.size() == 0);
   auto checkt1 = w.test_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(!checkt1);
   auto checkt2 = w.test_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(!checkt2);


   double length12 = (p1.second.position() - p2.second.position()).length();
   double r12 = p1.second.shape().radius() + p2.second.shape().radius();
   TINYTEST_ASSERT(length12 > r12);



   //TINYTEST_ASSERT(!boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p1.second.shape(), array_gen(p1.first))));
   //TINYTEST_ASSERT(!boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p2.second.shape(), array_gen(p2.first))));
   //TINYTEST_ASSERT(!w.check_particle_overlap(p1.second.shape(), p1.first));
   //TINYTEST_ASSERT(!w.check_particle_overlap(p2.second.shape(), p2.first));

   //auto p3(w.add_particle(s1.id(), wid, Vector3(.35, .32, .34)));
   auto p3(w.new_particle(Particle(s1, wid, Vector3(.35, .32, .34))));        // alternative construction
   TINYTEST_ASSERT(p3.first != p1.first);
   TINYTEST_ASSERT(p3.first != p2.first);

   auto check3 = w.check_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check3.size() == 0);
   auto check4 = w.check_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(check4.size() != 0);
   auto checkt3 = w.test_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(!checkt3);
   auto checkt4 = w.test_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(checkt4);

   double length13 = (p1.second.position() - p3.second.position()).length();
   double r13 = p1.second.shape().radius() + p3.second.shape().radius();
   TINYTEST_ASSERT(length13 > r13);
   double length23 = (p2.second.position() - p3.second.position()).length();
   double r23 = p2.second.shape().radius() + p3.second.shape().radius();
   TINYTEST_ASSERT(length23 < r23);


   //TINYTEST_ASSERT(boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p2.second.shape(), array_gen(p2.first))));
   //TINYTEST_ASSERT(boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p2.second.shape(), array_gen(p1.first, p2.first))));
   //TINYTEST_ASSERT(!boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p2.second.shape(), array_gen(p2.first, p3.first))));
   //TINYTEST_ASSERT(boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p3.second.shape(), array_gen(p3.first))));
   //TINYTEST_ASSERT(!boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p3.second.shape(), array_gen(p2.first, p3.first))));
   //TINYTEST_ASSERT(boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p3.second.shape(), array_gen(p1.first, p3.first))));

   //TINYTEST_ASSERT(!w.check_particle_overlap(p1.second.shape(), p1.first));
   //TINYTEST_ASSERT(w.check_particle_overlap(p2.second.shape(), p2.first));
   //TINYTEST_ASSERT(boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p2.second.shape(), array_gen(p1.first, p2.first))));
   //TINYTEST_ASSERT(!boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p2.second.shape(), array_gen(p3.first, p2.first))));
   //TINYTEST_ASSERT(w.check_particle_overlap(p3.second.shape(), p3.first));
   //TINYTEST_ASSERT(!boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p3.second.shape(), array_gen(p2.first, p3.first))));
   //TINYTEST_ASSERT(boost::scoped_ptr<particle_id_pair_and_distance_list>(w.check_particle_overlap(p3.second.shape(), array_gen(p1.first, p3.first))));*/


   auto range = w.get_particles();

   //auto b = w.particles_begin();
   //auto e = w.particles_end();

   //TINYTEST_EQUAL(std::distance(b, e), std::distance(range.begin(), range.end()));


   //   TINYTEST_ASSERT(std::is_sorted(b, e));
   //   TINYTEST_ASSERT(std::is_sorted(range.begin(), range.end()));

   for (auto pip : range)
   {
      TINYTEST_ASSERT(pip.first);
   }


   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldOverlap()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();

   SpeciesType s1("s1", sid, .3, 0.5);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .2, 0.6);
   m.add_species_type(s2);

   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   TINYTEST_ALMOST_EQUALS(10.0, w.world_size().X());
   TINYTEST_ALMOST_EQUALS(10.0 / CompileConfigSimulator::MatrixCellsX, w.cell_size());
   TINYTEST_EQUAL(CompileConfigSimulator::MatrixCellsX, w.matrix_size()[0]);

   auto p1(w.add_particle(s1.id(), wid, Vector3(4.0, 5.0, 5.0)));
   auto p2(w.add_particle(s1.id(), wid, Vector3(6.01, 5.0, 5.0)));

   auto check1 = w.check_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check1.size() == 0);
   auto check2 = w.check_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(check2.size() == 0);

   auto p3(w.add_particle(s2.id(), wid, Vector3(5.0, 5.0, 5.0)));

   auto check3 = w.check_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check3.size() != 0);
   TINYTEST_EQUAL(1, check3.size());
   TINYTEST_EQUAL(p3.first, (check3)[0].first.first);
   auto distance13 = (p1.second.position() - p3.second.position()).length() - (p1.second.radius() + p3.second.radius());        // distance takes the radius of the other particle into account, but not the radius of test particle
   TINYTEST_ALMOST_EQUALS(distance13, (check3)[0].second);

   auto check4 = w.check_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(check4.size() != 0);
   TINYTEST_EQUAL(1, check4.size());
   TINYTEST_EQUAL(p3.first, (check4)[0].first.first);
   auto distance23 = (p2.second.position() - p3.second.position()).length() - (p2.second.radius() + p3.second.radius());
   TINYTEST_ALMOST_EQUALS(distance23, (check4)[0].second);

   auto check5 = w.check_particle_overlap(p3.second.shape(), p3.first);
   TINYTEST_ASSERT(check5.size() != 0);
   TINYTEST_EQUAL(2, check5.size());
   TINYTEST_EQUAL(p1.first, (check5)[0].first.first);
   TINYTEST_EQUAL(p2.first, (check5)[1].first.first);
   auto distance31 = (p3.second.position() - p1.second.position()).length() - (p3.second.radius() + p1.second.radius());
   auto distance32 = (p3.second.position() - p2.second.position()).length() - (p3.second.radius() + p2.second.radius());
   TINYTEST_ALMOST_EQUALS(distance31, (check5)[0].second);
   TINYTEST_ALMOST_EQUALS(distance32, (check5)[1].second);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldOverlapCyclic()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();

   SpeciesType s1("s1", sid, .3, 0.5);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .2, 0.6);
   m.add_species_type(s2);

   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   TINYTEST_ALMOST_EQUALS(10.0, w.world_size().X());
   TINYTEST_ALMOST_EQUALS(10.0 / CompileConfigSimulator::MatrixCellsX, w.cell_size());
   TINYTEST_EQUAL(CompileConfigSimulator::MatrixCellsX, w.matrix_size()[0]);

   auto p1(w.add_particle(s1.id(), wid, Vector3(0.0, 5.0, 5.0)));
   auto p2(w.add_particle(s1.id(), wid, Vector3(8.0, 5.0, 5.0)));

   auto check1 = w.check_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check1.size() == 0);
   auto check2 = w.check_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(check2.size() == 0);

   auto p3(w.add_particle(s2.id(), wid, Vector3(9.0, 5.0, 5.0)));

   auto check3 = w.check_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check3.size() != 0);
   TINYTEST_EQUAL(1, check3.size());
   TINYTEST_EQUAL(p3.first, (check3)[0].first.first);
   auto distance13 = w.world_size().X() - (p1.second.position() - p3.second.position()).length() - (p1.second.radius() + p3.second.radius());
   TINYTEST_ALMOST_EQUALS(distance13, (check3)[0].second);

   auto check4 = w.check_particle_overlap(p2.second.shape(), p2.first);
   TINYTEST_ASSERT(check4.size() != 0);
   TINYTEST_EQUAL(1, check4.size());
   TINYTEST_EQUAL(p3.first, (check4)[0].first.first);
   auto distance23 = (p2.second.position() - p3.second.position()).length() - (p2.second.radius() + p3.second.radius());
   TINYTEST_ALMOST_EQUALS(distance23, (check4)[0].second);

   auto check5 = w.check_particle_overlap(p3.second.shape(), p3.first);
   TINYTEST_ASSERT(check5.size() != 0);
   TINYTEST_EQUAL(2, check5.size());
   TINYTEST_EQUAL(p2.first, (check5)[0].first.first);
   TINYTEST_EQUAL(p1.first, (check5)[1].first.first);
   auto distance32 = (p3.second.position() - p2.second.position()).length() - (p3.second.radius() + p2.second.radius());
   auto distance31 = w.world_size().X() - (p3.second.position() - p1.second.position()).length() - (p3.second.radius() + p1.second.radius());
   TINYTEST_ALMOST_EQUALS(distance32, (check5)[0].second);
   TINYTEST_ALMOST_EQUALS(distance31, (check5)[1].second);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldTransaction1()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();

   SpeciesType s1("s1", sid, .3, 0.5);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .2, 0.6);
   m.add_species_type(s2);

   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   TINYTEST_ALMOST_EQUALS(10.0, w.world_size().X());
   TINYTEST_ALMOST_EQUALS(10.0 / CompileConfigSimulator::MatrixCellsX, w.cell_size());
   TINYTEST_EQUAL(CompileConfigSimulator::MatrixCellsX, w.matrix_size()[0]);

   std::pair<ParticleID, Particle> p1, p2;
   {
      auto tx(w.create_transaction());

      p1 = tx->new_particle(Particle(s1, wid, Vector3(.2, .2, .2)));
      p2 = tx->new_particle(Particle(s2, wid, Vector3(.29, .27, .28)));

      TINYTEST_ASSERT(tx->get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(tx->get_particle(p2.first).second == p2.second);
      TINYTEST_ASSERT(w.get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(w.get_particle(p2.first).second == p2.second);

      auto ap = tx->get_added_particles();
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p1.first) != ap.end());
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p2.first) != ap.end());

      auto mp = tx->get_modified_particles();
      TINYTEST_ASSERT(std::distance(mp.begin(), mp.end()) == 0);

      auto rp = tx->get_removed_particles();
      TINYTEST_ASSERT(std::distance(rp.begin(), rp.end()) == 0);

      tx->rollback();
   }

   TINYTEST_CHECK_THROW(w.get_particle(p1.first), not_found);
   TINYTEST_CHECK_THROW(w.get_particle(p2.first), not_found);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldTransaction2()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();

   SpeciesType s1("s1", sid, .3, 0.5);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .3, 0.8);
   m.add_species_type(s2);

   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   Vector3 pos1 = Vector3(.2, .2, .2);
   Vector3 pos2 = Vector3(.3, .3, .3);

   auto p1(w.new_particle(Particle(s1, wid, pos1)));
   auto p2(w.new_particle(Particle(s2, wid, pos2)));

   {
      auto tx(w.create_transaction());

      TINYTEST_ASSERT(tx->get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(tx->get_particle(p2.first).second == p2.second);
      TINYTEST_ASSERT(w.get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(w.get_particle(p2.first).second == p2.second);

      {
         p1.second.position() += Vector3(5.0, 5.0, 5.0);
         tx->update_particle(p1);

         p2.second.position() += Vector3(2.0, 2.0, 2.0);
         tx->update_particle(p2);
      }

      TINYTEST_ASSERT(w.get_particle(p1.first).second.position() == Vector3(5.2, 5.2, 5.2));
      TINYTEST_ASSERT(w.get_particle(p2.first).second.position() == Vector3(2.3, 2.3, 2.3));

      auto ap = tx->get_added_particles();
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p1.first) == ap.end());
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p2.first) == ap.end());
      auto rp = tx->get_removed_particles();
      TINYTEST_ASSERT(std::find(rp.begin(), rp.end(), p1.first) == rp.end());
      TINYTEST_ASSERT(std::find(rp.begin(), rp.end(), p2.first) == rp.end());
      auto mp = tx->get_modified_particles();
      TINYTEST_ASSERT(std::find(mp.begin(), mp.end(), p1.first) != mp.end());
      TINYTEST_ASSERT(std::find(mp.begin(), mp.end(), p2.first) != mp.end());
      tx->rollback();
   }

   TINYTEST_ASSERT(w.get_particle(p1.first).second.position() == pos1);
   TINYTEST_ASSERT(w.get_particle(p2.first).second.position() == pos2);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldTransaction3()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();

   SpeciesType s1("s1", sid, .3, 0.5);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .3, 0.8);
   m.add_species_type(s2);

   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   auto p1(w.new_particle(Particle(s1, wid, Vector3(.2, .2, .2))));
   auto p2(w.new_particle(Particle(s2, wid, Vector3(.29, .27, .28))));
   World::particle_id_pair p3;

   {
      auto tx(w.create_transaction());

      p3 = tx->new_particle(Particle(s1, wid, Vector3(.4, .2, .1)));

      TINYTEST_ASSERT(tx->get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(tx->get_particle(p2.first).second == p2.second);
      TINYTEST_ASSERT(tx->get_particle(p3.first).second == p3.second);
      TINYTEST_ASSERT(w.get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(w.get_particle(p2.first).second == p2.second);
      TINYTEST_ASSERT(w.get_particle(p3.first).second == p3.second);

      tx->remove_particle(p1.first);
      TINYTEST_CHECK_THROW(tx->get_particle(p1.first), not_found);
      TINYTEST_ASSERT(tx->get_particle(p2.first).second == p2.second);
      TINYTEST_CHECK_THROW(w.get_particle(p1.first), not_found);

      auto ap = tx->get_added_particles();
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p1.first) == ap.end());
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p2.first) == ap.end());
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p3.first) != ap.end());
      auto rp = tx->get_removed_particles();
      TINYTEST_ASSERT(std::find(rp.begin(), rp.end(), p1.first) != rp.end());
      TINYTEST_ASSERT(std::find(rp.begin(), rp.end(), p2.first) == rp.end());
      TINYTEST_ASSERT(std::find(rp.begin(), rp.end(), p3.first) == rp.end());
      auto mp = tx->get_modified_particles();
      TINYTEST_ASSERT(std::find(mp.begin(), mp.end(), p1.first) == mp.end());
      TINYTEST_ASSERT(std::find(mp.begin(), mp.end(), p2.first) == mp.end());
      TINYTEST_ASSERT(std::find(mp.begin(), mp.end(), p3.first) == mp.end());
      tx->rollback();
   }

   TINYTEST_ASSERT(w.get_particle(p1.first).second == p1.second);
   TINYTEST_ASSERT(w.get_particle(p2.first).second == p2.second);
   TINYTEST_CHECK_THROW(w.get_particle(p3.first), not_found);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldTransaction4()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();

   SpeciesType s1("s1", sid, .3, 0.5);
   m.add_species_type(s1);
   SpeciesType s2("s2", sid, .3, 0.8);
   m.add_species_type(s2);

   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   auto p1(w.new_particle(Particle(s1, wid, Vector3(.2, .2, .2))));
   auto p2(w.new_particle(Particle(s2, wid, Vector3(.29, .27, .28))));

   {
      auto tx(w.create_transaction());

      TINYTEST_ASSERT(tx->get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(tx->get_particle(p2.first).second == p2.second);
      TINYTEST_ASSERT(w.get_particle(p1.first).second == p1.second);
      TINYTEST_ASSERT(w.get_particle(p2.first).second == p2.second);

      {
         auto new_p(p2.second);
         new_p.position() = Vector3(.35, .30, .38);
         tx->update_particle(World::particle_id_pair(p1.first, new_p));  // change p1
      }

      tx->remove_particle(p1.first);                                          // remove p1
      TINYTEST_CHECK_THROW(tx->get_particle(p1.first), not_found);
      TINYTEST_ASSERT(tx->get_particle(p2.first).second == p2.second);
      TINYTEST_CHECK_THROW(w.get_particle(p1.first), not_found);

      auto ap = tx->get_added_particles();
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p1.first) == ap.end());
      TINYTEST_ASSERT(std::find(ap.begin(), ap.end(), p2.first) == ap.end());
      auto rp = tx->get_removed_particles();
      TINYTEST_ASSERT(std::find(rp.begin(), rp.end(), p1.first) != rp.end());
      TINYTEST_ASSERT(std::find(rp.begin(), rp.end(), p2.first) == rp.end());
      auto mp = tx->get_modified_particles();
      TINYTEST_ASSERT(std::find(mp.begin(), mp.end(), p1.first) == mp.end());
      TINYTEST_ASSERT(std::find(mp.begin(), mp.end(), p2.first) == mp.end());
      tx->rollback();
   }

   TINYTEST_ASSERT(w.get_particle(p1.first).second == p1.second);
   TINYTEST_ASSERT(w.get_particle(p2.first).second == p2.second);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldThrowParticles()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();

   auto s1 = m.add_species_type(SpeciesType("A", sid, .3, 0.5));

   World w;
   w.initialize(10, m);

   TINYTEST_ALMOST_EQUALS(10.0, w.world_size().X());
   TINYTEST_ALMOST_EQUALS(10.0 / CompileConfigSimulator::MatrixCellsX, w.cell_size());
   TINYTEST_EQUAL(CompileConfigSimulator::MatrixCellsX, w.matrix_size()[0]);

   RandomNumberGenerator rng;
   w.throwInParticles(s1, 20, rng);

   return 1;
}


// --------------------------------------------------------------------------------------------------------------------------------

int testWorldThrowParticles2()
{
   Model m;
   auto s2 = m.add_structure_type(StructureType("plane"));
   auto s1 = m.add_species_type(SpeciesType("A", s2, .3, 0.05));
   auto sw = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), .3, 0.025));

   World w;
   w.initialize(10, m);
   auto wsid = w.get_def_structure_id();

   auto plane = std::make_shared<PlanarSurface>(PlanarSurface("plane", s2, wsid, Plane(Vector3(5, 5, 5), Vector3::ux, Vector3::uy, 5, 5, false)));
   w.add_structure(plane);

   //auto disk = std::make_shared<DiskSurface>(DiskSurface("disk", s2, wsid, Disk(Vector3(5,5,5), 1.0, Vector3::uz)));
   //w.add_structure(disk);

   //auto sphere = std::make_shared<SphericalSurface>(SphericalSurface("sphere", s2, wsid, Sphere(Vector3(5, 5, 5), 1.0)));
   //w.add_structure(sphere);

   //auto cyl = std::make_shared<CylindricalSurface>(CylindricalSurface("cylinder", s2, wsid, Cylinder(Vector3(5, 5, 5), 1.0, Vector3::uy, 2.0)));
   //w.add_structure(cyl);


   RandomNumberGenerator rng;
   w.throwInParticles(s1, 200, rng, false, Vector3(1, 1, 1), Vector3(9, 9, 9));
   w.throwInParticles(sw, 200, rng, false);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldOverlapSphere()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();
   auto s1 = m.add_species_type(SpeciesType("s1", sid, .3, 1.0));
   
   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   auto p1(w.add_particle(s1, wid, Vector3(5.0, 5.0, 5.0)));
   auto p2(w.add_particle(s1, wid, Vector3(7.01, 5.0, 5.0)));

   auto check1a = w.check_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check1a.size() == 0);
   auto check1b = w.test_particle_overlap(p1.second.shape(), p1.first);
   TINYTEST_ASSERT(check1b == false);

   Sphere sphere(Vector3(5, 5, 5), 5);
   auto check2 = w.check_particle_overlap(sphere);
   TINYTEST_ASSERT(check2.size() == 2);
   TINYTEST_ASSERT(check2[0].first.first == p1.first);
   TINYTEST_ALMOST_EQUALS(-6.0, check2[0].second);
   TINYTEST_ASSERT(check2[1].first.first == p2.first);
   TINYTEST_ALMOST_EQUALS(-3.99, check2[1].second);         // (5+5-7.01)-(1)

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testWorldOverlapCylinder()
{
   Model m;
   StructureTypeID sid = m.get_def_structure_type_id();
   auto s1 = m.add_species_type(SpeciesType("s1", sid, .3, 1.0));
   //auto s2 = m.add_species_type(SpeciesType("s2", sid, .2, 0.1));

   World w;
   w.initialize(10, m);
   StructureID wid = w.get_def_structure_id();

   auto p1(w.add_particle(s1, wid, Vector3(5.0, 5.0, 5.0)));
   auto p2(w.add_particle(s1, wid, Vector3(7.01, 5.0, 5.0)));

   Cylinder cyl(Vector3(5, 5, 5), 5, Vector3::uy, 2.5);
   auto check2 = w.check_particle_overlap(cyl);
   TINYTEST_ASSERT(check2.size() == 2);
   TINYTEST_ASSERT(check2[0].first.first == p1.first);
   TINYTEST_ALMOST_EQUALS(-3.5, check2[0].second);          // (-2.5)-(1)
   TINYTEST_ASSERT(check2[1].first.first == p2.first);
   TINYTEST_ALMOST_EQUALS(-3.5, check2[1].second);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------


TINYTEST_START_SUITE(World);
TINYTEST_ADD_TEST(testWorldCreate);
TINYTEST_ADD_TEST(testWorldNewParticles);
TINYTEST_ADD_TEST(testWorldOverlap);
TINYTEST_ADD_TEST(testWorldOverlapCyclic);
TINYTEST_ADD_TEST(testWorldTransaction1);
TINYTEST_ADD_TEST(testWorldTransaction2);
TINYTEST_ADD_TEST(testWorldTransaction3);
TINYTEST_ADD_TEST(testWorldTransaction4);
TINYTEST_ADD_TEST(testWorldThrowParticles);
TINYTEST_ADD_TEST(testWorldThrowParticles2);
TINYTEST_ADD_TEST(testWorldContainerCopyConstruct);
TINYTEST_ADD_TEST(testWorldOverlapSphere);
TINYTEST_ADD_TEST(testWorldOverlapCylinder);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
