// --------------------------------------------------------------------------------------------------------------------------------

#include <cstring>
#include <cctype>
#include "DefsEgfrd.hpp"
#include "Persistence.hpp"
#include "World.hpp"
#include "StructureContainer.hpp"
#include "ReactionRuleCollection.hpp"
#include "EventScheduler.hpp"
#include "EGFRDSimulator.hpp"
#include <GreensFunction3DAbsSym.hpp>
#include <GreensFunction3DRadAbs.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store(std::string filename)
{
   of_.open(filename, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
   of_.write(signature, sizeof(signature));
   write(version_number);
   int16_t platform = static_cast<int16_t>(sizeof(size_t));
   write(platform);
}

bool Persistence::retreive(std::string filename, callback_fn callback)
{
   // set CustomAction resolver callback
   callback_ = callback;

   if_.open(filename, std::ofstream::in | std::ofstream::binary);
   if (!if_.is_open()) return false;
   check("GFRD", "Expected GFRD signature.");
   read<int16_t>(version_);
   THROW_UNLESS_MSG(illegal_state, version_ <= version_number, "Cannot load state from future version.");
   int16_t platform = read<int16_t>();
   THROW_UNLESS_MSG(illegal_state, platform == static_cast<int16_t>(sizeof(size_t)), "Cannot load state from different platform (x86/x64).");
   return true;
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_egfrd(const EGFRDSimulator& gfrd)
{
   of_ << 'S' << 'I' << 'M';
   store_world(gfrd.world_);
   store_reactionrules(gfrd.reaction_rules_);
   store_rng(gfrd.rng_);

   write(gfrd.time_);
   write(gfrd.dt_);
   write(gfrd.num_steps_);
   write(gfrd.sidgen_.next_);
   write(gfrd.didgen_.next_);

   store_eventscheduler(gfrd.scheduler_);

   //domain_map;
   write(gfrd.domains_.size());
   for (const auto& d : gfrd.domains_)
   {
      write(d.first);
      std::string rtti = d.second->type_name();
      store_string(rtti);

      write(d.second->domainID_);
      write(d.second->eventID_);
      write(d.second->eventType_);
      write(d.second->last_time_);
      write(d.second->dt_);

      if (rtti == "SingleSpherical")
      {
         auto single = reinterpret_cast<Single*>(d.second.get());
         write(single->pid_pair_);
         write(single->sid_pair_);
         THROW_UNLESS_MSG(not_implemented, single->rrule_ == nullptr, "Rule drawn between states???");
         store_greensfunction(single->gf_);
      }
      else if (rtti == "PairSpherical")
      {
         auto pair = reinterpret_cast<PairSpherical*>(d.second.get());
         write(pair->pid_pair1_);
         write(pair->pid_pair2_);
         write(pair->sid_pair_);
         store_vector(pair->iv_);
         store_vector(pair->com_);
         write(pair->a_R_);
         write(pair->a_r_);
         write(pair->single1_ktotal_);
         write(pair->single2_ktotal_);
         write(pair->pid_reactingsingle_);
         THROW_UNLESS_MSG(not_implemented, pair->rrule_ == nullptr, "Rule drawn between states???");
         THROW_UNLESS_MSG(not_implemented, pair->gf_tmp_ == nullptr, "tmpGf set between states???");
         store_greensfunction(pair->gf_com_);
         store_greensfunction(pair->gf_iv_);
      }
      else if (rtti == "Multi")
      {
         auto multi = reinterpret_cast<Multi*>(d.second.get());
         size_t size = multi->particles_.size();
         write(size);
         of_.write(reinterpret_cast<const char*>(&multi->particles_[0]), size * sizeof(ParticleID));
         size = multi->shell_map_.size();
         write(size);
         for (const auto& m : multi->shell_map_)
         {
            write(m.first);
            write(m.second);
         }
         write(multi->reaction_length_);
         write(multi->start_time_);
      }
      else THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");

   }

   store_matrixspace(gfrd.shellmat_);
}

void Persistence::retreive_egfrd(EGFRDSimulator& gfrd)
{
   check("SIM", "Expected data for EGFRD Simulator.");
   retreive_world(gfrd.world_);
   retreive_reactionrules(gfrd.reaction_rules_);
   retreive_rng(gfrd.rng_);

   read(gfrd.time_);
   read(gfrd.dt_);
   read(gfrd.num_steps_);
   read(gfrd.sidgen_.next_);
   read(gfrd.didgen_.next_);

   retreive_eventscheduler(gfrd.scheduler_);
   //domain_map
   size_t size = read<size_t>();
   gfrd.domains_.clear();
   for (size_t i = 0; i < size; ++i)
   {
      DomainID key = read<DomainID>();
      std::string rtti;
      retreive_string(rtti);

      DomainID did = read<DomainID>();
      EventID eid = read<EventID>();
      Domain::EventType type = read<Domain::EventType>();
      double last_time = read<double>();
      double dt = read<double>();

      if (rtti == "SingleSpherical")
      {
         Single::particle_id_pair pid_pair;
         Single::shell_id_pair sid_pair;
         read(pid_pair);
         read(sid_pair);
         const auto& rr = gfrd.reaction_rules_.query_reaction_rules(pid_pair.second.sid());
         auto single = std::make_unique<SingleSpherical>(SingleSpherical(did, pid_pair, sid_pair, rr));
         single->eventID_ = eid;
         single->eventType_ = type;
         single->last_time_ = last_time;
         single->dt_ = dt;
         retreive_greensfunction(single->gf_);
         gfrd.domains_[key] = std::move(single);
      }
      else if (rtti == "PairSpherical")
      {
         Single::particle_id_pair pid1_pair, pid2_pair;
         Single::shell_id_pair sid_pair;
         read(pid1_pair);
         read(pid2_pair);
         read(sid_pair);
         const auto& rr = gfrd.reaction_rules_.query_reaction_rules(pid1_pair.second.sid(), pid2_pair.second.sid());
         auto pair = std::make_unique<PairSpherical>(PairSpherical(did, pid1_pair, pid2_pair, sid_pair, rr));
         pair->eventID_ = eid;
         pair->eventType_ = type;
         pair->last_time_ = last_time;
         pair->dt_ = dt;
         retreive_vector(pair->iv_);
         retreive_vector(pair->com_);
         read(pair->a_R_);
         read(pair->a_r_);
         read(pair->single1_ktotal_);
         read(pair->single2_ktotal_);
         read(pair->pid_reactingsingle_);
         retreive_greensfunction(pair->gf_com_);
         retreive_greensfunction(pair->gf_iv_);
         gfrd.domains_[key] = std::move(pair);
      }
      else if (rtti == "Multi")
      {
         auto multi = std::make_unique<Multi>(Multi(did, gfrd.world_, gfrd.shellmat_, gfrd.reaction_rules_, gfrd.rng_, gfrd.rrec_, gfrd.time_));
         multi->eventID_ = eid;
         multi->eventType_ = type;
         multi->last_time_ = last_time;
         multi->dt_ = dt;

         size_t size2 = read<size_t>();
         multi->particles_.resize(size2);
         if_.read(reinterpret_cast<char *>(&multi->particles_[0]), size2 * sizeof(ParticleID));
         size2 = read<size_t>();
         multi->shell_map_.clear();
         for (size_t j = 0; j < size2; ++j)
         {
            ParticleID pid = read<ParticleID>();
            ShellID value = read<ShellID>();
            multi->shell_map_[pid] = value;
         }
         read(multi->reaction_length_);
         read(multi->start_time_);
         gfrd.domains_[key] = std::move(multi);
      }
      else THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");
   }

   retreive_matrixspace(gfrd.shellmat_);
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_world(const World& w)
{
   of_ << 'W' << 'R' << 'L' << 'D';
   write(w.pidgen_.next_);
   store_matrixspace(w.pmat_);
   store_structurecontainer(w.structures_);
   //species_map
   write(w.species_map_.size());
   for (const auto& k : w.species_map_)
   {
      write(k.first);
      write(k.second.id_);
      store_string(k.second.name_);
      write(k.second.diffusion_coef_);
      write(k.second.drift_velocity_);
      write(k.second.radius_);
      write(k.second.structure_type_id_);
   }
   //structure_type_map
   write(w.structure_type_map_.size());
   for (const auto& k : w.structure_type_map_)
   {
      write(k.first);
      write(k.second.id_);
      store_string(k.second.name_);
   }

   store_pool(w.structure_pool_);
   store_pool(w.particle_pool_);
   store_pool(w.particleonstruct_pool_);

   write(w.default_structure_type_id_);
   write(w.default_structure_id_);
}

void Persistence::retreive_world(World& w)
{
   check("WRLD", "Expected data for World.");
   read(w.pidgen_.next_);
   retreive_matrixspace(w.pmat_);
   retreive_structurecontainer(w.structures_);
   //species_map
   size_t size = read<size_t>();
   w.species_map_.clear();
   for (size_t i = 0; i < size; ++i)
   {
      SpeciesTypeID key = read<SpeciesTypeID>();

      SpeciesTypeID id = read<SpeciesTypeID>();
      std::string name;
      retreive_string(name);

      double diffusion_coef = read<double>();
      double drift_velocity = read<double>();
      double radius = read<double>();
      StructureTypeID sid = read<StructureTypeID>();
      auto s = SpeciesType(name, sid, diffusion_coef, radius, drift_velocity);
      s.id_ = id;
      w.species_map_[key] = s;
   }
   //structure_type_map
   size = read<size_t>();
   w.structure_type_map_.clear();
   for (size_t i = 0; i < size; ++i)
   {
      StructureTypeID key = read<StructureTypeID>();

      StructureTypeID sid = read<StructureTypeID>();
      std::string name;
      retreive_string(name);
      auto s = StructureType(name);
      s.id_ = sid;
      w.structure_type_map_[key] = s;
   }

   retreive_pool(w.structure_pool_);
   retreive_pool(w.particle_pool_);
   retreive_pool(w.particleonstruct_pool_);

   read(w.default_structure_type_id_);
   read(w.default_structure_id_);
}

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TObj, typename TKey, uint NX, uint NY, uint NZ>
void Persistence::store_matrixspace(const MatrixSpace<TObj, TKey, NX, NY, NZ>& ms)
{
   of_ << 'M' << 'S';
   auto skey = sizeof(TKey);
   auto sobj = sizeof(TObj);
   write(skey);
   write(sobj);

   store_vector(ms.world_size_);
   write(ms.cell_size_);
   auto msize = ms.matrix_size();
   write(msize[0]);
   write(msize[1]);
   write(msize[2]);
   for (const auto &m : ms.matrix_)
   {
      size_t size = m.size();
      write(size);
      if (size) of_.write(reinterpret_cast<const char*>(&m[0]), size * sizeof(m[0]));
   }
   write(ms.rmap_.size());
   for (const auto& k : ms.rmap_)
   {
      write(k.first);
      write(k.second);
   }
   write(ms.values_.size());
   for (const auto& k : ms.values_)
   {
      write(k.first);
      write(k.second);
   }
}

template<typename TObj, typename TKey, uint NX, uint NY, uint NZ>
void Persistence::retreive_matrixspace(MatrixSpace<TObj, TKey, NX, NY, NZ>& ms)
{
   check("MS", "Expected data for MatrixSpace.");
   size_t skey = read<size_t>();
   size_t sobj = read<size_t>();
   THROW_UNLESS_MSG(illegal_state, skey == sizeof(TKey) && sobj == sizeof(TObj), "MatrixSpace types do not match.");

   retreive_vector(ms.world_size_);
   read(ms.cell_size_);
   typename MatrixSpace<TObj, TKey, NX, NY, NZ>::size_type nx, ny, nz;
   read(nx);
   read(ny);
   read(nz);
   auto msize = ms.matrix_size();
   THROW_UNLESS_MSG(illegal_state, nx == msize[0] && ny == msize[1] && nz == msize[2], "Cannot load matrix " << nx << "x" << ny << "x" << nz << " into matrix " << msize[0] << "x" << msize[1] << "x" << msize[2] << ".");
   for (auto &m : ms.matrix_)
   {
      size_t size = read<size_t>();
      m.resize(size);
      if (size) if_.read(reinterpret_cast<char*>(&m[0]), size * sizeof(typename MatrixSpace<TObj, TKey, NX, NY, NZ>::size_type));
   }
   size_t size = read<size_t>();
   ms.rmap_.clear();
   ms.rmap_.reserve(size);
   for (size_t i = 0; i < size; ++i)
   {
      TKey key = read<TKey>();
      typename MatrixSpace<TObj, TKey, NX, NY, NZ>::size_type value;
      read(value);
      ms.rmap_[key] = value;
   }
   size = read<size_t>();
   ms.values_.clear();
   ms.values_.reserve(size);
   for (size_t i = 0; i < size; ++i)
   {
      TKey id = read<TKey>();
      TObj value = read<TObj>();
      ms.values_.emplace_back(std::make_pair(id, value));
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_structurecontainer(const StructureContainer& s)
{
   of_ << 'S' << 'C';
   write(s.structidgen_.next_);
   write(s.structure_map_.size());
   for (const auto& k : s.structure_map_)
   {
      write(k.first);

      auto str = dynamic_cast<Structure*>(k.second.get());
      store_string(str->name_);
      write(str->sid_);
      write(str->id_);
      write(str->parent_id_);

      std::string rtti = k.second->type_name();
      store_string(rtti);
      if (rtti == "CuboidalRegion")
      {
         auto cr = dynamic_cast<CuboidalRegion*>(k.second.get());
         write(cr->shape());
      }
      else if (rtti == "CylindricalSurface")
      {
         auto cs = dynamic_cast<CylindricalSurface*>(k.second.get());
         write(cs->shape());
      }
      else if (rtti == "DiskSurface")
      {
         auto ds = dynamic_cast<DiskSurface*>(k.second.get());
         write(ds->shape());
      }
      else if (rtti == "PlanarSurface")
      {
         auto ps = dynamic_cast<PlanarSurface*>(k.second.get());
         write(ps->shape());
      }
      else if (rtti == "SphericalSurface")
      {
         auto ss = dynamic_cast<SphericalSurface*>(k.second.get());
         write(ss->shape());
      }
      else
         THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");
   }

   store_pool(s.structure_substructures_map_);

   write(s.planar_structs_bc_.neighbor_mapping_.size());
   for (const auto& k : s.planar_structs_bc_.neighbor_mapping_)
   {
      write(k.first);
      write(k.second.size());
      for (const auto& l : k.second)
      {
         write(l);
      }
   }
}

void Persistence::retreive_structurecontainer(StructureContainer& s)
{
   check("SC", "Expected data for StructureContainer.");
   read(s.structidgen_.next_);
   size_t size = read<size_t>();
   s.structure_map_.clear();
   for (size_t i = 0; i < size; ++i)
   {
      StructureID key = read<StructureID>();

      std::string name;
      retreive_string(name);
      StructureTypeID sid = read<StructureTypeID>();
      StructureID id = read<StructureID>();
      StructureID pid = read<StructureID>();

      std::string rtti;
      retreive_string(rtti);
      if (rtti == "CuboidalRegion")
      {
         Box b; read(b);
         auto cr = std::make_shared<CuboidalRegion>(CuboidalRegion(name, sid, pid, b));
         cr->id_ = id; s.structure_map_[key] = cr;
      }
      else if (rtti == "CylindricalSurface")
      {
         Cylinder c; read(c);
         auto cr = std::make_shared<CylindricalSurface>(CylindricalSurface(name, sid, pid, c));
         cr->id_ = id; s.structure_map_[key] = cr;
      }
      else if (rtti == "DiskSurface")
      {
         Disk d; read(d);
         auto cr = std::make_shared<DiskSurface>(DiskSurface(name, sid, pid, d));
         cr->id_ = id; s.structure_map_[key] = cr;
      }
      else if (rtti == "PlanarSurface")
      {
         Plane p; read(p);
         auto cr = std::make_shared<PlanarSurface>(PlanarSurface(name, sid, pid, p));
         cr->id_ = id; s.structure_map_[key] = cr;
      }
      else if (rtti == "SphericalSurface")
      {
         Sphere sp; read(sp);
         auto cr = std::make_shared<SphericalSurface>(SphericalSurface(name, sid, pid, sp));
         cr->id_ = id; s.structure_map_[key] = cr;
      }
      else
         THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");
   }

   retreive_pool(s.structure_substructures_map_);

   size = read<size_t>();
   s.planar_structs_bc_.neighbor_mapping_.clear();
   for (size_t i = 0; i < size; ++i)
   {
      StructureID key = read<StructureID>();

      size_t size2 = read<size_t>();
      THROW_UNLESS_MSG(illegal_state, size2 == 4, "Expected 4 data structures for planar surface boundary conditions.");

      ConnectivityContainer<4>::obj_data_pair_array_type tmp;
      if_.read(reinterpret_cast<char*>(&tmp[0]), sizeof(ConnectivityContainer<4>::obj_data_pair_array_type));
      s.planar_structs_bc_.neighbor_mapping_[key] = tmp;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_eventscheduler(const EventScheduler& es)
{
   of_ << 'E' << 'S';
   size_t size = es.queue_.items_.size();
   write(size);
   if (size) of_.write(reinterpret_cast<const char*>(&es.queue_.items_[0]), size * sizeof(es.queue_.items_[0]));

   size = es.queue_.heap_.size();
   write(size);
   if (size) of_.write(reinterpret_cast<const char*>(&es.queue_.heap_[0]), size * sizeof(es.queue_.heap_[0]));

   size = es.queue_.position_vector_.size();
   write(size);
   if (size) of_.write(reinterpret_cast<const char*>(&es.queue_.position_vector_[0]), size * sizeof(es.queue_.position_vector_[0]));

   write(es.queue_.index_map_.size());
   for (const auto& k : es.queue_.index_map_)
   {
      write(k.first);
      write(k.second);
   }

   write(es.queue_.idgen_.next_);
   write(es.time_);

   // fix for custom actions
   std::map<std::pair<std::string, double>, std::vector<EventID>> pool;
   for (const auto& e : es.queue_.items_)
   {
      if (e.second.action() == Event::actionType::CustomAction)
      {
         auto ca = e.second.custom_action();
         if (ca != nullptr)
         {
            std::string name = ca->type_name();
            double interval = ca->interval();
            pool[std::make_pair(name, interval)].emplace_back(e.first);
         }
      }
   }
   write(pool.size());
   for (const auto& k : pool)
   {
      store_string(k.first.first);
      write(k.first.second);

      size = k.second.size();
      write(size);
      if (size) of_.write(reinterpret_cast<const char*>(&k.second[0]), size * sizeof(k.second[0]));
   }
}

void Persistence::retreive_eventscheduler(EventScheduler& es)
{
   check("ES", "Expected data for EventScheduler");

   size_t size = read<size_t>();
   es.queue_.items_.resize(size);
   if (size) if_.read(reinterpret_cast<char*>(&es.queue_.items_[0]), size * sizeof(es.queue_.items_[0]));

   size = read<size_t>();
   es.queue_.heap_.resize(size);
   if (size) if_.read(reinterpret_cast<char*>(&es.queue_.heap_[0]), size * sizeof(es.queue_.heap_[0]));

   size = read<size_t>();
   es.queue_.position_vector_.resize(size);
   if (size) if_.read(reinterpret_cast<char*>(&es.queue_.position_vector_[0]), size * sizeof(es.queue_.position_vector_[0]));

   size = read<size_t>();
   es.queue_.index_map_.clear();
   es.queue_.index_map_.reserve(size);
   for (size_t i = 0; i < size; ++i)
   {
      EventID key = read<EventID>();
      size_t index = read<size_t>();
      es.queue_.index_map_[key] = index;
   }

   read(es.queue_.idgen_.next_);
   read(es.time_);

   // fix custom events (they have pointers, that are now pointing to wherever)
   size = read<size_t>();
   for (size_t i = 0; i < size; ++i)
   {
      std::string name;
      retreive_string(name);
      double interval = read<double>();

      CustomAction* ca = callback_ != nullptr ? callback_(name, interval) : nullptr;
      if (ca == nullptr) Log("Persistence").warn() << "CustomAction '"<< name << "' in EventScheduler cannot be restored. It will be ignored.";

      size_t size2 = read<size_t>();
      std::vector<EventID> tmp;
      tmp.resize(size2);
      if (size2)
      {
         if_.read(reinterpret_cast<char*>(&tmp[0]), size2 * sizeof(tmp[0]));
         for (auto& eid : tmp)
         {
            auto idx = es.queue_.findIndex(eid);
            es.queue_.items_[idx].second.u_.custom_action_ = ca;
         }
      }
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TKey, typename TObj>
void Persistence::store_pool(const std::map<TKey, std::set<TObj>>& pool)
{
   of_ << 'P';
   write(pool.size());
   for (const auto& k : pool)
   {
      write(k.first);
      write(k.second.size());
      for (auto l : k.second)
      {
         write(l);
      }
   }
}

template<typename TKey, typename TObj>
void Persistence::retreive_pool(std::map<TKey, std::set<TObj>>& pool)
{
   check("P", "Expected data for Pool.");
   size_t size = read<size_t>();
   pool.clear();
   for (size_t i = 0; i < size; ++i)
   {
      TKey key = read<TKey>();
      size_t size2 = read<size_t>();

      std::vector<TObj> tmp;
      tmp.resize(size2);
      if (size2) if_.read(reinterpret_cast<char*>(&tmp[0]), size2 * sizeof(TObj));
      std::set<TObj> value;
      for (auto s : tmp) value.emplace(s);
      pool[key] = value;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_vector(const Vector3& v)
{
   write(v.x_);
   write(v.y_);
   write(v.z_);
}

void Persistence::retreive_vector(Vector3& v)
{
   read(v.x_);
   read(v.y_);
   read(v.z_);
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_string(const std::string& text)
{
   of_ << 'S';
   uint16_t size = static_cast<uint16_t>(text.length());
   write(size);
   of_.write(reinterpret_cast<const char*>(&text[0]), size * sizeof(std::string::value_type));
}

void Persistence::retreive_string(std::string& text)
{
   check("S", "Expected data for String.");
   uint16_t size = read<uint16_t>();
   text.resize(size);
   if_.read(reinterpret_cast<char*>(&text[0]), sizeof(std::string::value_type) * size);
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_rng(const RandomNumberGenerator& rng)
{
   of_ << 'R' << 'N' << 'G';
   of_ << rng.rng_;
   of_ << 'G' << 'N' << 'R';
}

void Persistence::retreive_rng(RandomNumberGenerator& rng)
{
   check("RNG", "Expected data for RandomNumberGenerator");
   if_ >> rng.rng_;
   char dmy; do { if_ >> dmy; } while (dmy != 'G' && !if_.bad());
   check("NR", "Expected end for RandomNumberGenerator");
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_reactionrules(const ReactionRuleCollection& rrc)
{
   of_ << 'R' << 'R' << 'C';
   write(rrc.reaction_rules_map_.size());
   for (const auto& r : rrc.reaction_rules_map_)
   {
      write(r.first);
      write(r.second.size());
      for (const auto& rr : r.second)
      {
         write(rr.id_);
         write(rr.reactants_);
         size_t size = rr.products_.size();
         write(size);
         if (size) of_.write(reinterpret_cast<const char*>(&rr.products_[0]), size * sizeof(SpeciesTypeID));
         write(rr.k_);
      }
   }

   write(rrc.interaction_rules_map_.size());
   for (const auto& r : rrc.interaction_rules_map_)
   {
      write(r.first);
      write(r.second.size());
      for (const auto& ir : r.second)
      {
         write(ir.id_);
         write(ir.reactants_);
         size_t size = ir.products_.size();
         write(size);
         if (size) of_.write(reinterpret_cast<const char*>(&ir.products_[0]), size * sizeof(SpeciesTypeID));
         write(ir.k_);
      }
   }
   write(rrc.sgen_.next_);
}

void Persistence::retreive_reactionrules(ReactionRuleCollection& rrc)
{
   check("RRC", "Expected data for ReactionRuleCollection");
   size_t size = read<size_t>();
   rrc.reaction_rules_map_.clear();
   for (size_t i = 0; i < size; ++i)
   {
      ReactionRule::reactants key = read<ReactionRule::reactants>();
      size_t size2 = read<size_t>();

      ReactionRuleCollection::reaction_rule_set rrs;
      for (size_t j = 0; j < size2; ++j)
      {
         ReactionRule rr(SpeciesTypeID(1), 0.0, {});
         read(rr.id_);
         read(rr.reactants_);
         size_t size3 = read<size_t>();
         rr.products_.resize(size3);
         if (size3) if_.read(reinterpret_cast<char*>(&rr.products_[0]), size3 * sizeof(SpeciesTypeID));
         read(rr.k_);

         THROW_UNLESS_MSG(illegal_state, rr.reactants_ == key, "Illegal entry in ReactionRuleCollection");
         rrs.insert(rr);
      }
      rrc.reaction_rules_map_[key] = rrs;
   }
   size = read<size_t>();;
   rrc.interaction_rules_map_.clear();
   for (size_t i = 0; i < size; ++i)
   {
      InteractionRule::reactants key = read<InteractionRule::reactants>();
      size_t size2 = read<size_t>();

      ReactionRuleCollection::interaction_rule_set irs;
      for (size_t j = 0; j < size2; ++j)
      {
         InteractionRule ir(SpeciesTypeID(1), StructureTypeID(1), 0.0, {});
         read(ir.id_);
         read(ir.reactants_);
         size_t size3 = read<size_t>();
         ir.products_.resize(size3);
         if (size3) if_.read(reinterpret_cast<char*>(&ir.products_[0]), size3 * sizeof(SpeciesTypeID));
         read(ir.k_);

         THROW_UNLESS_MSG(illegal_state, ir.reactants_ == key, "Illegal entry in ReactionRuleCollection");
         irs.insert(ir);
      }
      rrc.interaction_rules_map_[key] = irs;
   }
   read(rrc.sgen_.next_);

}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::check(const char*chk, const char* message)
{
   size_t length = std::strlen(chk);
   ASSERT(length > 0 && length <= 8);
   char check[8];
   if_.read(&check[0], length);
   for (size_t i = 0; i < length; ++i)
      THROW_UNLESS_MSG(illegal_state, check[i] == chk[i], message);
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_greensfunction(const std::unique_ptr<GreensFunction>& gf)
{
   std::string rtti = gf != nullptr ? gf->type_name() : "NULL";
   store_string(rtti);
   if (gf == nullptr) return;

   write(gf->getD());

   if (rtti == "GreensFunction3DAbsSym")
   {
      auto gf2 = reinterpret_cast<GreensFunction3DAbsSym*>(gf.get());
      write(gf2->geta());
   }
   else
      THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");
}

void Persistence::retreive_greensfunction(std::unique_ptr<GreensFunction>& gf)
{
   std::string rtti;
   retreive_string(rtti);
   if (rtti == "NULL") { gf = nullptr; return; }

   double D = read<double>();

   if (rtti == "GreensFunction3DAbsSym")
   {
      double a = read<double>();
      gf = std::make_unique<GreensFunction3DAbsSym>(GreensFunction3DAbsSym(D, a));
   }
   else
      THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");
}

// --------------------------------------------------------------------------------------------------------------------------------

void Persistence::store_greensfunction(const std::unique_ptr<PairGreensFunction>& gf)
{
   std::string rtti = gf != nullptr ? gf->type_name() : "NULL";
   store_string(rtti);
   if (gf == nullptr) return;

   write(gf->getD());
   write(gf->getkf());
   write(gf->getr0());
   write(gf->getSigma());

   if (rtti == "GreensFunction3DRadAbs")
   {
      auto gf2 = reinterpret_cast<GreensFunction3DRadAbs*>(gf.get());
      write(gf2->geta());
   }
   else
      THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");
}

void Persistence::retreive_greensfunction(std::unique_ptr<PairGreensFunction>& gf)
{
   std::string rtti;
   retreive_string(rtti);
   if (rtti == "NULL") { gf = nullptr; return; }

   double D = read<double>();
   double kf = read<double>();
   double r0 = read<double>();
   double sigma = read<double>();

   if (rtti == "GreensFunction3DRadAbs")
   {
      double a = read<double>();
      gf = std::make_unique<GreensFunction3DRadAbs>(GreensFunction3DRadAbs(D, kf, r0, sigma, a));
   }
   else
      THROW_EXCEPTION(not_implemented, "Storage of type " << rtti << " is not implemented.");
}

// --------------------------------------------------------------------------------------------------------------------------------

