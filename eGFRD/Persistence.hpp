#ifndef PERSISTENCE_HPP
#define PERSISTENCE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <string>
#include <fstream>
#include <map>
#include <set>
#include <memory>

// --------------------------------------------------------------------------------------------------------------------------------

// forward declares
class PairGreensFunction;
class GreensFunction;
class CustomAction;
class World;
struct Vector3;
template<typename TObj, typename TKey, uint NX, uint NY, uint NZ> class MatrixSpace;
class StructureContainer;
class RandomNumberGenerator;
class ReactionRuleCollection;
class EGFRDSimulator;
class EventScheduler;

// --------------------------------------------------------------------------------------------------------------------------------

class Persistence
{
   const char signature[4] = { 'G','F','R','D' };
   const int16_t version_number = 0x0101;

public:
   using callback_fn = CustomAction* (*)(std::string name, double interval);

   Persistence() noexcept : if_(), of_(), version_(0), callback_(nullptr) {}
   ~Persistence() { if (of_.is_open()) of_.close(); if (if_.is_open()) if_.close(); }

   GFRD_EXPORT void store(std::string filename);
   GFRD_EXPORT bool retreive(std::string filename, callback_fn callback = nullptr);

   GFRD_EXPORT void store_egfrd(const EGFRDSimulator& gfrd);
   GFRD_EXPORT void retreive_egfrd(EGFRDSimulator& gfrd);

   GFRD_EXPORT void store_world(const World& w);
   GFRD_EXPORT void retreive_world(World& w);

   GFRD_EXPORT void store_rng(const RandomNumberGenerator& rng);
   GFRD_EXPORT void retreive_rng(RandomNumberGenerator& rng);

   GFRD_EXPORT void store_reactionrules(const ReactionRuleCollection& rng);
   GFRD_EXPORT void retreive_reactionrules(ReactionRuleCollection& rng);

protected:

   template<typename TObj, typename TKey, uint NX, uint NY, uint NZ>
   void store_matrixspace(const MatrixSpace<TObj, TKey, NX, NY, NZ>& w);
   template<typename TObj, typename TKey, uint NX, uint NY, uint NZ>
   void retreive_matrixspace(MatrixSpace<TObj, TKey, NX, NY, NZ>& w);

   void store_structurecontainer(const StructureContainer& v);
   void retreive_structurecontainer(StructureContainer& v);

   void store_eventscheduler(const EventScheduler& es);
   void retreive_eventscheduler(EventScheduler& es);

   void store_greensfunction(const std::unique_ptr<GreensFunction>& gf);
   void retreive_greensfunction(std::unique_ptr<GreensFunction>& gf);
   void store_greensfunction(const std::unique_ptr<PairGreensFunction>& gf);
   void retreive_greensfunction(std::unique_ptr<PairGreensFunction>& gf);

   template<typename TKey, typename TObj>
   void store_pool(const std::map<TKey, std::set<TObj>>& pool);
   template<typename TKey, typename TObj>
   void retreive_pool(std::map<TKey, std::set<TObj>>& pool);

   void store_vector(const Vector3& v);
   void retreive_vector(Vector3& v);

   void store_string(const std::string& v);
   void retreive_string(std::string& v);

private:

   void check(const char*chk, const char* message);
   
   template<typename T>
   void write(const T& t) { of_.write(reinterpret_cast<const char*>(&t), sizeof(T)); }
   
   template<typename T>
   T read() { T t; if_.read(reinterpret_cast<char*>(&t), sizeof(T)); return t; }

   template<typename T>
   void read(T& t) { if_.read(reinterpret_cast<char*>(&t), sizeof(T)); }

   std::ifstream if_;
   std::ofstream of_;
   int16_t version_;
   callback_fn callback_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* PERSISTENCE_HPP */
