#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <sstream>
#include <Logger.hpp>
#include "World.hpp"
#include "Event.hpp"
#include "ReactionRecorder.hpp"
#include "ReactionRuleCollection.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class CopyNumbersInst : public CustomAction
{
   // Computes the instantaneous number of particles per speciestype at the interval
   // Output to std::cout (default) or given ostream.

public:
   explicit CopyNumbersInst(const World& world, double interval = 1e-3) : CustomAction(interval), world_(world), log_(Log("CopyNummber"))
   {
      log_.set_stream(std::cout);
      log_.set_flags(Logger::logflags::None);
   }

   void do_action(double time) override
   {
      if (time == 0.0) print_header();
      print_row(time);
   }

   const char* type_name() const override { return "CopyNumbersInst"; }

   void set_output(std::ostream& stream) const { log_.set_stream(stream); }

private:
   void print_header()
   {
      std::stringstream line;
       max_name_length_ = 1;

       // Determine the column size based on the maximal species name length
       for (auto& s : world_.get_species())
       {
           auto size = s.second.name().length() + 1;
           if (size > max_name_length_)
           {
               max_name_length_ = size;
           }
       }

       line << std::setw(12) << "Time";
      for (auto& s : world_.get_species())
         line << std::setw(12) << s.second.name();
      log_.info() << line.str();
       log_.stream().flush(); // Force line to be written immediately, so output files can be monitored live during run
   }

   void print_row(double time) const
   {
      std::stringstream line;
      line << std::scientific << std::setprecision(6) << std::setw(12) << time << std::fixed;
      for (auto& s : world_.get_species())
         line << std::setw(max_name_length_) << world_.get_particle_ids(s.first).size();
      log_.info() << line.str();
       log_.stream().flush(); // Force line to be written immediately, so output files can be monitored live during run
   }

protected:
   const World& world_;
   Logger& log_;
   int max_name_length_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class CopyNumbersAvg : public CustomAction, public reaction_recorder
{
   // Computes the average number of particles per speciestype with the interval
   // Output to std::cout (default) or given ostream.

public:
   explicit CopyNumbersAvg(const World& world, const ReactionRuleCollection& reaction_rules, double interval = 1e-3) :
      CustomAction(interval), world_(world), reaction_rules_(reaction_rules), log_(Log("CopyNummber"))
   {
      log_.set_stream(std::cout);
      log_.set_flags(Logger::logflags::None);
   }

   void do_action(double time) override
   {
      if (time == 0.0) print_header();
      print_row(time);
   }

   const char* type_name() const override { return "CopyNumbersAvg"; }

   void set_output(std::ostream& stream) const { log_.set_stream(stream); }

protected:
   void StoreReaction(double time, ReactionRuleID rid, ParticleID r1, ParticleID r2, ParticleID p1, ParticleID p2) override
   {
      reaction_recorder::StoreReaction(time, rid, r1, r2, p1, p2);

      // get involved species
      const auto& rr = reaction_rules_.get_reaction_rule(rid);
      std::set<SpeciesTypeID> sids;
      auto ra = rr.get_reactants();
      sids.insert(ra.item1());
      if (ra.size() > 1) sids.insert(rr.get_reactants().item2());
      for (auto sid : rr.get_products())
         sids.insert(sid);

      // update change in count
      for (auto sid : sids)
      {
         size_t count = world_.get_particle_ids(sid).size();
         avg_map_[sid].update_count(time, count);
      }
   }

   struct avg_accu
   {
      avg_accu() : time_(0), count_(0), accu_time_(0), accu_count_(0) {}

      double time_;
      size_t count_;
      double accu_time_;
      double accu_count_;

      void reset(size_t count)
      {
         count_ = count;
         time_ = 0.0;
         accu_time_ = accu_count_ = 0.0;
      }

      void update_count(double time, size_t count)
      {
         double delta = time - time_;
         accu_time_ += delta;
         accu_count_ += count_ * delta;
         count_ = count;
         time_ = time;
      }

      double update_time(double time)
      {
         if (time == 0) return static_cast<double>(count_);
         update_count(time, count_);
         double avg = accu_count_ / accu_time_;
         accu_time_ = accu_count_ = 0.0;
         return avg;
      }
   };

private:

   void print_header()
   {
      std::stringstream line;
      max_name_length_ = 1;

      // Determine the column size based on the maximal species name length
      for (auto& s : world_.get_species())
      {
          auto size = s.second.name().length() + 1;
          if (size > max_name_length_)
          {
              max_name_length_ = size;
          }
      }

       line << std::setw(12) << "Time";

       for (auto& s : world_.get_species())
       {
         line << std::setw(max_name_length_) << s.second.name();

         size_t count = world_.get_particle_ids(s.first).size();
         avg_map_[s.first].reset(count);
      }
      log_.info() << line.str();
      log_.stream().flush(); // Force line to be written immediately, so output files can be monitored live during run
   }

   void print_row(double time)
   {
      std::stringstream line;
      line << std::scientific << std::setprecision(6) << std::setw(12) << time << std::fixed;
      for (auto& s : world_.get_species())
         line << std::setw(max_name_length_) << std::fixed << std::setprecision(3) << avg_map_[s.first].update_time(time);
      log_.info() << line.str();
      log_.stream().flush(); // Force line to be written immediately, so output files can be monitored live during run
   }

protected:
   std::unordered_map<SpeciesTypeID, avg_accu> avg_map_;
   const World& world_;
   const ReactionRuleCollection& reaction_rules_;
   Logger& log_;
   int max_name_length_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class ParticlePositions : public CustomAction
{
   // Computes the positions of all particles at the interval
   // Beware of the performance penalty, for particle positions it is required to burst the domains!
   // Output to std::cout (default) or given ostream.

public:
   explicit ParticlePositions(const std::unique_ptr<EGFRDSimulator>& sim, double interval = 1e-3) : CustomAction(interval), sim_(sim), log_(Log("ParticlePositions"))
   {
      log_.set_stream(std::cout);
      log_.set_flags(Logger::logflags::None);
   }

   void do_action(double time) override
   {
      print(time);
   }

   const char* type_name() const override { return "ParticlePositions"; }

   void set_output(std::ostream& stream) const { log_.set_stream(stream); }

private:

   void print(double time) const
   {
      sim_->burst_all();
      auto& world = sim_->world();

      log_.info() << "ParticlePositions T=" << std::scientific << std::setprecision(6) << std::setw(12) << time << " N=" << std::fixed << world.num_particles();
      log_.info() << std::setw(8) << "PID" << std::setw(8) << "SID" << std::setw(16) << "X" << std::setw(16) << "Y" << std::setw(16) << "Z";
      
      for (auto& pip : world.get_particles())
         log_.info() << std::fixed << std::setw(8) << pip.first() << std::setw(8) << pip.second.sid()() << 
                        std::scientific << std::setprecision(6) << std::setw(16) << pip.second.position().X() << std::setw(16) << pip.second.position().Y() << std::setw(16) << pip.second.position().Z();
      log_.info();
      log_.stream().flush(); // Force line to be written immediately, so output files can be monitored live during run
   }

protected:
   const std::unique_ptr<EGFRDSimulator>& sim_;
   Logger& log_;
};

// --------------------------------------------------------------------------------------------------------------------------------
