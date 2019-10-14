#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "ReactionRuleID.hpp"
#include "ParticleID.hpp"
#include "Logger.hpp"
#include "Model.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class reaction_recorder
{
public:
   explicit reaction_recorder() : next_(nullptr) { }
   virtual ~reaction_recorder() { }

   virtual const char * type_name() const { return "ReactionRecorder"; };     // type implementation name (for Persistence)

   void StoreDecayReaction(double time, ReactionRuleID rid, ParticleID r)
   {
      StoreReaction(time, rid, r, ParticleID(0), ParticleID(0), ParticleID(0));
   }

   void StoreUniMolecularReaction(double time, ReactionRuleID rid, ParticleID r, ParticleID p)
   {
      StoreReaction(time, rid, r, ParticleID(0), p, ParticleID(0));
   }

   void StoreAnnihilationReaction(double time, ReactionRuleID rid, ParticleID r1, ParticleID r2)
   {
      StoreReaction(time, rid, r1, r2, ParticleID(0), ParticleID(0));
   }

   void StoreBindingReaction(double time, ReactionRuleID rid, ParticleID r1, ParticleID r2, ParticleID p)
   {
      StoreReaction(time, rid, r1, r2, p, ParticleID(0));
   }

   void StoreUnbindingReaction(double time, ReactionRuleID rid, ParticleID r, ParticleID p1, ParticleID p2)
   {
      StoreReaction(time, rid, r, ParticleID(0), p1, p2);
   }

   void link(reaction_recorder* next) { next_ = next; }
   reaction_recorder* next() const { return next_; }

protected:
   virtual void StoreReaction(double time, ReactionRuleID rid, ParticleID r1, ParticleID r2, ParticleID p1, ParticleID p2)
   {
      if (next_) next_->StoreReaction(time, rid, r1, r2, p1, p2);
   }

private:
   reaction_recorder* next_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class reaction_recorder_log : public reaction_recorder
{
public:
   explicit reaction_recorder_log(const ReactionRuleCollection& reaction_rules, const Model& model) : reaction_rules_(reaction_rules), model_(model), log_(Log("ReactionRecord")), first_(true)
   {
      log_.set_stream(std::clog);
      log_.set_flags(Logger::logflags::None);
   }

   const char * type_name() const override { return "ReactionRecorderLog"; };     // type implementation name (for Persistence)

   void set_output(std::ostream& stream) const { log_.set_stream(stream); }

   void StoreReaction(double time, ReactionRuleID rid, ParticleID r1, ParticleID r2, ParticleID p1, ParticleID p2) override
   {
      if (first_) print_header();
      log_.info() << std::scientific << std::setprecision(16) << std::setw(2 * c) << std::setfill(' ') <<
         time << std::setw(c / 2) << static_cast<idtype>(rid) << std::setw(c) << f(r1) << std::setw(c) << f(r2) << std::setw(c) << f(p1) << std::setw(c) << f(p2);
      log_.stream().flush(); // Force line to be written immediately, so output files can be monitored live during run
      reaction_recorder::StoreReaction(time, rid, r1, r2, p1, p2);
   }

private:
   const int c = 12;    // column width

   std::string f(ParticleID p)
   {
      std::stringstream s;
      idtype id = static_cast<idtype>(p);
      if (id != 0) s << id;
      return s.str();
   }

   void print_header()
   {
      log_.info() << std::setw(c / 2) << std::setfill(' ') << "Rule" << std::setw(2*c) << "k" << std::setw(2 * c) << "Reactant(s)" << std::setw(2 * c) << "Product(s)";

      // build rules table
      std::map<ReactionRuleID, std::string> table;
      for (const auto& r : reaction_rules_.get_reaction_rules())
      {
         const ReactionRule::reactants& ra = r.first;
         const ReactionRuleCollection::reaction_rule_set& rules = r.second;

         for (const auto &rule : rules)
         {
            if (rule.getK() == 0) continue;  // skip repulsive rules

            std::stringstream io;
            io << std::setw(c / 2) << static_cast<idtype>(rule.id()) << std::setw(2*c)  << std::scientific << std::setprecision(8) << rule.getK();
            std::stringstream tmp1;
            if (ra.size() > 0) tmp1 << model_.get_species_type_by_id(ra.item1()).name();
            if (ra.size() > 1) tmp1 << " + " << model_.get_species_type_by_id(ra.item2()).name();
            io << std::setw(2 * c) << tmp1.str();

            std::stringstream tmp2;
            const auto& products = rule.get_products();
            if (products.size() > 0) tmp2 << model_.get_species_type_by_id(products[0]).name();
            if (products.size() > 1) tmp2 << " + " << model_.get_species_type_by_id(products[1]).name();
            io << std::setw(2 * c) << tmp2.str();
            table.emplace(rule.id(), io.str());
         }
      }

      // show rules in order of id
      for (auto& line : table)
         log_.info() << line.second;

      // show header
      first_ = false;
      log_.info() << "\n\n" << std::setw(2 * c) << "Time [s]" << std::setw(c / 2) << "Rule" << std::setw(c) << "oldP1" << std::setw(c) << "oldP2" << std::setw(c) << "newP1" << std::setw(c) << "newP2";
      log_.stream().flush(); // Force line to be written immediately, so output files can be monitored live during run
   }

protected:
   friend class Persistence;

   const ReactionRuleCollection& reaction_rules_;
   const Model& model_;

   Logger& log_;
   bool first_;
};

// --------------------------------------------------------------------------------------------------------------------------------
