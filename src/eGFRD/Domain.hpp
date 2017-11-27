#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "DomainID.hpp"
#include "ShellID.hpp"
#include "EventScheduler.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Shell;

// --------------------------------------------------------------------------------------------------------------------------------

class Domain
{
public:
   using shell_id_pair_ref = std::pair<ShellID, std::reference_wrapper<const Shell>>;

   enum class Multiplicity
   {
       SINGLE = 1,
       PAIR = 2,
       MULTI = 3,
   };

   enum class EventType
   {
      INIT,
      SINGLE_REACTION,
      SINGLE_ESCAPE,
      SINGLE_INTERACTION,
      COM_ESCAPE,
      IV_EVENT,
      IV_ESCAPE,
      IV_REACTION,
      IV_INTERACTION,
      BURST,
      MULTI_ESCAPE,
      MULTI_UNIMOLECULAR_REACTION,
      MULTI_BIMOLECULAR_REACTION,
      MULTI_DIFFUSION,
   };


   explicit Domain(const DomainID id) noexcept : domainID_(id), eventID_(0), eventType_(EventType::INIT), last_time_(0.), dt_(0.) {}

   virtual ~Domain() = default;

   DomainID id() const { return domainID_; }

   EventID eventID() const { return eventID_; }

   void set_event_id(const EventID eid) { eventID_ = eid; }

   EventType eventType() const { return eventType_; }

   double last_time() const { return last_time_; }

   double dt() const { return dt_; }

   virtual size_t num_shells() const = 0;

   virtual Multiplicity multiplicity() const = 0;          // 1=single, 2=pair, 3=multi

   virtual std::string as_string() const = 0;      // runtime parameters (particle pos, radius, shells, etc)

   virtual const char * type_name() const = 0;     // type implementation name (SingleSpherical, Multi, PairCyclindrical, etc)

   virtual shell_id_pair_ref get_shell() const = 0;     // return the shell when num_shells == 1, otherwise throw exception

   virtual std::vector<shell_id_pair_ref> get_shell_list() const = 0;     // return all shells in this domain

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_last_time(double t) { last_time_ = t; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_burst_time(double time)
   {
      // burst now, update dt_
      dt_ = time - last_time_;
      eventType_ = EventType::BURST;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:
   friend class Persistence;

   DomainID domainID_;
   EventID eventID_;
   EventType eventType_;
   double last_time_;
   double dt_;
};


// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Domain& d)
{
   stream << d.type_name() << "{" << d.id() << ", " << d.eventID() << ", last=" << d.last_time() << ", dt=" << d.dt() << d.as_string() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------
