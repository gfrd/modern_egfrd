#ifndef EVENT_HPP
#define EVENT_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "DomainID.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class CustomAction
{
public:
   explicit CustomAction(double interval) noexcept : interval_(interval) {}
   virtual ~CustomAction() = default;

   double interval() const { return interval_; }
   virtual void do_action(double) {}

   virtual const char * type_name() const = 0;     // type implementation name (CopyNumber, Progress, etc)

private:
   const double interval_;
};

// --------------------------------------------------------------------------------------------------------------------------------

struct Event
{
   explicit Event() noexcept : time_(0), action_(actionType::Unspecified), u_(nullptr) {} // default constructor (use only for Persistence)

   explicit Event(double time, const DomainID did) noexcept : time_(time), action_(actionType::DomainUpdate), u_(did) {}

   explicit Event(double time, CustomAction* action) noexcept : time_(time), action_(actionType::CustomAction), u_(action) {}

   explicit Event(double time) noexcept : time_(time), action_(actionType::Unspecified), u_(nullptr) {}

   enum class actionType { Unspecified, DomainUpdate, CustomAction, };

   double time() const { return time_; }
   actionType action() const { return action_; }

   DomainID id() const { return u_.domainID_; }
   CustomAction* custom_action() const { return u_.custom_action_; }

   struct comparator
   {
      bool operator()(const Event& lhs, const Event& rhs) const
      {
         return lhs.time() <= rhs.time();
      }
   };

private:
   friend class Persistence;

   double time_;
   actionType action_;

   union uu       // unite domainID and CustomAction since there mutually exclusive, saves memory and time
   {
      uu(DomainID did) : domainID_(did) {}
      uu(CustomAction* ca) : custom_action_(ca) {}

      DomainID domainID_;
      CustomAction* custom_action_;
   } u_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif // EVENT_HPP