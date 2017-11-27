#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "DynamicPriorityQueue.hpp"
#include "EventID.hpp"
#include "Event.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class EventScheduler
{
public:

   EventScheduler() : queue_(), time_(0.0) {}

   double time() const { return time_; }

   size_t size() const { return queue_.size(); }

   bool empty() const { return queue_.empty(); }

   std::pair<EventID, const Event&> top() const { return queue_.top(); }

   std::pair<EventID, Event> pop()
   {
      if (queue_.empty()) throw std::out_of_range("queue is empty");
      std::pair<EventID, Event> top(queue_.top());      // value copy event (no ref), because pop throws it out of the Queue
      queue_.pop();
      time_ = top.second.time();
      return top;
   }

   std::pair<EventID, const Event&> second() const { return queue_.second(); }

   Event get(const EventID id) const { return queue_.get(id); }

   EventID add(const Event& event) { return queue_.push(event); }

   void remove(const EventID id) { queue_.pop(id); }

   void update(const EventID id, Event value) { queue_.replace(id, value); }

   void clear() { time_ = 0.0; queue_.clear(); }

   bool check() const { return queue_.check(); }

   bool has(EventID id) const { return queue_.has(id); }

private:
   friend class Persistence;

   DynamicPriorityQueue<EventID, Event> queue_;
   double time_;
};

// --------------------------------------------------------------------------------------------------------------------------------
