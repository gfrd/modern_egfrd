#include "EventScheduler_test.hpp"
#include "EventScheduler.hpp"
#include "randomNumberGenerator.hpp"
#include <algorithm>

// --------------------------------------------------------------------------------------------------------------------------------

int testEventSchedulerEmpty()
{
   EventScheduler scheduler;

   TINYTEST_ASSERT(scheduler.size() == 0);
   TINYTEST_EQUAL(scheduler.time(), 0.0);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testEventSchedulerOneEvent()
{
   EventScheduler scheduler;

   Event ev1(1.0);

   auto id = scheduler.add(ev1);
   TINYTEST_EQUAL(0.0, scheduler.time());
   TINYTEST_EQUAL(1.0, scheduler.top().second.time());
   TINYTEST_EQUAL(id, scheduler.top().first);

   auto data = scheduler.pop();
   TINYTEST_EQUAL(data.first, id);
   TINYTEST_EQUAL(data.second.time(), ev1.time());

   TINYTEST_EQUAL(0, scheduler.size());
   TINYTEST_EQUAL(1.0, scheduler.time());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testEventSchedulerSecondEvent()
{
   EventScheduler scheduler;

   Event eventa(1.0);
   Event eventb(0.5);

   auto id1 = scheduler.add(eventa);
   auto id2 = scheduler.add(eventb);
   TINYTEST_EQUAL(2, scheduler.size());

   auto sec = scheduler.second();
   TINYTEST_EQUAL(1.0, sec.second.time());
   TINYTEST_EQUAL(id1, sec.first);

   auto first = scheduler.top();
   TINYTEST_EQUAL(0.5, first.second.time());
   TINYTEST_EQUAL(id2, first.first);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testEventSchedulerReplaceEvent()
{
   EventScheduler scheduler;

   Event ev1(0.1);
   Event ev2(0.5);
   Event ev3(0.9);

   auto id1 = scheduler.add(ev1);
   auto id2 = scheduler.add(ev2);
   auto id3 = scheduler.add(ev3);
   TINYTEST_EQUAL(3, scheduler.size());

   scheduler.update(id3, Event(0.2));

   auto data1 = scheduler.pop();
   TINYTEST_EQUAL(data1.first, id1);
   auto data3 = scheduler.pop();
   TINYTEST_EQUAL(data3.first, id3);
   auto data2 = scheduler.pop();
   TINYTEST_EQUAL(data2.first, id2);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testEventSchedulerRemoveEvent()
{
   EventScheduler scheduler;

   Event ev1(0.1);
   Event ev2(0.5);
   Event ev3(0.9);

   auto id1 = scheduler.add(ev1);
   auto id2 = scheduler.add(ev2);
   auto id3 = scheduler.add(ev3);
   TINYTEST_EQUAL(3, scheduler.size());

   scheduler.remove(id2);
   TINYTEST_EQUAL(2, scheduler.size());

   auto data1 = scheduler.pop();
   TINYTEST_EQUAL(data1.first, id1);
   auto data3 = scheduler.pop();
   TINYTEST_EQUAL(data3.first, id3);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testEventSchedulerMoreEvents()
{
   EventScheduler scheduler;

   Event event1(9.0);
   Event event2(8.0);
   Event event3(7.0);
   Event event4(6.0);
   Event event5(5.0);
   Event event6(4.0);
   Event event7(3.0);
   Event event8(2.0);
   Event event9(1.0);

   auto id1 = scheduler.add(event1);
   auto id2 = scheduler.add(event2);
   auto id3 = scheduler.add(event3);
   auto id4 = scheduler.add(event4);
   auto id5 = scheduler.add(event5);
   auto id6 = scheduler.add(event6);
   auto id7 = scheduler.add(event7);
   auto id8 = scheduler.add(event8);
   auto id9 = scheduler.add(event9);

   TINYTEST_EQUAL(9, scheduler.size());
   TINYTEST_ASSERT(scheduler.check());

   
   auto data9 = scheduler.pop();
   TINYTEST_EQUAL(data9.first, id9);
   TINYTEST_EQUAL(data9.second.time(), event9.time());

   auto data8 = scheduler.pop();
   TINYTEST_EQUAL(data8.first, id8);
   TINYTEST_EQUAL(data8.second.time(), event8.time());

   auto data7 = scheduler.pop();
   TINYTEST_EQUAL(data7.first, id7);
   TINYTEST_EQUAL(data7.second.time(), event7.time());

   auto data6 = scheduler.pop();
   TINYTEST_EQUAL(data6.first, id6);
   TINYTEST_EQUAL(data6.second.time(), event6.time());

   auto data5 = scheduler.pop();
   TINYTEST_EQUAL(data5.first, id5);
   TINYTEST_EQUAL(data5.second.time(), event5.time());

   auto data4 = scheduler.pop();
   TINYTEST_EQUAL(data4.first, id4);
   TINYTEST_EQUAL(data4.second.time(), event4.time());

   auto data3 = scheduler.pop();
   TINYTEST_EQUAL(data3.first, id3);
   TINYTEST_EQUAL(data3.second.time(), event3.time());

   auto data2 = scheduler.pop();
   TINYTEST_EQUAL(data2.first, id2);
   TINYTEST_EQUAL(data2.second.time(), event2.time());

   auto data1 = scheduler.pop();
   TINYTEST_EQUAL(data1.first, id1);
   TINYTEST_EQUAL(data1.second.time(), event1.time());

   TINYTEST_EQUAL(0, scheduler.size());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testEventSchedulerMoreEvents2()
{
   RandomNumberGenerator rng;
   EventScheduler scheduler;
   const uint N = 400;

   std::vector<double> testTimes;
   testTimes.reserve(N);
   for (uint i = 0; i < N; ++i)
   {
      double t = rng.uniform(0.0, 10.0);
      testTimes.emplace_back(t);
      scheduler.add(Event(t));
   }

   std::sort(testTimes.begin(), testTimes.end());

   TINYTEST_EQUAL(N, scheduler.size());
   TINYTEST_ASSERT(scheduler.check());

   double time = 0.0;
   for (uint i = 0; i < N; ++i)
   {
      auto evnt = scheduler.pop();
      TINYTEST_ASSERT(evnt.second.time() >= time);
      time = evnt.second.time();
      TINYTEST_EQUAL(evnt.second.time(), testTimes[i]);

      TINYTEST_EQUAL(N-1-i, scheduler.size());
      TINYTEST_ASSERT(scheduler.check());
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(EventScheduler);
TINYTEST_ADD_TEST(testEventSchedulerEmpty);
TINYTEST_ADD_TEST(testEventSchedulerOneEvent);
TINYTEST_ADD_TEST(testEventSchedulerSecondEvent);
TINYTEST_ADD_TEST(testEventSchedulerReplaceEvent);
TINYTEST_ADD_TEST(testEventSchedulerRemoveEvent);
TINYTEST_ADD_TEST(testEventSchedulerMoreEvents);
TINYTEST_ADD_TEST(testEventSchedulerMoreEvents2);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------