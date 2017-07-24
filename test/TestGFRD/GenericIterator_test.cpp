#include "GenericIterator_test.hpp"
#include "genericIterator.hpp"
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "DefsEgfrd.hpp"
#include <makeString.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

int testGIVector()
{
   std::vector<int>     c;
   for (uint i = 1; i <= 10; ++i)
      c.push_back(i);

   TINYTEST_EQUAL(10, c.size());
   TINYTEST_ASSERT(c.begin() != c.end());

   auto range = agi::iteratorRange<int>(c);
   for (const auto &ii : range) std::cout << ii;

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testGISet()
{
   std::set<int>     c;
   for (uint i = 1; i <= 10; ++i)
      c.insert(i * 3);

   TINYTEST_EQUAL(10, c.size());
   TINYTEST_ASSERT(c.begin() != c.end());

   auto range = agi::iteratorRange<const int>(c);
   for (auto i : range) std::cout << i;

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testGIPairVector()
{
   std::vector<std::pair<int, int>>     c;
   for (uint i = 1; i <= 10; ++i)
      c.push_back(std::make_pair(i, i + 1));

   TINYTEST_EQUAL(10, c.size());
   TINYTEST_ASSERT(c.begin() != c.end());

   auto range = agi::iteratorRange<std::pair<int, int>>(c);

   for (auto i = range.begin(); i != range.end(); ++i)
   {
      (*i).first++;   // its const
      std::cout << "F=" << (*i).first << " S=" << (*i).second;
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testGISecondSelect()
{
   std::vector<std::pair<int, double>>     c;
   for (uint i = 0; i <= 7; ++i)
      c.push_back(std::make_pair(i, std::sin(i * 2.0 * 3.1415 / 8)));

   agi::iteratorRange<double, gi::SelectSecond<double>> range(c);
   for (auto i = range.begin(); i != range.end(); ++i)
   {
      std::cout << (*i);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testGITiming()
{
   std::vector<int> c;
   c.reserve(1000);
   for (uint i = 1; i <= 1000; ++i)
      c.push_back(i);

   TINYTEST_TIME_START();
   for (const auto &i : c) TINYTEST_ASSERT(i > 0 && i <= 1000);
   TINYTEST_TIME_STOP("VectorItr:");

   TINYTEST_TIME_START();
   auto range = agi::iteratorRange<int>(c);
   for (const auto &i : range) TINYTEST_ASSERT(i > 0 && i <= 1000);
   TINYTEST_TIME_STOP("AbstractItr:");

   TINYTEST_TIME_START();
   auto range = gi::iteratorRange<std::vector<int>>(c);
   for (const auto &i : range) TINYTEST_ASSERT(i > 0 && i <= 1000);
   TINYTEST_TIME_STOP("GenericItr:");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testGIGeneric()
{
   typedef std::vector<std::pair<int, double>> container;

   container     c;
   for (uint i = 0; i <= 7; ++i)
      c.push_back(std::make_pair(i, std::sin(i * 2.0 * 3.1415 / 8)));

   gi::iteratorRange<container, double, gi::SelectSecond<double>> range(c);
   for (auto i = range.begin(); i != range.end(); ++i)
   {
      std::cout << (*i);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testGIFind1()
{
   typedef std::vector<uint> container;

   container     c;
   for (uint i = 1; i <= 9; ++i)
      c.push_back(i);

   gi::iteratorRange<container> range(c);
   auto f1 = std::find(range.begin(), range.end(), 7u);
   TINYTEST_ASSERT(f1 != range.end());

   auto f2 = std::find(range.begin(), range.end(), 9u);
   TINYTEST_ASSERT(f2 != range.end());

   auto f3 = std::find(range.begin(), range.end(), 19u);
   TINYTEST_ASSERT(f3 == range.end());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testGIFind2()
{
   typedef std::vector<std::pair<int, std::string>> container;

   container     c;
   for (uint i = 1; i <= 9; ++i)
   {
      std::string name(make_string() << "Test:" << i);
      c.push_back(std::make_pair(i, name));
   }

   gi::iteratorRange<container, std::string, gi::SelectSecond<std::string>> range(c);
   for (auto i = range.begin(); i != range.end(); ++i)
   {
      std::cout << (*i);
   }

   auto f1 = std::find(range.begin(), range.end(), "Test:7");
   TINYTEST_ASSERT(f1 != range.end());

   auto f2 = std::find(range.begin(), range.end(), "Test:9");
   TINYTEST_ASSERT(f2 != range.end());

   auto f3 = std::find(range.begin(), range.end(), "Test:100");
   TINYTEST_ASSERT(f3 == range.end());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(GenericIterator);
TINYTEST_ADD_TEST(testGIVector);
TINYTEST_ADD_TEST(testGISet);
TINYTEST_ADD_TEST(testGIPairVector);
TINYTEST_ADD_TEST(testGISecondSelect);
TINYTEST_ADD_TEST(testGITiming);
TINYTEST_ADD_TEST(testGIGeneric);
TINYTEST_ADD_TEST(testGIFind1);
TINYTEST_ADD_TEST(testGIFind2);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
