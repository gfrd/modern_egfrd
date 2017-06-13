#ifndef REACTANTS_HPP
#define REACTANTS_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "DefsEgfrd.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TID1, typename TID2>
class Reactants
{
   using value1_type = TID1;
   using value2_type = TID2;

   explicit Reactants() noexcept {}       // default constructor (private, available only to Persistence)

public:
   explicit Reactants(const value1_type one) { ASSERT(one);  item1_ = one; item2_ = value2_type(); }

   template<typename T1 = value1_type, typename T2 = value2_type>
   explicit Reactants(T1 one, T1 two, typename std::enable_if<std::is_same<T1, T2>::value>::type* = nullptr)
   {
      ASSERT(one);
      ASSERT(two);

      if (one <= two)
      {
         item1_ = one;
         item2_ = two;
      }
      else
      {
         item1_ = two;
         item2_ = one;
      }
   }

   template<typename T1 = value1_type, typename T2 = value2_type>
   explicit Reactants(T1 one, T2 two, typename std::enable_if<!std::is_same<T1, T2>::value>::type* = nullptr)
   {
      ASSERT(one);
      ASSERT(two);

      item1_ = one;
      item2_ = two;
   }

   std::size_t size() const { return !item1_ ? 0 : !item2_() ? 1 : 2; }

   const value1_type& item1() const { ASSERT(0 < size()); return item1_; }
   const value2_type& item2() const { ASSERT(1 < size()); return item2_; }

   bool operator<(const Reactants& rhs) const
   {
      if (size() < rhs.size()) return true;
      if (size() > rhs.size()) return false;
      switch (size())
      {
      case 0: return false;
      case 1: return item1_ < rhs.item1_;
      case 2: return item1_ == rhs.item1_ ? (item2_ < rhs.item2_) : (item1_ < rhs.item1_);
      default: throw illegal_state("nop");
      }
   }

   bool operator==(const Reactants& rhs) const
   {
      if (rhs.size() != size()) return false;
      switch (size())
      {
      case 0: return true;
      case 1: return item1_ == rhs.item1_;
      case 2: return item1_ == rhs.item1_ && item2_ == rhs.item2_;
      default: throw illegal_state("nop");
      }
   }

   bool operator!=(const Reactants& rhs) const { return !operator==(rhs); }

protected:
   friend class Persistence;

   value1_type item1_;
   value2_type item2_;
};

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TID1, typename TID2>
inline std::ostream& operator<<(std::ostream& stream, const Reactants<TID1, TID2>& r)
{
   stream << "Reactants{";
   switch (r.size())
   {
   case 1: stream << r.item1(); break;
   case 2: stream << r.item1() << "+" << r.item2(); break;
   default: ;
   }
   stream << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* REACTANTS_HPP */

