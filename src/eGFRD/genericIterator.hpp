#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <utility>
#include <iterator>
#include <typeinfo>

// --------------------------------------------------------------------------------------------------------------------------------

namespace gi
{
   // --------------------------------------------------------------------------------------------------------------------------------

   // default selector, just select E
   template<typename E>
   struct Select
   {
      E& operator()(E& e) const { return e; }
   };

   // pair first selector, select E from std::pair<E,x>
   template<typename E>
   struct SelectFirst
   {
      template <typename T1>
      E& operator()(std::pair<E, T1>& t) const { return t.first; }
   };

   // pair second selector, select E from std::pair<x,E>
   template<typename E>
   struct SelectSecond
   {
      template <typename T0>
      E& operator()(std::pair<T0, E>& t) const { return t.second; }
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   /// template wrapper class for iterators holding type E, independent of container/implementation
   template<typename TContainer, typename E, typename TSelector>
   class Itr2 : public std::iterator<std::forward_iterator_tag, typename std::remove_cv<E>::type, std::ptrdiff_t, E*, E&>
   {
   public:
      Itr2(typename TContainer::iterator i, const TSelector sel) : citr_(i), sel_(sel) { }

      void operator++() { ++citr_; }
      E& operator*() const { return sel_(*citr_); }
      Itr2& operator=(Itr2& o) { citr_ = o.citr_; return *this; }

      bool operator==(const Itr2& o) const { return (citr_ == o.citr_); }
      bool operator!=(const Itr2& o) const { return !(*this == o); }

   private:
      typename TContainer::iterator citr_;
      const TSelector sel_;
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   // helper class that holds begin and end iterators of a container (using known implementation/container, which is faster)
   //template<typename TContainer, typename E = std::remove_cv_t<typename TContainer::value_type>, typename TSelector = Select<E>>
   template<typename TContainer, typename E = typename TContainer::value_type, typename TSelector = Select<E>>
   class iteratorRange
   {
   public:

      iteratorRange(const TContainer& c) : begin_(std::begin(const_cast<TContainer&>(c)), TSelector()), end_(std::end(const_cast<TContainer&>(c)), TSelector()) { }

      Itr2<TContainer, E, TSelector> begin() const { return begin_; }
      Itr2<TContainer, E, TSelector> end() const { return end_; }

   private:
      const Itr2<TContainer, E, TSelector>      begin_;
      const Itr2<TContainer, E, TSelector>      end_;
   };
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace agi
{
   // --------------------------------------------------------------------------------------------------------------------------------

   // abstract base class for implementing Itr
   template<typename E>
   class ItrBase
   {
   public:
      ItrBase() = default;
      virtual ~ItrBase() = default;
      virtual void operator++() {}
      virtual E& operator*() const = 0; // { return E(); }
      virtual ItrBase* clone() const = 0; // { return new ItrBase(*this); }
      bool operator==(const ItrBase& o) const { return typeid(*this) == typeid(o) && equal(o); }

   protected:
      virtual bool equal(const ItrBase&) const { return true; }
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   /// dynamic generic wrapper class for iterators holding type E, independent of container/implementation
   template<typename E>
   class Itr : public std::iterator<std::forward_iterator_tag, typename std::remove_cv<E>::type, std::ptrdiff_t, E*, E&>
   {
   public:
      Itr() : itr_(nullptr) {}
      explicit Itr(ItrBase<E>* impl) : itr_(impl) {}

      ~Itr() { delete itr_; }
      Itr(const Itr& o) : itr_(o.itr_->clone()) {}
      Itr& operator=(const Itr& o)
      {
         if (itr_ != o.itr_) { delete itr_; itr_ = o.itr_->clone(); }
         return *this;
      }

      Itr& operator++() { ++(*itr_); return *this; }
      E& operator*() const { return *(*itr_); }
      bool operator==(const Itr& o) const { return (itr_ == o.itr_) || (*itr_ == *o.itr_); }
      bool operator!=(const Itr& o) const { return !(*this == o); }
   private:
      ItrBase<E>* itr_;
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   // implementation of ItrBase for all elements (no filter/adapters) of given container
   template<typename TContainer, typename E, typename TSelector>
   class ItrAll : public ItrBase<E>
   {
   public:
      ItrAll(const TContainer& c, typename TContainer::iterator i, const TSelector sel = TSelector()) : ItrBase<E>(), c_(c), citr_(i), sel_(sel) { }
      virtual void operator++() override { ++citr_; }
      virtual E& operator*() const override { return sel_(*citr_); }
      virtual ItrBase<E>* clone() const override { return new ItrAll(*this); }
   protected:
      virtual bool equal(const ItrBase<E>& o) const override
      {
         const ItrAll& o2 = static_cast<const ItrAll&>(o);
         return &c_ == &o2.c_ && citr_ == o2.citr_;
      }
   private:
      const TContainer& c_;
      typename TContainer::iterator citr_;
      const TSelector sel_;
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   // helper methods for getting a begin iterator of a container (using generic iterator implementation)
   template<typename TContainer, typename E, typename TSelector = gi::Select<E>>
   Itr<E> begin(const TContainer& c)
   {
      return Itr<E>(new ItrAll<TContainer, E, TSelector>(c, const_cast<TContainer&>(c).begin()));
   }

   // helper methods for getting a end iterator of a container (using generic iterator implementation)
   template<typename TContainer, typename E, typename TSelector = gi::Select<E>>
   Itr<E> end(const TContainer& c)
   {
      return Itr<E>(new ItrAll<TContainer, E, TSelector>(c, const_cast<TContainer&>(c).end()));
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // helper class that holds begin and end iterators of a container (using generic iterator implementation)
   template<typename E, typename TSelector = gi::Select<E>>
   class iteratorRange
   {
   public:

      template<typename TContainer>
      iteratorRange(const TContainer& c) : begin_(agi::begin<TContainer, E, TSelector>(c)), end_(agi::end<TContainer, E, TSelector>(c)) { }

      Itr<E> begin() const { return begin_; }
      Itr<E> end() const { return end_; }

   private:
      const Itr<E>      begin_;
      const Itr<E>      end_;
   };

   // --------------------------------------------------------------------------------------------------------------------------------
}
