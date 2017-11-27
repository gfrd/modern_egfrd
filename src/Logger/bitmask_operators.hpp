#pragma once

// (C) Copyright 2015 Just Software Solutions Ltd (Boost Software License)

#include <type_traits>

template<typename E>
struct enable_bitmask_operators {
   static const bool enable = false;
};

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable, E>::type
operator|(E lhs, E rhs) {
   typedef typename std::underlying_type<E>::type underlying;
   return static_cast<E>(
      static_cast<underlying>(lhs) | static_cast<underlying>(rhs));
}

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable, E>::type
operator&(E lhs, E rhs) {
   typedef typename std::underlying_type<E>::type underlying;
   return static_cast<E>(
      static_cast<underlying>(lhs) & static_cast<underlying>(rhs));
}

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable, E>::type
operator^(E lhs, E rhs) {
   typedef typename std::underlying_type<E>::type underlying;
   return static_cast<E>(
      static_cast<underlying>(lhs) ^ static_cast<underlying>(rhs));
}

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable, E>::type
operator~(E lhs) {
   typedef typename std::underlying_type<E>::type underlying;
   return static_cast<E>(
      ~static_cast<underlying>(lhs));
}

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable, E&>::type
operator|=(E& lhs, E rhs) {
   typedef typename std::underlying_type<E>::type underlying;
   lhs = static_cast<E>(
      static_cast<underlying>(lhs) | static_cast<underlying>(rhs));
   return lhs;
}

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable, E&>::type
operator&=(E& lhs, E rhs) {
   typedef typename std::underlying_type<E>::type underlying;
   lhs = static_cast<E>(
      static_cast<underlying>(lhs) & static_cast<underlying>(rhs));
   return lhs;
}

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable, E&>::type
operator^=(E& lhs, E rhs) {
   typedef typename std::underlying_type<E>::type underlying;
   lhs = static_cast<E>(
      static_cast<underlying>(lhs) ^ static_cast<underlying>(rhs));
   return lhs;
}
