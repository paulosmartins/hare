#ifndef __TMP_HPP__
#define __TMP_HPP__
#include <functional>

template<size_t i>
struct static_log2
{
  static constexpr size_t value = 1 + static_log2<i/2>::value;
};

template<>
struct static_log2<1>
{
  static constexpr size_t value = 0;
};

#endif
