#ifndef NFL_OPS_HPP
#define NFL_OPS_HPP

#include "nfl/meta.hpp"
#include "nfl/arch.hpp"
#include <tuple>
#include <iostream>
#include <limits>

// *********************************************************
//     Operations over polynomial coefficients
// *********************************************************

#include "nfl/arch/common.hpp"

namespace nfl {

#define DECLARE_UNARY_OPERATOR(SYM, NAME)\
  template<class T, size_t Degree, size_t NbModuli, size_t IndMod0>			\
  auto SYM(nfl::poly<T, Degree, NbModuli, IndMod0> const& op) -> decltype(ops::make_op<ops::NAME<nfl::poly<T, Degree, NbModuli, IndMod0>, CC_SIMD>>(op)) { \
  return ops::make_op<ops::NAME<nfl::poly<T, Degree, NbModuli, IndMod0>, CC_SIMD>>(op);\
}\
template<class Op, class...Args>\
auto SYM(nfl::ops::expr<Op, Args...> const& op) -> decltype(ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::poly_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op)) {\
  return ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::poly_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op);\
}

#define DECLARE_BINARY_OPERATOR(SYM, NAME)\
  template<class T, size_t Degree, size_t NbModuli, size_t IndMod0>			\
  auto SYM(nfl::poly<T, Degree, NbModuli, IndMod0> const& op0, nfl::poly<T, Degree, NbModuli, IndMod0> const& op1) -> decltype(ops::make_op<ops::NAME<nfl::poly<T, Degree, NbModuli, IndMod0>, CC_SIMD>>(op0, op1)) { \
  return ops::make_op<ops::NAME<nfl::poly<T, Degree, NbModuli, IndMod0>, CC_SIMD>>(op0, op1);\
}\
  template<class T, size_t Degree, size_t NbModuli, size_t IndMod0, class Op, class...Args> \
  auto SYM(nfl::poly<T, Degree, NbModuli, IndMod0> const& op0, nfl::ops::expr<Op, Args...> const& op1) -> decltype(ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::poly_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1)) { \
  return ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::poly_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1);\
}\
  template<class T, size_t Degree, size_t NbModuli, size_t IndMod0, class Op, class...Args> \
  auto SYM(nfl::ops::expr<Op, Args...> const& op0, nfl::poly<T, Degree, NbModuli, IndMod0> const& op1) -> decltype(ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::poly_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1)){ \
  return ops::make_op<ops::NAME<typename nfl::ops::expr<Op, Args...>::poly_type, typename nfl::ops::expr<Op, Args...>::simd_mode>>(op0, op1);\
}\
template<class Op0, class...Args0, class Op1, class... Args1>\
auto SYM(nfl::ops::expr<Op0, Args0...> const& op0,  nfl::ops::expr<Op1, Args1...> const& op1) -> decltype(ops::make_op<ops::NAME<typename nfl::ops::expr<Op0, Args0...>::poly_type, typename common_mode<typename Op0::simd_mode, typename Op1::simd_mode>::type>>(op0, op1)) {\
  static_assert(std::is_same<typename nfl::ops::expr<Op0, Args0...>::poly_type, typename nfl::ops::expr<Op1, Args1...>::poly_type>::value, "correct type combination");\
  return ops::make_op<ops::NAME<typename nfl::ops::expr<Op0, Args0...>::poly_type, typename common_mode<typename Op0::simd_mode, typename Op1::simd_mode>::type>>(op0, op1);\
}

namespace ops {

template<class Arg, class... Args>
  constexpr Arg first_of(Arg arg, Args... args) { return arg; }

template <class Op, class... Args>
struct expr {
  using simd_mode = typename Op::simd_mode;
  using poly_type = typename Op::poly_type;

  std::tuple<Args const&...> args;
  expr(Args const&... args) : args{args...}
  {}

  using value_type = typename poly_type::value_type;

  static constexpr size_t degree = first_of(Args::degree...);
  static constexpr size_t nmoduli = first_of(Args::nmoduli...);
  static constexpr size_t indmod0 = first_of(Args::indmod0...);

  template<class M, size_t... I>
  auto _load(size_t cm, size_t i, seq<I...>) const -> decltype(Op{}(std::get<I>(args).template load<M>(cm, i)..., cm+indmod0)) const
  {
    return Op{}(std::get<I>(args).template load<M>(cm, i)..., cm+indmod0);
  }

  template<class M>
  auto load(size_t cm, size_t i) const -> decltype(this->_load<M>(cm, i, typename gens<sizeof...(Args)>::type{}))
  {
    return _load<M>(cm, i, typename gens<sizeof...(Args)>::type{});
  }

  operator bool() const {
    for(size_t cm = 0; cm < nmoduli; ++cm) {
      constexpr size_t vector_size = simd_mode::template elt_count<value_type>::value;
      constexpr size_t vector_bound = degree / vector_size * vector_size;
      static_assert(vector_bound == degree, "no need for a footer");
      for(size_t j = 0; j < vector_bound; j+= vector_size) {
        alignas(32) value_type tmp[vector_size];
        simd_mode::store(tmp, load<simd_mode>(cm, j));
        for(size_t k = 0; k < vector_size; ++k)
          if(tmp[k])
            return true;
      }
    }
    return false;
  }

};

template<class T, class tag>
struct eqmod {
  using simd_mode = tag;
  using poly_type = T;
  template<class A>
  A operator()(A x, A y, size_t) const
  {
    return x == y;
  }
};

template<class T, class tag>
struct neqmod {
  using simd_mode = tag;
  using poly_type = T;
  template<class A>
  A operator()(A x, A y, size_t) const
  {
    return x != y;
  }
};

// Additions
// Modular addition plain and simple : add then test if sum is superior
// to the modulus in which case substract the modulus to the sum
// ASSUMPTION: x < p && y < p
// OUTPUT: x + y mod p
template<class T, class tag> struct addmod;
template<class T>
struct addmod<T, simd::serial> {
  using simd_mode = simd::serial;
  using poly_type = T;
  using value_type = typename T::value_type;
  value_type operator()(value_type x, value_type y, size_t cm) const
  {
    auto const p = poly_type::base.params.P[cm];
    ASSERT_STRICTMOD((__int128_t)x + y < 2*p);
    const value_type z = x + y;
    return z - ((z >= p) ? p : 0);
  }
};


// Modular subtract trick: if y is smaller than p return x+p-y else return x+p-(y%p)
// ASSUMPTION: x < p && y < p
// OUTPUT: x + y mod p
template<class T, class tag> struct submod;
template<class T>
struct submod<T, simd::serial> {
  using simd_mode = simd::serial;
  using poly_type = T;
  using value_type = typename T::value_type;  
  value_type operator()(value_type x, value_type y, size_t cm) const
  {
    auto const p = poly_type::base.params.P[cm];
    ASSERT_STRICTMOD(x<p && y<p);
    return addmod<T, simd_mode>{}(x, static_cast<value_type>(p-y), cm);
  }
};

template<class T, class tag> struct shoup;
template<class T>
struct shoup<T, simd::serial> {
  using simd_mode = simd::serial;
  using poly_type = T;
  using value_type = typename T::value_type;  
  template<class CM>
  value_type operator()(value_type x, value_type y, CM cm) const
  {
    static_assert(CM::omg_this_should_not_happen, "this should not happen");
    return value_type{};
  }
};

template<class T, class tag> struct compute_shoup;
template<class T>
struct compute_shoup<T, simd::serial> {
  using simd_mode = simd::serial;
  using poly_type = T;
  using value_type = typename T::value_type;
  using greater_value_type = typename T::greater_value_type;
  
  value_type operator()(value_type x, size_t cm) const
  {
    auto const p = poly_type::base.params.P[cm];
    while(x>=p) x-=p;
    return (((greater_value_type) x << poly_type::base.params.kModulusRepresentationBitsize) / p);
  }
};


// Multiplications
// Trivial (and very costly) modular multiplication with division
// OUTPUT: x * y mod p
template<class T, class tag> struct mulmod;
template<class T>
struct mulmod<T, simd::serial> {
  using simd_mode = simd::serial;
  using poly_type = T;
  using value_type = typename T::value_type;
  using greater_value_type = typename T::greater_value_type;
  static_assert(std::is_same<value_type, uint64_t>(), "this code has only been tested for uint64_t");
  value_type operator()(value_type x, value_type y, size_t cm) const
  {
    auto const p = poly_type::base.params.P[cm];
    auto const shift = poly_type::base.params.kModulusRepresentationBitsize;
    auto const s0 = shift-poly_type::base.params.kModulusBitsize;
    ASSERT_STRICTMOD((x<p) && (y<p));
    greater_value_type res = (greater_value_type)x * y;
    greater_value_type q = ((greater_value_type)poly_type::base.params.Pn[cm] * (res >> shift)) + (res<<s0) ;
    value_type r  = res - (q>>shift) * p;
    if (r >= p) { r -= p; }
    ASSERT_STRICTMOD(r == ((greater_value_type)(x) * y) % p);
    return r;
  }
};

// Modular multiplication: much faster alternative
// Works only if yshoup = ((uint128_t) y << 64) / p (for T=uint64_t)
// has been pre-computed AND p is small enough (2 bits less than limb).
// OUTPUT: x * y lazymod p (if x and y are below p the result is x * y mod p)
template<class T, class tag> struct mulmod_shoup;

template<class T>
struct mulmod_shoup<T, simd::serial> {
  using simd_mode = simd::serial;
  using poly_type = T;
  using value_type = typename T::value_type;
  using greater_value_type = typename T::greater_value_type;
  value_type operator()(value_type x, value_type y, value_type yprime, size_t cm) const
  {
    auto const p = poly_type::base.params.P[cm];
    ASSERT_STRICTMOD(x<p && y<p);
    greater_value_type res;
    value_type q = ((greater_value_type) x * yprime) >> poly_type::base.params.kModulusRepresentationBitsize;
    res = x * y - q * p;
    ASSERT_STRICTMOD(res - ((res>=p) ? p : 0) == ((greater_value_type)(x%p)*(y%p))%p);
    return res - ((res>=p) ? p : 0);
  }
};

//
// Detect various patterns
//

// generic implem to specialize
template<class Op, class... Args>
struct _make_op {
  expr<Op, Args...> operator()(Args const&... args) const {
    return expr<Op, Args...>{args...};
  }
};

template<class Op, class... Args>
auto make_op(Args const&... args) -> decltype(_make_op<Op, Args...>{}(args...))
{
  return _make_op<Op, Args...>{}(args...);
}

template<class... Args>
using retag = typename common_mode<typename Args::simd_mode...>::type;


// detect mul shoup
template<class tag0, class tag1, class type, class Arg0, class Arg1, class Arg2>
struct _make_op<shoup<type, tag0>, expr<mulmod<type, tag1>,Arg0, Arg1>, Arg2> {
  expr<mulmod_shoup<type, retag<Arg0, Arg1, Arg2>>, Arg0, Arg1, Arg2>
  operator()(expr<mulmod<type, tag1>, Arg0, Arg1> const& from0, Arg2 const& from1) const {
    return expr<mulmod_shoup<type, retag<Arg0, Arg1, Arg2>>, Arg0, Arg1, Arg2>{
      std::get<0>(from0.args),
      std::get<1>(from0.args),
      from1
    };
  }
};


} // ops

} // nfl

#ifdef NFL_OPTIMIZED

#include "nfl/opt/ops.hpp"

// #if defined __SSE4_2__ && defined NTT_SSE
// #include "nfl/opt/arch/sse.hpp"
// #endif

// #if defined __AVX2__ && defined NTT_AVX2
// #include "nfl/opt/arch/avx2.hpp"
// #endif

#endif

#endif
