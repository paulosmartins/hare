#ifndef common_hpr_h
#define common_hpr_h

// #include "params.hpp"
#include "tmp.hpp"

namespace hpr
{
  namespace commons
  {
    namespace tools
    {
      // to test equality of elements of R
      // not critical: only used to check correctness of decryption
      template<class P>
      bool equal(const P &p1, const P &p2)
      {
	auto* i2 = p2.begin();
	for(auto* i1 = p1.begin(); i1 < p1.end(); ++i1, ++i2)
	  if(*i1 != *i2)
	    return false;
	return true;
      }

      template<class C> C MIN(C a, C b)
      {
	return (a > b) ? b : a;
      }

      // time measurement function ; from NFL
      template <class T>
      double get_time_us(T const& start, T const& end, uint32_t N)
      {
	auto diff = end-start;
	return (long double)(std::chrono::duration_cast<std::chrono::microseconds>
			     (diff).count())/N;
      }

      // apply modulo t to an element of R (with residues in [0,t)
      // not critical: only used for testing correctness of decryption
      template<class P, class parameters>
      void modt(P &p)
      {
	auto* iter = p.begin();
	for(size_t cm = 0; cm < P::nmoduli; cm++)
	  {
	    typename P::value_type mod = p.get_modulus(cm);
	    for(size_t i = 0; i < P::degree; i++, ++iter)
	      *iter = (*iter > ((mod-1)>>1)) ? parameters::t-((mod-*iter)&parameters::mmask_t)
		: *iter&parameters::mmask_t;
	  }
      }

      // set value v as coefficient of degree d of element p in R
      // not critical: only used when doing precomputations
      template<class P>
      void set_1coeff_mpz(P &p, const mpz_class &v, size_t d)
      {
	auto* iter = p.begin()+d;
	for(size_t cm = 0; cm < P::nmoduli; cm++)
	  {
	    *iter = static_cast<typename P::value_type>(mpz_fdiv_ui(v.get_mpz_t(),
								    p.get_modulus(cm)));
	    iter += P::degree;
	  }
      }

      // multiply polynomial p by v
      // not critical: only used when precomputing relinearization key rlk
      template<class P>
      void mul_mpz_class(P &p, const mpz_class &v)
      {
	mpz_class tmp;
	nfl::ops::mulmod<P, nfl::simd::serial> op_mul;
	auto* iter = p.begin();
	for(size_t cm = 0; cm < P::nmoduli; cm++)
	  {
	    auto mod = p.get_modulus(cm);
	    typename P::value_type val = static_cast<typename P::value_type>(mpz_mod_ui(
											tmp.get_mpz_t(), v.get_mpz_t(), mod));
	    for(size_t i = 0; i < P::degree; i++)
	      {
		*iter = op_mul(val, *iter, cm+P::indmod0);
		++iter;
	      }
	  }
      }

      // multiply coefficient p(cm,i) by v
      // not critical: only used when precomputing rlk
      template<class P>
      void mul_val(P &p, typename P::value_type v, size_t cm, size_t i)
      {
	nfl::ops::mulmod<typename P::value_type, nfl::simd::serial> op_mul;
	p.data()[cm*P::degree+i] = op_mul(p(cm, i), v, cm+P::indmod0);
      }

      // neg coefficients
      // not critical so far
      template<class P>
      void neg(P &p)
      {
	constexpr size_t degree = P::degree;
	auto* it = p.begin();
	for(uint cm = 0; cm < P::nmoduli; cm++)
	  {
	    typename P::value_type mod = p.get_modulus(cm);
	    for(uint n = 0; n < degree; n++, ++it)
	      *it = mod-(*it);
	  }
      }

      // first to last: DEGREE coefficients, ie c_0,...,c_(deg-1)
      // put these into p, ie p <- (c_0,...,c_(deg-1)) mod q1, ..., (c_0,...,c_(deg-1)) mod qk
      // critical: used in homomorphic mult.
      template <class P, class It>
      void set_partial(P &p, It first, It last)
      {
	auto* it = p.begin();
	for (uint i = 0; i < P::nmoduli; i++)// it != p.end())
	  it = std::copy(first, last, it);
      }

      template <class P, class It>
      inline void set_partial_signed(P &p, It first, It last)
      {
	using main_type = typename P::value_type;
	using main_greater_type = typename P::greater_value_type;
	using main_signed_type = typename P::signed_value_type;

	auto* it = p.begin();
	for(uint cm = 0; cm < P::nmoduli; cm++)
	  {
	    main_signed_type mod_s = static_cast<main_signed_type>(p.get_modulus(cm));
	    for(It it2 = first; it2 != last; ++it2, ++it)
	      {
		assert(mod_s+(*it2)>=0);
		*it = (*it2 >= 0) ? static_cast<typename P::value_type>(*it2) :
		  static_cast<main_type>(mod_s+(*it2));
	      }
	  }
      }

      // critical: used in homomorphic mult.
      template <class P>
      inline void spread_value(P &p, typename P::value_type v)
      {
	std::fill(p.begin(), p.end(), v);
      }

    }

    namespace mem_manag
    {
      // convenient function from NFL
      template <class T, size_t Align, class... Args>
      T* alloc_aligned(size_t n, Args&& ... args)
      {
	T* ret;
	if (posix_memalign((void**) &ret, Align, sizeof(T)*n) != 0)
	  {
	    throw std::bad_alloc();
	  }
	for (size_t i = 0; i < n; i++)
	  {
	    new (&ret[i]) T(std::forward<Args>(args)...);
	  }
	return ret;
      }

      // convenient function from NFL
      template <class T>
      void free_aligned(size_t n, T* p)
      {
	for (size_t i = 0; i < n; i++)
	  {
	    p[i].~T();
	  }
	free(p);
      }
    }

    namespace reductions
    {
      //  Inner Modular Reduction with output in [0,2*mod)
      template<typename poly_type>
      inline typename poly_type::value_type IMR
      (typename poly_type::greater_value_type in,
       typename poly_type::value_type w0,
       typename poly_type::value_type mod)
      {
	using main_type = typename poly_type::value_type;
	using main_greater_type = typename poly_type::greater_value_type;
	using main_signed_type = typename poly_type::signed_value_type;

	size_t main_bs = poly_type::base.params.kModulusBitsize;
	size_t main_bs_full = poly_type::base.params.kModulusRepresentationBitsize;
	size_t s0 = main_bs_full - main_bs;

	main_greater_type q = w0*(in>>main_bs_full)+(in<<s0);
	main_type res = static_cast<main_type>(in)-static_cast<main_type>
	  (q>>main_bs_full)*mod;
	return res;
      }

      //  full IMR: output in [0,mod)
      template<typename poly_type>
      inline typename poly_type::value_type IMR_full
      (typename poly_type::greater_value_type in,
       typename poly_type::value_type w0,
       typename poly_type::value_type mod)
      {
	using main_type = typename poly_type::value_type;
	using main_greater_type = typename poly_type::greater_value_type;
	using main_signed_type = typename poly_type::signed_value_type;
	size_t main_bs = poly_type::base.params.kModulusBitsize;
	size_t main_bs_full = poly_type::base.params.kModulusRepresentationBitsize;
	size_t s0 = main_bs_full - main_bs;

	main_greater_type q = w0*(in>>main_bs_full)+(in<<s0);
	main_type res = static_cast<main_type>(in)-static_cast<main_type>
	  (q>>main_bs_full)*mod;
	return (res >= mod) ? res-mod : res;
      }
    }


    namespace key
    {
      template<class poly_hpr, class parameters>
      struct key_chain
      {
	typedef typename poly_hpr::_poly_digit poly_digit;
	typedef typename poly_hpr::Pq Pq;
	typedef typename poly_hpr::Pb Pb;
	typedef typename poly_hpr::Psk Psk;

	poly_digit* _sk; // private key (poly with small coeff)
	poly_digit* _s2;
	poly_hpr* _pk; // (b = -a*s+e, a), a: uniform poly, e: small noise

	void allocate();
	void copy_from(const key_chain &other);
	key_chain();
	key_chain(const key_chain &other);
	key_chain& operator=(const key_chain &other);
	~key_chain();

	template<typename be_tool_set>
	void generate(be_tool_set*);
      };

      template<class poly_hpr, class parameters>
      void key_chain<poly_hpr, parameters>::allocate()
      {
	_s2 = new poly_digit();
	_sk = new poly_digit();
	_pk = new poly_hpr();
      }

      template<class poly_hpr, class parameters>
      void key_chain<poly_hpr, parameters>::copy_from
      (const key_chain<poly_hpr, parameters> &other)
      {
	*_s2 = *other._s2;
	*_sk = *other._sk;
	*_pk = *other._pk;
      }
      
      template<class poly_hpr, class parameters>
      key_chain<poly_hpr, parameters>::key_chain()
      {
	allocate();
      }

      template<class poly_hpr, class parameters>
      key_chain<poly_hpr, parameters>::key_chain
      (const key_chain<poly_hpr, parameters> &other)
      {
	allocate();
	copy_from(other);
      }

      template<class poly_hpr, class parameters>
      key_chain<poly_hpr, parameters>&
      key_chain<poly_hpr, parameters>::operator=
      (const key_chain<poly_hpr, parameters> &other)
      {
	copy_from(other);
	return *this;
      }
      
      template<class poly_hpr, class parameters>
      key_chain<poly_hpr, parameters>::~key_chain()
      {
	delete _s2;
	delete _sk;
	delete _pk;
      }


      template<class poly_hpr, class parameters>
      template<typename be_tool_set>
      void key_chain<poly_hpr, parameters>::generate(be_tool_set *ts)
      {
	_pk->set_mat_set_conv(ts);
	_sk->set_mat_set_conv(ts);
	_s2->set_mat_set_conv(ts);

	// creates keys in ntt mode

	nfl::uniform unif;
	nfl::non_uniform nunif(parameters::B_key+1);

	_pk->set0();

	// generates noise for first component
	_pk->__digits[0]->gen_noise(0); // pk = (e, 0)

	// generates uniform poly for second component
	int  i=0;
	for(typename poly_hpr::iterator it = _pk->begin(); it != _pk->end(); ++it)
	  {
	    (*it)->__cq[1] = unif;
	    if(it != _pk->end()-1)
	      (*it)->exact_conv_gmp_q_to_bsk(1);
	  }

	// pk = (e, a)

	// generating secret key
	_sk->set0();
	_sk->__cq[0] = nunif;

	// copying from q to bsk
	_sk->spread_smallpoly_q2bsk(0);

	_sk->ntt(0);
	_sk->sub(1, _sk, 0); // sk = (s, -s)

	_s2->set0();
	_s2->mul(0, _sk, 0, _sk, 0);  // s2 = (s^2, 0)

	_pk->ntt();
	_pk->mulDigit_acc(0, _pk, 1, _sk, 1); // pk = (e - a*s, a)

	_pk->inv_ntt(0);
	_pk->carry_propagation_1side(0, 0, _pk->__d-1);
	_pk->ntt(0);

	_pk->fast_modq();
      }
    }

    namespace raw_pol
    {
      // struct for polynomials modulo m_tilde (NTT not required: no need for polynomial product here)
      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      struct PolMT
      {
	using iterator = type_mtilde*;
	using const_iterator = type_mtilde const*;
	static constexpr size_t _bs_mt = bs_mt;

	static const type_mtilde mtilde_half = static_cast<type_mtilde>(1)<<(bs_mt-1);
	static const main_type mtilde = static_cast<main_type>(1)<<bs_mt;

	type_mtilde* __data;

	PolMT()
	{
	  __data = (type_mtilde*)_mm_malloc(degree*sizeof(type_mtilde),32);
	  std::fill(__data, __data + degree, 0);
	};
	
	PolMT(const PolMT &other)
	{
	  __data = (type_mtilde*)_mm_malloc(degree*sizeof(type_mtilde),32);
	  std::copy(other.__data, other.__data+degree, __data);
	}

	PolMT& operator=(const PolMT &other)
	{
	  std::copy(other.__data, other.__data+degree, __data);
	  return *this;
	}

	~PolMT()
	{
	  _mm_free(__data);
	}

	iterator begin()
	{
	  return __data;
	}
	const_iterator begin() const
	{
	  return __data;
	}
	iterator end()
	{
	  return __data+degree;
	}
	const_iterator end() const
	{
	  return __data+degree;
	}

	void set(const PolMT &pol);
	void set0()
	{
	  std::fill(__data, __data + degree, 0);
	}

	// base conversion
	template <class P>
	//#if (BS_MT_FAST == 0)
	//            void fast_base_conversion(const P &p, const std::array<type_mtilde, P::nmoduli> &v);
	//#else
	void fast_base_conversion(const P &p);
	//#endif

	// centered remainders
	void get_modc(std::array<main_signed_type, degree> &v) const;

	// print for debug
	void print() const;
      };

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      void PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::set(
												   const PolMT &pol)
      {
	auto iter = pol.begin();
	for(auto it = begin(); it < end(); ++it, ++iter)
	  *it = static_cast<type_mtilde>(*iter);
      }


      // for base conversion, it will be from q
      // the precomputations will be (-1/q1, ..., -1/qk) mod mt
      /*#if (BS_MT_FAST == 0)
	template<class P>
	void PolMT::fast_base_conversion(const P &p, const std::array<type_mtilde, P::nmoduli> &v)
	{
	auto citer = p.begin();
	for(auto itv = v.begin(); itv != v.end(); ++itv)
	for(auto it = begin(); it != end(); ++it, ++citer)
	*it += (*itv)*static_cast<type_mtilde>(*citer);
	}
	#else
      */
      // here, qi = 1 mod mt for all i
      // then, (-1/q1, ..., -1/qk) mod mt = (-1,...)
      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      template<class P>
      inline void
      PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::fast_base_conversion(
													       const P &p)
      {
	auto citer = p.begin();
	for(auto citer = p.begin(); citer != p.end();)
	  for(auto it = begin(); it != end(); ++it, ++citer)
	    *it += (~ static_cast<type_mtilde>(*citer))+1;
      }
      //#endif

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      inline void
      PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::get_modc(
												   std::array<main_signed_type, degree> &v) const
      {
	// if residue > mod/2, return residue - mod
	// handled through simple conversion to signed type (except for mtilde_half)
	auto iter = begin();
	for(auto iterv = v.begin(); iterv != v.end(); ++iterv)
	  {
	    type_mtilde val = static_cast<signed_type_mtilde>(*iter++);
	    *iterv = (val == mtilde_half) ? static_cast<main_signed_type>
	      (val) : static_cast<main_signed_type>(static_cast<signed_type_mtilde>(val));
	  }
      }

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      void PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::print()
	const
      {
	std::cout<<"{ ";
	for(uint i = 0; i < degree; i++)
	  {
	    std::cout<<(uint32_t)(__data[i]);
	    if(i<degree-1)
	      std::cout<<", ";
	    else
	      std::cout<<" }";
	  }
	std::cout<<std::endl;
      }



      // Polynomial class for messages: only used for convenient checkings during tests
      template<typename parameters, size_t degree>
      struct PolT
      {
	using type_t = typename parameters::type_t;
	using iterator = typename std::array<type_t, degree>::iterator;
	using const_iterator = typename std::array<type_t, degree>::const_iterator;

	static constexpr size_t _dim = degree;

	std::array<type_t, _dim> __data; // contains polynomial coeff.

	PolT()
	{
	  __data.fill(0);
	};

	~PolT() {};

	PolT(const PolT &other)
	{
	  std::copy(&other.__data[0], &other.__data[0]+degree, &__data[0]);
	}

	PolT& operator=(const PolT &other)
	{
	  std::copy(&other.__data[0], &other.__data[0]+degree, &__data[0]);
	  return *this;
	}

	iterator begin()
	{
	  return std::begin(__data);
	}
	iterator end()
	{
	  return std::end(__data);
	}

	const_iterator begin() const
	{
	  return std::begin(__data);
	}
	const_iterator end() const
	{
	  return std::end(__data);
	}

	// setters
	void set0()
	{
	  __data.fill(0);
	};
	void set_unif();
	template <class P>
	void set(const P &p);
	template <class T>
	void set_val(T);
	template <class P>
	void get(P &p) const;

	// ope
	bool operator==(const PolT &other) const;
	void mul(const PolT *other);
	void add(const PolT *other);
	void modt();

	// print for debug
	void print(size_t i0, size_t i1) const;
      };

      // fill current PolT with p
      template<typename parameters, size_t degree>
      template <class P>
      void PolT<parameters, degree>::get(P &p) const
      {
	for(uint i = 0; i < _dim; i++)
	  for(uint cm = 0; cm < P::nmoduli; cm++)
	    *(p.begin()+(cm<<static_log2<degree>::value)+i) =
	      static_cast<typename P::value_type>(__data[i]);
      }

      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::set_unif()
      {
	nfl::fastrandombytes((unsigned char *)__data.data(), sizeof(*this));
	this->modt();
      }

      // fill p with current PolT
      template<typename parameters, size_t degree>
      template<class P>
      void PolT<parameters, degree>::set(const P &p)
      {
	auto iter = p.begin();
	typename P::value_type mod = p.get_modulus(0);
	for(auto it = begin(); it != end(); ++it, ++iter)
	  *it = (*iter < (mod>>1)) ? static_cast<type_t>(*iter)&parameters::mask_t :
	    parameters::t-static_cast<type_t>(mod-*iter)&parameters::mask_t;
      }

      template<typename parameters, size_t degree>
      template<class T>
      void PolT<parameters, degree>::set_val(T v)
      {
	this->set0();
	__data[0] = (type_t)v;
      }

      template<typename parameters, size_t degree>
      bool PolT<parameters, degree>::operator==(const PolT<parameters, degree> &other)
	const
      {
	auto iter_ot = other.begin();
	for(auto it = begin(); it != end(); ++it, ++iter_ot)
	  {
	    if (*it != *iter_ot)
	      return false;
	  }
	return true;
      }

      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::modt()
      {
	for(auto it = begin(); it != end(); ++it)
	  *it = (*it)&parameters::mask_t;
      }

      // quadratic multiplication: only used for checking
      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::mul(const PolT<parameters, degree> *other)
      {
	std::array<type_t, _dim> tmp;
	tmp.fill(0);
	for(size_t k = 0; k <= _dim-1; ++k)
	  for(size_t i = 0; i <= k; ++i)
	    tmp[k] += __data[i] * other->__data[k-i];
	for(size_t k = _dim; k <= 2*_dim-1; ++k)
	  for(size_t i = k-_dim+1; i <= _dim-1; ++i)
	    tmp[k-_dim] -= __data[i]*other->__data[k-i];
	for(size_t k = 0; k <= _dim-1; ++k)
	  __data[k] = tmp[k];

	this->modt();
      }

      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::add(const PolT<parameters, degree> *other)
      {
	for(size_t k = 0; k <= _dim-1; ++k)
	  __data[k] += other->__data[k];

	this->modt();
      }

      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::print(size_t i0, size_t i1) const
      {
	std::cout<<"{ ";
	for(auto it = this->begin()+i0; it != this->begin()+i1; ++it)
	  {
	    std::cout<<(uint32_t)(*it);
	    if(it+1 != this->end())
	      std::cout<<", ";
	    else
	      std::cout<<" }";
	  }
	std::cout<<std::endl;
      }

    }
  }
}
#endif /* common_hpr_h */
