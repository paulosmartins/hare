#ifndef fv_common_h
#define fv_common_h

#include "tmp.hpp"
// #include "params.hpp"
#include <mm_malloc.h>
#include <gmpxx.h>
#include <chrono>
#include <nfl.hpp>

namespace rns
{
  namespace commons
  {
    namespace tools
    {

      size_t bit_rev(size_t x, size_t d)
      {
	size_t r = 0;
	for(size_t i = 0; i < d; ++i)
	  if(x&(1<<(d-i-1)))
	    r = r+(1<<i);
	return r;
      }

      template<class P>
      void CRT(mpz_class &r, const P &p, size_t d)
      {
	r = 0;
	mpz_class M = 1;
	mpz_class tmp;
	for(uint i = 0; i < P::nmoduli; i++)
	  {
	    mpz_set_ui(tmp.get_mpz_t(), p.get_modulus(i));
	    M = M*tmp;
	  }
	for(uint i = 0; i < P::nmoduli; i++)
	  {
	    mpz_set_ui(tmp.get_mpz_t(), p.get_modulus(i));
	    mpz_class tmp2 = M/tmp;
	    mpz_invert(tmp2.get_mpz_t(), tmp2.get_mpz_t(), tmp.get_mpz_t());
	    mpz_class tmp3;
	    mpz_set_ui(tmp3.get_mpz_t(), p(i, d));
	    r = r+((tmp2*tmp3)%tmp)*M/tmp;
	  }
	r = r%M;
      }

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
	nfl::ops::mulmod<P, nfl::simd::serial> op_mul;
	p.data()[cm*P::degree+i] = op_mul(p(cm, i), v, cm+P::indmod0);
      }

      // first to last: DEGREE coefficients, ie c_0,...,c_(deg-1)
      // put these into p, ie p <- (c_0,...,c_(deg-1)) mod q1, ..., (c_0,...,c_(deg-1)) mod qk
      // critical: used in homomorphic mult.
      template <class P, class It>
      inline void set_partial(P &p, It first, It last)
      {
	using main_type = typename P::value_type;
	using greater_main_type = typename P::greater_value_type;
	using main_signed_type = typename P::signed_value_type;

	auto* it = p.begin();
	for(uint cm = 0; cm < P::nmoduli; cm++)
	  {
	    main_type mod = p.get_modulus(cm);
	    for(It it2 = first; it2 < last; ++it2, ++it)
	      *it = (*it2 >= mod) ? *it2-mod : *it2;
	  }
      }

      // critical: used in homomorphic mult.
      template <class P>
      inline void spread_value(P &p, typename P::value_type v)
      {
	//for(auto* it = p.begin(); it < p.end(); ++it)
	//    *it = v;
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
      inline typename poly_type::value_type
      IMR(typename poly_type::greater_value_type in,
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
      inline typename poly_type::value_type
      IMR_full(typename poly_type::greater_value_type in,
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

#if AVX2
      // standard reduction in AVX2, with output in [0,mod)
      inline void modOpt_AVX2 (__m256i &u, const __m256i &mod, const __m256i &v0)
      {
	__m256i u1 = _mm256_srli_epi64(u, 32);
	__m256i q = _mm256_add_epi64(_mm256_mul_epu32(u1, v0), _mm256_slli_epi64(u, 2));
	q = _mm256_srli_epi64(q, 32);
	__m256i r = _mm256_sub_epi64(u, _mm256_mul_epu32(q, mod));
	r = _mm256_and_si256(r, parameters::mask32_avx);
	u = _mm256_min_epu32(r, _mm256_sub_epi32(r, mod));
      }
#endif

      // function SmMRq; cf. Alg. 2 in paper
      // variable acc contains the centered residues mod m_tilde
      // v contains precomputed data
      // computes p <- p + acc*v
      template<class P>
      void SmMRq(P &p, const std::array<typename P::signed_value_type, P::degree>
		 &acc, const std::array<typename P::value_type, P::nmoduli> &v)
      {
	using main_type = typename P::value_type;
	using main_greater_type = typename P::greater_value_type;
	using main_signed_type = typename P::signed_value_type;

	auto iter = p.begin();
	auto iter_v = v.begin();
	for(uint cm = 0; cm < P::nmoduli; cm++)
	  {
	    main_type mod = P::base.params.P[cm+P::indmod0];
	    main_greater_type w0 = (main_greater_type)P::base.params.Pn[cm+P::indmod0];
	    for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
	      {
		main_greater_type tmp = (*iter_acc >= 0) ? (main_greater_type)(*iter_v)*
		  (*iter_acc) : (main_greater_type)(*iter_v)*(-*iter_acc);
		main_type tmpv = IMR_full<P>(tmp, w0, mod);

		*iter += (*iter_acc >= 0) ? tmpv : mod-tmpv;
		if(*iter >= mod)
		  *iter -= mod;
		++iter;
	      }
	    ++iter_v;
	  }
      }
    }

    namespace base_conversions
    {

      template<class R, class E>
      void fast_div_round (R &r, const E &e,
			   const std::array<typename E::main_double, E::nmoduli> &thetas,
			   const std::array<typename R::value_type, R::nmoduli> &lambdas,
			   const std::array<std::array<typename R::value_type, E::nmoduli>, R::nmoduli>
			   &omegas)
      {
	using main_type = typename R::value_type;
	using main_greater_type = typename R::greater_value_type;
	using main_signed_type = typename R::signed_value_type;
	using main_double = typename R::main_double;
	using main_greater_double = typename R::main_greater_double;
	constexpr size_t degree = R::degree;

	constexpr size_t nmoduli = E::nmoduli;
	constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 * R::nbits);

	std::array<main_greater_type, degree> accm;
	std::array<main_double, degree> acc1;
	std::array<main_signed_type, degree> accs;
	acc1.fill(0.0);
	accs.fill(0);

	auto iter_omegas = omegas.begin();
	auto iter_thetas = thetas.begin();
	auto iter_lambdas = lambdas.begin();
	auto iter_r = r.begin();
	auto iter_acc1 = acc1.begin();
	auto iter_accs = accs.begin();
	auto iter_accm = accm.begin();

	for(uint cm = 0; cm < R::nmoduli; cm++, ++iter_lambdas, ++iter_omegas)
	  {
	    accm.fill(0);

	    auto iter_omegas_bis = iter_omegas->begin();
	    main_type mod_R = R::base.params.P[cm+R::indmod0];
	    main_greater_type w0_R = static_cast<main_greater_type>(R::base.params.Pn[cm
										      +R::indmod0]);

	    auto* iter_e = e.begin();
	    size_t cm_e = 0;
	    uint cpt = 0;

	    while(iter_e != e.end())
	      {
		if constexpr (nmoduli >= gap_lzr)
		    {
		      cpt++;
		    }

		if(cm == 0)
		  iter_acc1 = acc1.begin();

		for(iter_accm = accm.begin(); iter_accm != accm.end(); ++iter_accm, ++iter_e)
		  {
		    main_greater_type coeff_e = static_cast<main_greater_type>(*iter_e);

		    if(cm == 0)
		      {
			main_double coeff_e_dbl = static_cast<main_double>(coeff_e);
			if (coeff_e > (E::base.params.P[cm_e+E::indmod0]>>1))
			  coeff_e_dbl -= static_cast<main_greater_type>(E::base.params.P[cm_e
											 +E::indmod0]);
			*iter_acc1++ += coeff_e_dbl * (*iter_thetas);
		      }

		    if (coeff_e > (E::base.params.P[cm_e+E::indmod0]/2))
		      {
			main_type mpj = (E::base.params.P[cm_e+E::indmod0] > R::base.params.P[cm
											      +R::indmod0] ?
					 (R::base.params.P[cm+R::indmod0]<<1) - E::base.params.P[cm_e+E::indmod0] :
					 R::base.params.P[cm+R::indmod0] - E::base.params.P[cm_e+E::indmod0]);
			coeff_e += mpj;
			if (coeff_e >= R::base.params.P[cm+R::indmod0])
			  coeff_e -= R::base.params.P[cm+R::indmod0];
		      }
		    *iter_accm += coeff_e * (*iter_omegas_bis);
		  }

		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   {
				     for(auto iter_accm = accm.begin(); iter_accm != accm.end(); ++iter_accm)
				       *iter_accm = commons::reductions::IMR<R>(*iter_accm, w0_R, mod_R);
				     cpt = 0;
				   }
			       }

		if(cm == 0)
		  ++iter_thetas;
		++iter_omegas_bis;
		++cm_e;
	      }

	    if(cm == 0)
	      iter_acc1 = acc1.begin();
	    iter_accs = accs.begin();
	    iter_accm = accm.begin();
	    for(iter_accm = accm.begin(); iter_accm != accm.end();
		++iter_accs, ++iter_accm, ++iter_r)
	      {
		if(cm == 0)
		  *iter_accs = std::round(*iter_acc1++);

		main_greater_double tmp = *iter_accs;

		while (tmp < 0)
		  tmp += R::base.params.P[cm+R::indmod0];
		while (tmp >= R::base.params.P[cm+R::indmod0])
		  tmp -= R::base.params.P[cm+R::indmod0];

		*iter_accm += static_cast<main_greater_type>(tmp) * (*iter_lambdas);
		*iter_r += commons::reductions::IMR_full<R>(*iter_accm, w0_R, mod_R);

		if(*iter_r >= mod_R)
		  *iter_r -= mod_R;
	      }
	  }
      }

      /*template<class R, class E>
	void fast_div_round (R &r, const E &e,
	const std::array<typename R::main_double, E::nmoduli> &thetas,
	const std::array<typename R::value_type, R::nmoduli> &lambdas,
	const std::array<std::array<typename E::value_type, E::nmoduli>, R::nmoduli> &omegas)
	{
        using main_type = typename R::value_type;
        using main_greater_type = typename R::greater_value_type;
        using main_signed_type = typename R::signed_value_type;
        using main_double = typename R::main_double;
        using main_greater_double = typename R::main_greater_double;

        std::array<main_double, DEGREE> v_dbl;
        v_dbl.fill (0.0);
        std::array<main_signed_type, DEGREE> v;

        for (size_t i = 0; i < R::nmoduli; i++) {
	std::array<main_greater_type, DEGREE> acc;
	acc.fill (0);
	main_greater_type cpt = 0;

	for (size_t j = 0; j < E::nmoduli; j++) {
	cpt++;

	for (size_t k = 0; k < DEGREE; k++) {
	main_greater_type ejk = e (j, k);
	if (i == 0) {
	main_double d_ejk = e (j, k);
	if (ejk > E::base.params.P[j+E::indmod0]/2)
	d_ejk -= E::base.params.P[j+E::indmod0];

	v_dbl[k] += d_ejk * thetas[j];
	}

	if (ejk > E::base.params.P[j+E::indmod0] / 2) {

	main_type mpj = (E::base.params.P[j+E::indmod0] > R::base.params.P[i+R::indmod0] ?
	(R::base.params.P[i+R::indmod0]<<1) - E::base.params.P[j+E::indmod0] :
	R::base.params.P[i+R::indmod0] - E::base.params.P[j+E::indmod0]);

	ejk += mpj;

	if (ejk >= R::base.params.P[i+R::indmod0])
	ejk -= R::base.params.P[i+R::indmod0];
	}

	acc[k] += ejk * omegas[i][j];

	if (cpt == GAP_LZR) {
	acc[k] = reductions::IMR<R> (acc[k],
	R::base.params.Pn[i+R::indmod0],
	R::base.params.P[i+R::indmod0]);
	}
	}

	if (cpt == GAP_LZR)
	cpt = 0;
	}

	for (size_t j = 0; j < DEGREE; j++) {
	if (i == 0) {
	v[j] = std::round (v_dbl[j]);
	}

	acc[j] += static_cast<main_greater_type>(r (i, j)) * lambdas[i];

	r (i, j) = reductions::IMR_full<R> (acc[j],
	R::base.params.Pn[i+R::indmod0],
	R::base.params.P[i+R::indmod0]);


	main_greater_double vj_mod_pi = v[j];
	while (vj_mod_pi < 0)
	vj_mod_pi += R::base.params.P[i+R::indmod0];
	while (vj_mod_pi >= R::base.params.P[i+R::indmod0])
	vj_mod_pi -= R::base.params.P[i+R::indmod0];
	main_type vj = static_cast<main_type>(vj_mod_pi);

	r (i, j) += vj;

	if (r (i, j) >= R::base.params.P[i+R::indmod0])
	r (i, j) -= R::base.params.P[i+R::indmod0];
	}
        }
	}*/

      // fast base conversion from (1 poly with NB_MOD moduli) towards (1 poly with many moduli)
      // handles lazy reduction
      // variable coeffs contains precomputed data: matrix k1*k2
      // computes R(cm, d) <- R(cm, d) + sum_(cm2 = 1 to k2) E(cm2, d)*coeffs[cm][cm2] (ie matrix-vector product)
      template<class R, class E>
      void fast_1m_to_1m(R &r, const E &e,
			 const std::array<std::array<typename E::value_type, E::nmoduli>, R::nmoduli>
			 &coeffs)
      {
	using main_type = typename R::value_type;
	using main_greater_type = typename R::greater_value_type;
	using main_signed_type = typename R::signed_value_type;
	std::array<main_greater_type, R::degree> acc;

	constexpr size_t nmoduli = 2*E::nmoduli;
	constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 * R::nbits);

	auto iter_c = coeffs.begin();
	auto iter_r = r.begin();

	for(uint cm = 0; cm < R::nmoduli; cm++)
	  {
	    main_type mod = R::base.params.P[cm+R::indmod0];
	    main_greater_type w0 = static_cast<main_greater_type>(R::base.params.Pn[cm
										    +R::indmod0]);

	    acc.fill(0);
	    auto iter_c2 = iter_c->begin();
	    auto iter_e = e.begin();
	    uint cpt =
	      0; // counting the number of terms in the sum; goal: reduce after every GAP_LZR terms; only useful if NB_MOD >= GAP_LZR
	    while(iter_e < e.end())
	      {
		if constexpr (nmoduli >= gap_lzr)
		    {
		      ++cpt;
		    }

		for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
		  {
		    *iter_acc += (*iter_c2)*static_cast<main_greater_type>(*iter_e++);
		    if constexpr (nmoduli >= gap_lzr)
				   {
				     if(cpt == gap_lzr)
				       {
					 main_type tmp = reductions::IMR<R>(*iter_acc, w0, mod);
					 *iter_acc = static_cast<main_greater_type>(tmp);
				       }
				   }
		  }
		if constexpr (nmoduli >= gap_lzr)
			      {
				if(cpt == gap_lzr)
				  cpt = 0;
			      }
		++iter_c2;
	      }

	    for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_r)
	      {
		*iter_acc += static_cast<main_greater_type>(*iter_r);
		*iter_r = reductions::IMR_full<R>(*iter_acc, w0, mod);
	      }
	    ++iter_c;
	  }
      }


      template<class R, class E>
      void fast_1dm_to_1m(R &r, const E &e,
			  const std::array<std::array<typename E::value_type, E::nmoduli>, R::nmoduli>
			  &coeffs)
      {
	using main_type = typename R::value_type;
	using main_greater_type = typename R::greater_value_type;
	using main_signed_type = typename R::signed_value_type;

	constexpr size_t nmoduli = E::nmoduli;
	constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 * R::nbits);

	std::array<main_greater_type, R::degree> acc;

	auto iter_c = coeffs.begin();
	auto iter_r = r.begin();

	for(uint cm = 0; cm < R::nmoduli; cm++)
	  {
	    main_type mod = R::base.params.P[cm+R::indmod0];
	    main_greater_type w0 = static_cast<main_greater_type>(R::base.params.Pn[cm
										    +R::indmod0]);

	    acc.fill(0);
	    auto iter_c2 = iter_c->begin();
	    auto iter_e = e.begin();
	    uint cpt =
	      0; // counting the number of terms in the sum; goal: reduce after every GAP_LZR terms; only useful if NB_MOD >= GAP_LZR

	    while(iter_e < e.end())
	      {
		if constexpr (nmoduli >= gap_lzr)
		    {
		      ++cpt;
		    }

		for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
		  {
		    *iter_acc += (*iter_c2)*static_cast<main_greater_type>(*iter_e++);

		    if constexpr (nmoduli >= gap_lzr)
				   {
				     if(cpt == gap_lzr)
				       {
					 main_type tmp = reductions::IMR<R>(*iter_acc, w0, mod);
					 *iter_acc = static_cast<main_greater_type>(tmp);
				       }
				   }
		  }
		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   cpt = 0;
			       }
		++iter_c2;
	      }

	    for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_r)
	      {
		*iter_acc += static_cast<main_greater_type>(*iter_r);
		*iter_r = reductions::IMR_full<R>(*iter_acc, w0, mod);
	      }
	    ++iter_c;
	  }
      }

      template<class Rm, class E>
      void fast_1m_to_1m_float (Rm &rm, const E &e,
				const std::array<typename E::main_double, E::nmoduli> &q_i_inv,
				const std::array<std::array<typename E::value_type, E::nmoduli>, Rm::nmoduli>
				&q_i,
				const std::array<typename Rm::value_type, Rm::nmoduli> &q)
      {
	using main_type = typename Rm::value_type;
	using main_greater_type = typename Rm::greater_value_type;
	using main_signed_type = typename Rm::signed_value_type;
	using main_double = typename Rm::main_double;
	using main_greater_double = typename Rm::main_greater_double;

	constexpr size_t nmoduli = E::nmoduli;
	constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					  Rm::nbits);
	constexpr size_t degree = Rm::degree;

	std::array<main_greater_type, degree> accm;
	std::array<main_double, degree> acc1;
	std::array<main_signed_type, degree> accs;
	acc1.fill(0.0);
	accs.fill(0);

	auto iter_c_E_Rm = q_i.begin();
	auto iter_c_E_R1 = q_i_inv.begin();
	auto iter_q = q.begin();
	auto iter_rm = rm.begin();
	auto iter_acc1 = acc1.begin();
	auto iter_accs = accs.begin();
	auto iter_accm = accm.begin();

	for(uint cm = 0; cm < Rm::nmoduli; cm++, ++iter_q, ++iter_c_E_Rm)
	  {
	    accm.fill(0);

	    auto iter_c_E_Rm_bis = iter_c_E_Rm->begin();
	    main_type mod_Rm = Rm::base.params.P[cm+Rm::indmod0];
	    main_greater_type w0_Rm = static_cast<main_greater_type>
	      (Rm::base.params.Pn[cm+Rm::indmod0]);

	    auto iter_e = e.begin();
	    size_t cm_e = 0;
	    uint cpt = 0;
	    while(iter_e != e.end())
	      {
		if constexpr (nmoduli >= gap_lzr)
			       {
				 cpt++;
			       }
		if(cm == 0)
		  iter_acc1 = acc1.begin();

		for(iter_accm = accm.begin(); iter_accm != accm.end(); ++iter_accm, ++iter_e)
		  {
		    main_greater_type coeff_e = static_cast<main_greater_type>(*iter_e);

		    if(cm == 0)
		      {
			main_double coeff_e_dbl = static_cast<main_double>(coeff_e);
			if (coeff_e > (E::base.params.P[cm_e+E::indmod0]/2))
			  coeff_e_dbl -= static_cast<main_greater_type>(E::base.params.P[cm_e
											 +E::indmod0]);
			*iter_acc1++ += coeff_e_dbl * (*iter_c_E_R1);
		      }

		    if (coeff_e > (E::base.params.P[cm_e+E::indmod0]/2))
		      {
			main_type mpj = (E::base.params.P[cm_e+E::indmod0] > Rm::base.params.P[cm
											       +Rm::indmod0] ?
					 (Rm::base.params.P[cm+Rm::indmod0]<<1) - E::base.params.P[cm_e+E::indmod0] :
					 Rm::base.params.P[cm+Rm::indmod0] - E::base.params.P[cm_e+E::indmod0]);
			coeff_e += mpj;
			if (coeff_e >= Rm::base.params.P[cm+Rm::indmod0])
			  coeff_e -= Rm::base.params.P[cm+Rm::indmod0];
		      }
		    *iter_accm += coeff_e * (*iter_c_E_Rm_bis);
		  }

		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   {
				     for(auto iter_accm = accm.begin(); iter_accm != accm.end(); ++iter_accm)
				       *iter_accm = commons::reductions::IMR<Rm>(*iter_accm, w0_Rm, mod_Rm);
				     cpt = 0;
				   }
			       }

		if(cm == 0)
		  ++iter_c_E_R1;
		++iter_c_E_Rm_bis;
		++cm_e;
	      }
	    if(cm == 0)
	      iter_acc1 = acc1.begin();
	    iter_accs = accs.begin();
	    for(iter_accm = accm.begin(); iter_accm != accm.end();
		++iter_accs, ++iter_accm, ++iter_rm)
	      {
		if(cm == 0)
		  *iter_accs = -std::round(*iter_acc1++);
		if(*iter_accs < 0)
		  *iter_accm += static_cast<main_greater_type>((*iter_accs) + mod_Rm)  *
		    (*iter_q);
		else
		  *iter_accm += static_cast<main_greater_type>(*iter_accs) * (*iter_q);
		*iter_rm = commons::reductions::IMR_full<Rm>(*iter_accm, w0_Rm, mod_Rm);
	      }
	  }
      }


      /*template<class Rm, class E>
	void fast_1m_to_1m_float (Rm &rm, const E &e,
	const std::array<typename E::main_double, E::nmoduli> &q_i_inv,
	const std::array<std::array<typename E::value_type, E::nmoduli>, Rm::nmoduli> &q_i,
	const std::array<typename Rm::value_type, Rm::nmoduli> &q)
	{
        using main_type = typename Rm::value_type;
        using main_greater_type = typename Rm::greater_value_type;
        using main_signed_type = typename Rm::signed_value_type;
        using main_double = typename Rm::main_double;
        using main_greater_double = typename Rm::main_greater_double;

        std::array<main_double, DEGREE> v_dbl;
        v_dbl.fill (0.0);
        std::array<main_signed_type, DEGREE> v;

        for (size_t i = 0; i < Rm::nmoduli; i++) {
	std::array<main_greater_type, DEGREE> acc;
	acc.fill (0);
	main_greater_type cpt = 0;

	for (size_t j = 0; j < E::nmoduli; j++) {
	cpt++;

	for (size_t k = 0; k < DEGREE; k++) {
	main_greater_type ejk = e (j, k);

	if (i == 0) {
	main_double d_ejk = e (j, k);

	if (ejk > E::base.params.P[j+E::indmod0] / 2)
	d_ejk -= E::base.params.P[j+E::indmod0];

	v_dbl[k] += d_ejk * q_i_inv[j];
	}

	if (ejk > E::base.params.P[j+E::indmod0] / 2) {

	main_type mpj = (E::base.params.P[j+E::indmod0] > Rm::base.params.P[i+Rm::indmod0] ?
	(Rm::base.params.P[i+Rm::indmod0]<<1) - E::base.params.P[j+E::indmod0] :
	Rm::base.params.P[i+Rm::indmod0] - E::base.params.P[j+E::indmod0]);
	ejk += mpj;

	if (ejk >= Rm::base.params.P[i+Rm::indmod0])
	ejk -= Rm::base.params.P[i+Rm::indmod0];
	}

	acc[k] += ejk * q_i[i][j];
	if (cpt == GAP_LZR) {
	acc[k] = reductions::IMR<Rm> (acc[k],
	Rm::base.params.Pn[i+Rm::indmod0],
	Rm::base.params.P[i+Rm::indmod0]);
	}
	}

	if (cpt == GAP_LZR)
	cpt = 0;

	}

	for (size_t j = 0; j < DEGREE; j++) {
	if (i == 0) {
	v[j] = -std::round (v_dbl[j]);
	}

	acc[j] += static_cast<main_greater_type>
	((v[j] < 0 ?
	v[j] + Rm::base.params.P[i+Rm::indmod0] :
	v[j])) *
	q[i];

	rm (i, j) = reductions::IMR_full<Rm> (acc[j],
	Rm::base.params.Pn[i+Rm::indmod0],
	Rm::base.params.P[i+Rm::indmod0]);
	}
        }
	}*/
      // fast base conversion from (1 poly+NB_MOD moduli: variable e) towards (1 poly+many moduli: variable rm, and 1 poly+1 modulus: variable r1)
      template<class Rm, class R1, class E>
      void fast_1m_to_1m_11(Rm &rm, R1 &r1, const E &e,
			    const std::array<std::array<typename E::value_type, E::nmoduli>, Rm::nmoduli>
			    &coeffs_E_Rm, const std::array<typename E::value_type, E::nmoduli> &coeffs_E_R1)
      {
	using main_type = typename Rm::value_type;
	using main_greater_type = typename Rm::greater_value_type;
	using main_signed_type = typename Rm::signed_value_type;
	std::array<main_greater_type, Rm::degree> acc, acc2;

	constexpr size_t nmoduli = E::nmoduli;
	constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					  Rm::nbits);

	acc2.fill(0);

	auto iter_c_E_Rm = coeffs_E_Rm.begin();
	auto iter_c_E_R1 = coeffs_E_R1.begin();
	auto iter_rm = rm.begin();
	auto iter_acc2 = std::begin(acc2);

	main_type mod_R1 = R1::base.params.P[R1::indmod0];
	main_greater_type w0_R1 = static_cast<main_greater_type>
	  (R1::base.params.Pn[R1::indmod0]);

	for(uint cm = 0; cm < Rm::nmoduli; cm++)
	  {
	    acc.fill(0);
	    auto iter_c_E_Rm_bis = iter_c_E_Rm->begin();
	    main_type mod_Rm = Rm::base.params.P[cm+Rm::indmod0];
	    main_greater_type w0_Rm = static_cast<main_greater_type>
	      (Rm::base.params.Pn[cm+Rm::indmod0]);

	    auto* iter_e = e.begin();
	    uint cpt = 0;
	    while(iter_e < e.end())
	      {
		if(cm == 0)
		  iter_acc2 = acc2.begin();
		if constexpr (nmoduli >= gap_lzr)
		    {
		      ++cpt;
		    }

		for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
		  {
		    main_greater_type tmp = static_cast<main_greater_type>(*iter_e++);
		    *iter_acc += (*iter_c_E_Rm_bis)*tmp;
		    if(cm == 0)
		      *iter_acc2++ += (*iter_c_E_R1)*tmp;
		    if constexpr (nmoduli >= gap_lzr)
				   {
				     if(cpt == gap_lzr)
				       {
					 main_type tmp = reductions::IMR<Rm>(*iter_acc, w0_Rm, mod_Rm);
					 *iter_acc = static_cast<main_greater_type>(tmp);
					 if(cm == 0)
					   {
					     --iter_acc2;
					     main_type tmp = reductions::IMR<R1>(*iter_acc2, w0_R1, mod_R1);
					     *iter_acc2++ = static_cast<main_greater_type>(tmp);
					   }
				       }
				   }
		  }
		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   cpt = 0;
			       }

		++iter_c_E_Rm_bis;
		if(cm == 0)
		  ++iter_c_E_R1;
	      }

	    for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_rm)
	      {
		*iter_acc += static_cast<main_greater_type>(*iter_rm);
		*iter_rm = reductions::IMR_full<Rm>(*iter_acc, w0_Rm, mod_Rm);
	      }
	    if(cm == 0)
	      {
		auto iter_r1 = r1.begin();
		for(auto iter_acc2 = acc2.begin(); iter_acc2 < acc2.end();
		    ++iter_acc2, ++iter_r1)
		  {
		    *iter_acc2 += static_cast<main_greater_type>(*iter_r1);
		    *iter_r1 = reductions::IMR_full<R1>(*iter_acc2, w0_R1, mod_R1);
		  }
	      }
	    ++iter_c_E_Rm;
	  }
      }

      // fast base conversion from (2 poly with NB_MOD moduli: e1 and e2) towards (1 poly+1 modulus: r)
      template<class R, class Em1, class Em2>
      void fast_2m_to_11(R &r, const Em1 &e1, const Em2 &e2,
			 const std::array<typename Em1::value_type, Em1::nmoduli> &coeffs_Em1_R,
			 const std::array<typename Em2::value_type, Em2::nmoduli> &coeffs_Em2_R)
      {
	using main_type = typename R::value_type;
	using main_greater_type = typename R::greater_value_type;
	using main_signed_type = typename R::signed_value_type;

	std::array<main_greater_type, R::degree> acc;
	acc.fill(0);

	constexpr size_t nmoduli = Em1::nmoduli + Em2::nmoduli;
	constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 * R::nbits);

	main_type mod_R = R::base.params.P[R::indmod0];
	main_greater_type w0_R = static_cast<main_greater_type>
	  (R::base.params.Pn[R::indmod0]);

	uint cpt = 0;

	auto iter_c_Em1_R = coeffs_Em1_R.begin();
	auto iter_e1 = e1.begin();
	while(iter_e1 < e1.end())
	  {
	    if constexpr (nmoduli >= gap_lzr)
		{
		  ++cpt;
		}

	    main_greater_type tmpv = static_cast<main_greater_type>(*iter_c_Em1_R++);
	    for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
	      {
		*iter_acc += tmpv*static_cast<main_greater_type>(*iter_e1++);
		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   {
				     main_type tmp = reductions::IMR<R>(*iter_acc, w0_R, mod_R);
				     *iter_acc = static_cast<main_greater_type>(tmp);
				   }
			       }
	      }
	    if constexpr (nmoduli >= gap_lzr)
			   {
			     if(cpt == gap_lzr)
			       cpt = 0;
			   }
	  }

	auto iter_c_Em2_R = coeffs_Em2_R.begin();
	auto iter_e2 = e2.begin();
	while(iter_e2 < e2.end())
	  {
	    if constexpr (nmoduli  >= gap_lzr)
		{
		  ++cpt;
		}


	    main_greater_type tmpv = static_cast<main_greater_type>(*iter_c_Em2_R++);
	    for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
	      {
		*iter_acc += tmpv*static_cast<main_greater_type>(*iter_e2++);
		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   {
				     main_type tmp = reductions::IMR<R>(*iter_acc, w0_R, mod_R);
				     *iter_acc = static_cast<main_greater_type>(tmp);
				   }
			       }
	      }
	    if constexpr (nmoduli >= gap_lzr)
			   {
			     if(cpt == gap_lzr)
			       cpt = 0;
			   }
	  }

	auto iter_r = r.begin();
	for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_r)
	  {
	    *iter_acc += static_cast<main_greater_type>(*iter_r);
	    *iter_r = reductions::IMR_full<R>(*iter_acc, w0_R, mod_R);
	  }
      }

      // fast conversion from e to r
      // and Shenoy and Kumaresan like base conversion from esk to r (ie centered residue modulo m_sk)
      template<class R, class E, class Esk>
      void sk(R &r, const E &e, const Esk &esk,
	      const std::array<std::array<typename E::value_type, E::nmoduli>, R::nmoduli>
	      &coeffs_E_R, const std::array<typename R::value_type, R::nmoduli> &coeffs_Esk_R)
      {
	using main_type = typename R::value_type;
	using main_greater_type = typename R::greater_value_type;
	using main_signed_type = typename R::signed_value_type;
	std::array<main_greater_type, R::degree> acc;

	constexpr size_t nmoduli = E::nmoduli+1;
	constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 * R::nbits);

	auto iter_r = r.begin();
	auto iter_c_E_R = coeffs_E_R.begin();
	auto iter_c_Esk_R = coeffs_Esk_R.begin();

	main_type m_sk = Esk::base.params.P[Esk::indmod0];
	main_type m_sk_ov2 = (m_sk-1)>>1;

	for(uint cm = 0; cm < R::nmoduli; cm++)
	  {
	    main_type mod_R = R::base.params.P[cm+R::indmod0];
	    main_greater_type w0_R = static_cast<main_greater_type>(R::base.params.Pn[cm
										      +R::indmod0]);

	    acc.fill(0);
	    uint cpt = 0;
	    auto iter_c_E_R2 = iter_c_E_R->begin();
	    auto iter_e = e.begin();
	    while(iter_e < e.end())
	      {
		if constexpr (nmoduli >= gap_lzr)
		    {
		      ++cpt;
		    }
		
		for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
		  {
		    *iter_acc += (*iter_c_E_R2)*static_cast<main_greater_type>(*iter_e++);
		    if constexpr (nmoduli >= gap_lzr)
				   {
				     if(cpt == gap_lzr)
				       {
					 main_type tmp = reductions::IMR<R>(*iter_acc, w0_R, mod_R);
					 *iter_acc = static_cast<main_greater_type>(tmp);
				       }
				   }
		  }
		
		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   cpt = 0;
			       }

		++iter_c_E_R2;
	      }
	    ++iter_c_E_R;
	    auto iter_esk = esk.begin();
	    while(iter_esk < esk.end())
	      {
		if constexpr (nmoduli >= gap_lzr)
		    {
		      ++cpt;
		    }

		for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
		  {
		    // centering the residue mod m_sk
		    *iter_acc += (*iter_esk <= m_sk_ov2) ? (*iter_c_Esk_R)
		      *static_cast<main_greater_type>(*iter_esk) : (*iter_c_Esk_R)
		      *static_cast<main_greater_type>(mod_R+*iter_esk-m_sk);
		    if constexpr (nmoduli >= gap_lzr)
				   {
				     if(cpt == gap_lzr)
				       {
					 main_type tmp = reductions::IMR<R>(*iter_acc, w0_R, mod_R);
					 *iter_acc = static_cast<main_greater_type>(tmp);
				       }
				   }
		    ++iter_esk;
		  }
		if constexpr (nmoduli >= gap_lzr)
			       {
				 if(cpt == gap_lzr)
				   cpt = 0;
			       }

		++iter_c_Esk_R;
	      }

	    for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_r)
	      {
		*iter_acc += static_cast<main_greater_type>(*iter_r);
		*iter_r = reductions::IMR_full<R>(*iter_acc, w0_R, mod_R);
	      }
	  }
      }
    }

    namespace key_gen
    {
      template<class P, class parameters>
      void key_gen(P* keys, int h = 0)
      {
	using main_type = typename P::value_type;
	using main_greater_type = typename P::greater_value_type;
	using main_signed_type = typename P::signed_value_type;

	// creates keys in ntt mode
	P* a = &keys[0];
	P* b = &keys[1];
	P* s = &keys[2];

	nfl::uniform unif;
	nfl::non_uniform nunif(parameters::B_key+1);
	P* e = mem_manag::alloc_aligned<P, 32, nfl::gaussian<uint8_t, main_type, 2>>(1,
										     nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng));
	*s = nunif;
	*a = unif;

	e[0].ntt_pow_phi();
	s->ntt_pow_phi();
	a->ntt_pow_phi();
	*b = *a**s+e[0];
	e[0] = 0;
	*b = e[0]-*b;

	mem_manag::free_aligned<P>(1, e);
      }
    }

    namespace fast_pol
    {
      // struct for polynomials modulo m_tilde (NTT not required: no need for polynomial product here)
      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      struct PolMT
      {
	using iterator = type_mtilde*;
	using const_iterator = type_mtilde const*;
	using type_mtilde_ = type_mtilde;
	using signed_type_mtilde_ = signed_type_mtilde;

	static constexpr type_mtilde mtilde_half = static_cast<type_mtilde>(1)<<
          (bs_mt-1);
	static constexpr main_type mtilde = static_cast<main_type>(1)<<bs_mt;

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

	bool operator==(const PolMT &other) const;

	void set(const PolMT &pol);
	void set0()
	{
	  std::fill(__data, __data + degree, 0);
	}

	// base conversion
	template <class P>
	void fast_base_conversion(const P &p,
				  const std::array<type_mtilde, P::nmoduli> &v);

	// centered remainders
	inline void get_modc(std::array<main_signed_type, degree> &v) const;

	// print for debug
	void print() const;
      };

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      constexpr type_mtilde
      PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::mtilde_half;
      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      constexpr main_type
      PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::mtilde;

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      void PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::set(
												   const PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>
												   &pol)
      {
	auto iter = pol.begin();
	for(auto it = begin(); it < end(); ++it, ++iter)
	  *it = static_cast<type_mtilde>(*iter);
      }

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      bool PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::operator==
      (const PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>
       &other) const
      {
	auto iter_ot = other.begin();
	for(auto it = begin(); it < end(); ++it, ++iter_ot)
	  {
	    if (*it != *iter_ot)
	      return false;
	  }
	return true;
      }

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      template<class P>
      void PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::fast_base_conversion(
														    const P &p, const std::array<type_mtilde, P::nmoduli> &v)
      {
	auto citer = p.begin();
	for(auto itv = v.begin(); itv < v.end(); ++itv)
	  for(auto it = begin(); it < end(); ++it, ++citer)
	    *it += (*itv)*static_cast<type_mtilde>(*citer);
      }

      template<typename main_type, typename main_signed_type, typename type_mtilde, typename signed_type_mtilde, size_t bs_mt, size_t degree>
      inline void
      PolMT<main_type, main_signed_type, type_mtilde, signed_type_mtilde, bs_mt, degree>::get_modc(
												   std::array<main_signed_type, degree> &v) const
      {
	// if residue > mod/2, return residue - mod
	// handled through simple conversion to signed type (except for mtilde_half)
	auto iter = begin();
	for(auto iterv = v.begin(); iterv < v.end(); ++iterv)
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
	using iterator = type_t*;
	using const_iterator = type_t const*;

	type_t* __data; // contains polynomial coeff.

	PolT()
	{
	  __data = (type_t*)_mm_malloc(degree*sizeof(type_t),32);
	  std::fill(__data, __data + degree, 0);
	};

	~PolT()
	{
	  _mm_free(__data);
	}

	PolT(const PolT &other)
	{
	  __data = (type_t*)_mm_malloc(degree*sizeof(type_t),32);
	  std::copy(&other.__data[0], &other.__data[0]+degree, &__data[0]);
	}

	PolT& operator=(const PolT &other)
	{
	  std::copy(&other.__data[0], &other.__data[0]+degree, &__data[0]);
	  return *this;
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

	// setters
	void set0()
	{
	  std::fill(__data, __data + degree, 0);
	};
	template <class P>
	void set(const P &p);
	template <class P>
	void get(P &p) const;

	// ope
	bool operator==(const PolT &other) const;
	void mul(const PolT &other);
	void modt();

	// print for debug
	void print() const;
      };

      // fill current PolT with p
      template<typename parameters, size_t degree>
      template <class P>
      void PolT<parameters, degree>::get(P &p) const
      {
	for(uint i = 0; i < degree; i++)
	  for(uint cm = 0; cm < P::nmoduli; cm++)
	    *(p.begin()+(cm<<static_log2<degree>::value)+i) =
	      static_cast<typename P::value_type>(__data[i]);
      }

      // fill p with current PolT
      template<typename parameters, size_t degree>
      template<class P>
      void PolT<parameters, degree>::set(const P &p)
      {
	auto iter = p.begin();
	typename P::value_type mod = p.get_modulus(0);
	for(auto it = begin(); it < end(); ++it, ++iter)
	  *it = (*iter < (mod>>1)) ? static_cast<type_t>(*iter)&parameters::mask_t :
	    parameters::t-static_cast<type_t>(mod-*iter)&parameters::mask_t;
      }

      template<typename parameters, size_t degree>
      bool PolT<parameters, degree>::operator==(const PolT<parameters, degree> &other)
	const
      {
	auto iter_ot = other.begin();
	for(auto it = begin(); it < end(); ++it, ++iter_ot)
	  {
	    if (*it != *iter_ot)
	      return false;
	  }
	return true;
      }

      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::modt()
      {
	for(auto it = begin(); it < end(); ++it)
	  *it = (*it)&parameters::mask_t;
      }

      // quadratic multiplication: only used for checking
      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::mul(const PolT<parameters, degree> &other)
      {
	std::array<type_t, degree> tmp;
	tmp.fill(0);
	for(uint k = 0; k <= degree-1; k++)
	  for(uint i = 0; i <= k; i++)
	    tmp[k] += __data[i] * other.__data[k-i];
	for(uint k = degree; k <= 2*degree-1; k++)
	  for(uint i = k-degree+1; i <= degree-1; i++)
	    tmp[k-degree] -= __data[i]*other.__data[k-i];
	for(uint k = 0; k <= degree-1; k++)
	  __data[k] = tmp[k];
      }

      template<typename parameters, size_t degree>
      void PolT<parameters, degree>::print() const
      {
	std::cout<<"{ ";
	for(uint i = 0; i < degree; ++i)
	  {
	    std::cout<<(uint32_t)(__data[i]);
	    if(i<degree-1)
	      std::cout<<", ";
	    else
	      std::cout<<" }";
	  }
	std::cout<<std::endl;
      }

    }

  }
}
#endif
