#ifndef fvrns_be_hpr_h
#define fvrns_be_hpr_h

#include "tmp.hpp"
// #include "params.hpp"
#include "hpr/common_hpr.hpp"

namespace hpr
{
  namespace base_conversions
  {
    template<class Rm, class E>
    void fast_1m_to_1m_float (Rm &rm, const E &e,
			      const std::array<typename Rm::main_double, E::nmoduli> &q_i_inv,
			      const std::array<std::array<typename E::value_type, E::nmoduli>, Rm::nmoduli>
			      &q_i,
			      const std::array<typename Rm::value_type, Rm::nmoduli> &q)
    {
      using main_type = typename Rm::value_type;
      using main_greater_type = typename Rm::greater_value_type;
      using main_signed_type = typename Rm::signed_value_type;
      using main_double = typename Rm::main_double;
      using main_greater_double = typename Rm::main_greater_double;
      constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					Rm::nbits);
      constexpr size_t degree = Rm::degree;

      std::array<main_double, degree> v_dbl;
      v_dbl.fill (0.0);
      std::array<main_signed_type, degree> v;

      for (size_t i = 0; i < Rm::nmoduli; i++)
	{
	  std::array<main_greater_type, degree> acc;
	  std::array<main_greater_type, degree> acc_sk;
	  acc.fill (0);
	  if (i == 0) acc_sk.fill (0);
	  main_greater_type cpt = 0;

	  for (size_t j = 0; j < E::nmoduli; j++)
	    {
	      cpt++;

	      for (size_t k = 0; k < degree; k++)
		{
		  main_greater_type ejk = e (j, k);

		  if (i == 0)
		    {
		      main_double d_ejk = e (j, k);
		      if (ejk > E::base.params.P[j+E::indmod0]/2)
			d_ejk -= E::base.params.P[j+E::indmod0];

		      v_dbl[k] += d_ejk * q_i_inv[j];
		    }

		  if (ejk > E::base.params.P[j+E::indmod0] / 2)
		    {

		      main_type mpj = (E::base.params.P[j+E::indmod0] > Rm::base.params.P[i
											  +Rm::indmod0] ?
				       (Rm::base.params.P[i+Rm::indmod0]<<1) - E::base.params.P[j+E::indmod0] :
				       Rm::base.params.P[i+Rm::indmod0] - E::base.params.P[j+E::indmod0]);
		      ejk += mpj;

		      if (ejk >= Rm::base.params.P[i+Rm::indmod0])
			ejk -= Rm::base.params.P[i+Rm::indmod0];
		    }

		  acc[k] += ejk * q_i[i][j];
		  if (cpt == gap_lzr)
		    {
		      acc[k] = commons::reductions::IMR<Rm> (acc[k],
							     Rm::base.params.Pn[i+Rm::indmod0],
							     Rm::base.params.P[i+Rm::indmod0]);
		    }
		}

	      if (cpt == gap_lzr)
		cpt = 0;

	    }

	  if (i == 0)
	    {
	      for (size_t j = 0; j < degree; j++)
		{
		  v[j] = -std::round (v_dbl[j]);
		}
	    }

	  for (size_t j = 0; j < degree; j++)
	    {
	      acc[j] += static_cast<main_greater_type>
		((v[j] < 0 ?
		  v[j] + Rm::base.params.P[i+Rm::indmod0] :
		  v[j])) *
		q[i];

	      rm (i, j) = commons::reductions::IMR_full<Rm> (acc[j],
							     Rm::base.params.Pn[i+Rm::indmod0],
							     Rm::base.params.P[i+Rm::indmod0]);
	    }
	}
    }

    template<class Rm, class R1, class E>
    void fast_1m_to_1m_1s_float (Rm &rm, R1 &r1, const E &e,
				 const std::array<typename Rm::main_double, E::nmoduli> &q_i_inv,
				 const std::array<std::array<typename E::value_type, E::nmoduli>, Rm::nmoduli>
				 &q_i,
				 const std::array<typename E::value_type, E::nmoduli> &q_i_sk,
				 const std::array<typename E::value_type, Rm::nmoduli> &q,
				 const typename Rm::value_type &q_sk)
    {
      using main_type = typename Rm::value_type;
      using main_greater_type = typename Rm::greater_value_type;
      using main_signed_type = typename Rm::signed_value_type;
      using main_double = typename Rm::main_double;
      using main_greater_double = typename Rm::main_greater_double;
      constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					Rm::nbits);
      constexpr size_t degree = Rm::degree;

      std::array<main_double, degree> v_dbl;
      std::array<main_signed_type, degree> v;
      v_dbl.fill (0.0);

      for (size_t i = 0; i < Rm::nmoduli; i++)
	{
	  std::array<main_greater_type, degree> acc;
	  acc.fill (0);
	  std::array<main_greater_type, degree> acc_sk;
	  acc_sk.fill (0);
	  main_greater_type cpt = 0;

	  for (size_t j = 0; j < E::nmoduli; j++)
	    {
	      cpt++;

	      for (size_t k = 0; k < degree; k++)
		{
		  main_greater_type ejk = e (j, k);

		  if (i == 0)
		    {
		      main_double d_ejk = e (j, k);
		      if (ejk > E::base.params.P[j+E::indmod0]/2)
			d_ejk -= E::base.params.P[j+E::indmod0];

		      v_dbl[k] += d_ejk * q_i_inv[j];
		    }

		  if (ejk > E::base.params.P[j+E::indmod0] / 2)
		    {

		      main_type mpj = (E::base.params.P[j+E::indmod0] > Rm::base.params.P[i
											  +Rm::indmod0] ?
				       (Rm::base.params.P[i+Rm::indmod0]<<1) - E::base.params.P[j+E::indmod0] :
				       Rm::base.params.P[i+Rm::indmod0] - E::base.params.P[j+E::indmod0]);
		      ejk += mpj;

		      if (ejk >= Rm::base.params.P[i+Rm::indmod0])
			ejk -= Rm::base.params.P[i+Rm::indmod0];
		    }

		  acc[k] += ejk * q_i[i][j];
		  if (cpt == gap_lzr)
		    {
		      acc[k] = commons::reductions::IMR<Rm> (acc[k],
							     Rm::base.params.Pn[i+Rm::indmod0],
							     Rm::base.params.P[i+Rm::indmod0]);
		    }

		  if (i == 0)
		    {
		      main_greater_type ejk = e (j, k);
		      if (ejk > E::base.params.P[j+E::indmod0] / 2)
			{

			  main_type mpj = (E::base.params.P[j+E::indmod0] > R1::base.params.P[R1::indmod0]
					   ?
					   (R1::base.params.P[R1::indmod0]<<1) - E::base.params.P[j+E::indmod0] :
					   R1::base.params.P[R1::indmod0] - E::base.params.P[j+E::indmod0]);
			  ejk += mpj;

			  if (ejk >= R1::base.params.P[R1::indmod0])
			    ejk -= R1::base.params.P[R1::indmod0];
			}

		      acc_sk[k] += ejk * q_i_sk[j];
		      if (cpt == gap_lzr)
			{
			  acc_sk[k] = commons::reductions::IMR<R1> (acc_sk[k],
								    R1::base.params.Pn[R1::indmod0],
								    R1::base.params.P[R1::indmod0]);
			}
		    }
		}

	      if (cpt == gap_lzr)
		cpt = 0;

	    }

	  if (i == 0)
	    {
	      for (size_t j = 0; j < degree; j++)
		{
		  v[j] = -std::round (v_dbl[j]);
		}
	    }

	  for (size_t j = 0; j < degree; j++)
	    {
	      acc[j] += static_cast<main_greater_type>
		((v[j] < 0 ?
		  v[j] + Rm::base.params.P[i+Rm::indmod0] :
		  v[j])) *
		q[i];

	      rm (i, j) = commons::reductions::IMR_full<Rm> (acc[j],
							     Rm::base.params.Pn[i+Rm::indmod0],
							     Rm::base.params.P[i+Rm::indmod0]);

	      if (i == 0)
		{
		  acc_sk[j] += static_cast<main_greater_type>
		    ((v[j] < 0 ?
		      v[j] + R1::base.params.P[R1::indmod0] :
		      v[j])) *
		    q_sk;

		  r1 (0, j) = commons::reductions::IMR_full<R1> (acc_sk[j],
								 R1::base.params.Pn[R1::indmod0],
								 R1::base.params.P[R1::indmod0]);
		}
	    }
	}
    }


    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    struct base_conv_tool_set
    {
      static constexpr size_t nmQ = Pq::nmoduli;
      static constexpr size_t nmB = Pb::nmoduli;
      using main_type = typename Pq::value_type;
      using main_greater_type = typename Pq::greater_value_type;
      using main_signed_type = typename Pq::signed_value_type;
      using main_double = typename Pq::main_double;
      using main_greater_double = typename Pq::main_greater_double;

      typedef std::array<std::array<main_type, nmQ>, nmB> matBEq2b;
      typedef std::array<main_type, nmQ> vecBEq2sk;

      typedef std::array<main_type, nmB> vecBEb2sk;

      typedef std::array<std::array<main_type, nmB>, nmQ> matBEb2q;
      typedef std::array<main_type, nmQ> vecBEsk2q;

      typedef std::array<main_type, nmB+1> vecBE2bsk_mt;

      matBEq2b _Mq2b; // (q/q1,...,q/qk) mod (b1,...,bl)
      matBEb2q _Mb2q; // (b/b1,...,b/bl) mod (q1,...,qk)
      vecBEq2sk _Vq2sk; // (q/q1,...,q/qk) mod msk
      vecBEb2sk _Vb2sk; // (b/b1,...,b/bl) mod msk
      vecBEsk2q _Vsk2q; // (-b,...,-b) mod (q1,...,qk)

      matBEq2b _Mq2b_mt; // (q/q1/mt,...,q/qk/mt) mod (b1,...,bl)
      vecBEq2sk _Vq2sk_mt; // (q/q1/mt,...,q/qk/mt) mod msk
      vecBE2bsk_mt _V2bsk_mt; // (q/mt,...,q/mt, q/mt) mod (b1,...,bl, msk)

      Pq* _xi_q; // (q1/q mod q1, ....), (q/q1 mod q1, ....)
      Pb* _xi_b; // (b1/b mod b1, ....), (b/b1 mod b1, ....)
      Psk* _xi_sk; // -1/b mod msk
      Pq* _xi_q_gt; // (q1/q*gamma*t mod q1, ....)
      Pq* _xi_q_mt; // (mt*q1/q mod q1, ...), (q/q1/mt mod q1, ...)


      // for carry propag.
      Pb* _qinv_modB;
      Psk* _qinv_modmsk;

      void allocate();
      void copy_from(const base_conv_tool_set &other);
      base_conv_tool_set();
      base_conv_tool_set(const base_conv_tool_set &other);
      base_conv_tool_set& operator=(const base_conv_tool_set &other);
      ~base_conv_tool_set();

      void generate();

      void print_Mq2b();
      void print_Mb2q();
      void print_Vq2sk();
      void print_Vb2sk();
      void print_Vsk2q();

      template<class Rm, class R1, class E, typename Mat_E_Rm, typename Vec_E_R1>
      void fbe_1m_to_1m_1s(Rm&, R1&, const E&, const Mat_E_Rm&,
			   const Vec_E_R1&) const;

      template<class Rm, class E, typename Mat_E_Rm>
      void fbe_1m_to_1m(Rm&, const E&, const Mat_E_Rm&) const;

      template<class R1, class Em, typename Vec_Em_R1>
      void fbe_1m_to_1s(R1&, const Em&, const Vec_Em_R1&) const;

      void fbe_q_to_b_msk(Pb&, Psk&, const Pq&) const;
      void fbe_q_to_b(Pb&, const Pq&) const;
      void fbe_b_to_q_msk(Pq&, Psk&, const Pb&) const;
      void fbe_q_to_msk(Psk&, const Pq&) const;

      void sk(Pq&, Psk&, const Pb&) const;

      void SmMRq(Pb&, Psk&, const Pq&,
		 const std::array<main_signed_type, Pq::degree>&) const;

      std::array<main_double, nmQ> __M_q_i_inv;
      std::array<std::array<main_type, nmQ>, nmB> __M_q_i;
      std::array<main_type, nmB> __M_q;
      std::array<main_type, nmQ> __M_q_i_msk;
      main_type __M_q_msk;

      std::array<main_double, nmB> __M_b_i_inv;
      std::array<std::array<main_type, nmB>, nmQ> __M_b_i;
      std::array<main_type, nmQ> __M_b;

      void float_q_to_b (Pb&, const Pq&) const;
      void float_q_to_b_msk (Pb&, Psk&, const Pq&) const;
      void float_b_to_q (Pq&, const Pb&) const;

    };

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::allocate()
    {
      _qinv_modB = commons::mem_manag::alloc_aligned<Pb, 32>(1);
      _qinv_modmsk = commons::mem_manag::alloc_aligned<Psk, 32>(1);
      _xi_sk = commons::mem_manag::alloc_aligned<Psk, 32>(1);
      _xi_q = commons::mem_manag::alloc_aligned<Pq, 32>(2);
      _xi_b = commons::mem_manag::alloc_aligned<Pb, 32>(2);
      _xi_q_gt = commons::mem_manag::alloc_aligned<Pq, 32>(1);
      _xi_q_mt = commons::mem_manag::alloc_aligned<Pq, 32>(2);
    }
    
    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::copy_from
    (const base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>& other)
    {
      _Mq2b = other._Mq2b;
      _Mb2q = other._Mb2q;
      _Vq2sk = other._Vq2sk;
      _Vb2sk = other._Vb2sk;
      _Vsk2q = other._Vsk2q;
      _Mq2b_mt = other._Mq2b_mt;
      _Vq2sk_mt = other._Vq2sk_mt;
      _V2bsk_mt = other._V2bsk_mt;
      _xi_q[0] = other._xi_q[0];
      _xi_q[1] = other._xi_q[1];
      _xi_b[0] = other._xi_b[0];
      _xi_b[1] = other._xi_b[1];
      _xi_sk[0] = other._xi_sk[0];
      _xi_q_gt[0] = other._xi_q_gt[0];
      _xi_q_mt[0] = other._xi_q_mt[0];
      _xi_q_mt[1] = other._xi_q_mt[1];
      _qinv_modB[0] = other._qinv_modB[0];
      _qinv_modmsk[0] = other._qinv_modmsk[0];
      __M_q_i_inv = other.__M_q_i_inv;
      __M_q_i = other.__M_q_i;
      __M_q = other.__M_q;
      __M_q_i_msk = other.__M_q_i_msk;
      __M_q_msk = other.__M_q_msk;
      __M_b_i_inv = other.__M_b_i_inv;
      __M_b_i = other.__M_b_i;
      __M_b = other.__M_b;
    }
    
    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::base_conv_tool_set()
    {
      allocate();
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::base_conv_tool_set
    (const base_conv_tool_set &other)
    {
      allocate();
      copy_from(other);
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>&
    base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::operator=
    (const base_conv_tool_set &other)
    {
      copy_from(other);
      return *this;
    }    

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::~base_conv_tool_set()
    {
      commons::mem_manag::free_aligned<Psk>(1, _xi_sk);
      commons::mem_manag::free_aligned<Pq>(2, _xi_q);
      commons::mem_manag::free_aligned<Pb>(2, _xi_b);
      commons::mem_manag::free_aligned<Pq>(1, _xi_q_gt);
      commons::mem_manag::free_aligned<Pq>(2, _xi_q_mt);
      commons::mem_manag::free_aligned<Pb>(1, _qinv_modB);
      commons::mem_manag::free_aligned<Psk>(1, _qinv_modmsk);
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::generate()
    {
      constexpr size_t degree = Pq::degree;

      mpz_class q = 1;
      {
	mpz_class tmp;
	for (uint cmQ = 0; cmQ < nmQ; cmQ++)
	  {
	    mpz_set_ui(tmp.get_mpz_t(), Pq::base.params.P[Pq::indmod0+cmQ]);
	    q = q*tmp;
	  }
      }

      mpz_class b = 1;
      {
	mpz_class tmp;
	for (uint cmB = 0; cmB < nmB; cmB++)
	  {
	    mpz_set_ui(tmp.get_mpz_t(), Pb::base.params.P[Pb::indmod0+cmB]);
	    b = b*tmp;
	  }
      }

      mpz_class tmp_mt = 1;
      tmp_mt <<= Pmt::_bs_mt;

      // _xi_q: (q1/q mod q1, ....), (q/q1 mod q1, ....)
      // _xi_q_gt: (q1/q*gamma*t mod q1, ....)
      // _xi_q_mt: (q1/q*mt mod q1, ....)
      {
	mpz_class tmp_gt = 1, tmp, tmp2, tmp_q;
	tmp_gt <<= (parameters::_bs_t+parameters::_bs_g);
	auto *it_mt0 = _xi_q_mt[0].begin(), *it_mt1 = _xi_q_mt[1].begin();
	auto *it0 = _xi_q[0].begin(), *it1 = _xi_q[1].begin(),
	  *it_gt = _xi_q_gt->begin();
	for (uint cmQ = 0; cmQ < nmQ;
	     cmQ++, it0 += degree, it1 += degree, it_gt += degree
               , it_mt0 += degree, it_mt1 += degree
	     )
	  {
	    mpz_set_ui(tmp_q.get_mpz_t(), Pq::base.params.P[Pq::indmod0+cmQ]);
	    tmp = (q/tmp_q)%tmp_q;
	    std::fill(it1, it1+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

	    mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(), tmp_q.get_mpz_t());
	    std::fill(it0, it0+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

	    tmp2 = (tmp_gt*tmp)%tmp_q;  // qi/q*gt mod qi
	    std::fill(it_gt, it_gt+degree,
		      static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t())));
	    tmp2 = (tmp_mt*tmp)%tmp_q;  // qi/q*mt mod qi
	    std::fill(it_mt0, it_mt0+degree,
		      static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t())));

	    // q/qi/mt mod qi
	    mpz_invert(tmp2.get_mpz_t(), tmp2.get_mpz_t(), tmp_q.get_mpz_t());
	    std::fill(it_mt1, it_mt1+degree,
		      static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t())));
	  }
      }

      // _xi_b: (b1/b mod b1, ....), (b/b1 mod b1, ....)
      // _qinv_modB: (1/q mod b1, ...)
      {
	*_xi_b = 0;
	mpz_class tmp_b, tmp;
	auto *it0 = _xi_b[0].begin(), *it1 = _xi_b[1].begin(),
	  *it = _qinv_modB->begin();
	for (uint cmB = 0; cmB < nmB; cmB++, it0 += degree, it1 += degree, it += degree)
	  {
	    mpz_set_ui(tmp_b.get_mpz_t(), Pb::base.params.P[Pb::indmod0+cmB]);
	    tmp = (b/tmp_b)%tmp_b;
	    std::fill(it1, it1+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

	    mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(), tmp_b.get_mpz_t());
	    std::fill(it0, it0+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

	    mpz_invert(tmp.get_mpz_t(), q.get_mpz_t(), tmp_b.get_mpz_t());
	    std::fill(it, it+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
	  }
      }

      // _xi_sk: -1/b mod msk
      // _qinv_modmsk: 1/q mod msk
      {
	mpz_class tmp, mod;
	mpz_set_ui(mod.get_mpz_t(), Psk::base.params.P[Psk::indmod0]);
	mpz_invert(tmp.get_mpz_t(), b.get_mpz_t(), mod.get_mpz_t());
	tmp = mod-tmp;
	std::fill(_xi_sk->begin(), _xi_sk->end(),
		  static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

	mpz_invert(tmp.get_mpz_t(), q.get_mpz_t(), mod.get_mpz_t());
	std::fill(_qinv_modmsk->begin(), _qinv_modmsk->end(),
		  static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
      }

      // _Mq2b = (q/q1,...,q/qk) mod (b1,...,bl)
      // _Vq2sk = (q/q1,...,q/qk) mod msk
      // _Vq2mt = (-1/q1,...,-1/qk) mod gamma*t
      // _V2bsk_mt = (q/mt,...,q/mt, q/mt) mod (b1,...,bl, msk)
      // _Mq2b_mt = (q/q1/mt,...,q/qk/mt) mod (b1,...,bl)
      // _Vq2sk_mt = (q/q1/mt,...,q/qk/mt) mod msk
      {
	mpz_class tmp_sk;
	mpz_set_ui(tmp_sk.get_mpz_t(), Psk::base.params.P[Psk::indmod0]);

	for (uint cmQ = 0; cmQ < nmQ; cmQ++)
	  {
	    mpz_class tmp_q, mod, tmp;
	    mpz_set_ui(mod.get_mpz_t(), Pq::base.params.P[Pq::indmod0+cmQ]);
	    tmp_q = q/mod;

	    for (uint cmB = 0; cmB < nmB; cmB++)
	      {
		mpz_class tmp_b;
		mpz_set_ui(tmp_b.get_mpz_t(), Pb::base.params.P[Pb::indmod0+cmB]);
		tmp = tmp_q%tmp_b;

		_Mq2b[cmB][cmQ] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
		mpz_invert(tmp.get_mpz_t(), tmp_mt.get_mpz_t(), tmp_b.get_mpz_t());
		tmp = (tmp_q*tmp) % tmp_b;
		_Mq2b_mt[cmB][cmQ] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
	      }

	    tmp = tmp_q%tmp_sk;
	    _Vq2sk[cmQ] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
	    mpz_invert(tmp.get_mpz_t(), tmp_mt.get_mpz_t(), tmp_sk.get_mpz_t());
	    tmp = (tmp_q*tmp) % tmp_sk;
	    _Vq2sk_mt[cmQ] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
	  }

	mpz_class tmp;
	for (uint cmB = 0; cmB < nmB; cmB++)
	  {
	    mpz_class tmp_b;
	    mpz_set_ui(tmp_b.get_mpz_t(), Pb::base.params.P[Pb::indmod0+cmB]);
	    mpz_invert(tmp.get_mpz_t(), tmp_mt.get_mpz_t(), tmp_b.get_mpz_t());
	    tmp = (q*tmp)%tmp_b;

	    _V2bsk_mt[cmB] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
	  }

	mpz_invert(tmp.get_mpz_t(), tmp_mt.get_mpz_t(), tmp_sk.get_mpz_t());
	tmp = (tmp*q) % tmp_sk;
	_V2bsk_mt[nmB] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
      }

      // _Mb2q = (b/b1,...,b/bl) mod (q1,...,qk)
      // _Vsk2q = (-b,...,-b) mod (q1,...,qk)
      // _Vb2sk = (1/b1,...,1/bl) mod msk
      {
	mpz_class mod_sk;
	mpz_set_ui(mod_sk.get_mpz_t(), Psk::base.params.P[Psk::indmod0]);

	for (uint cmB = 0; cmB < nmB; cmB++)
	  {
	    mpz_class tmp, mod_b;
	    mpz_set_ui(mod_b.get_mpz_t(), Pb::base.params.P[Pb::indmod0+cmB]);
	    mpz_class tmp_b = b/mod_b;
	    for (uint cmQ = 0; cmQ < nmQ; cmQ++)
	      {
		mpz_class mod_q;
		mpz_set_ui(mod_q.get_mpz_t(), Pq::base.params.P[Pq::indmod0+cmQ]);

		if(cmB == 0)
		  {
		    mpz_class tmp_v = mod_q - (b%mod_q);
		    _Vsk2q[cmQ] = static_cast<main_type>(mpz_get_ui(tmp_v.get_mpz_t()));
		  }

		tmp = tmp_b%mod_q;
		_Mb2q[cmQ][cmB] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
	      }
	    mpz_invert(tmp.get_mpz_t(), mod_b.get_mpz_t(), mod_sk.get_mpz_t());
	    _Vb2sk[cmB] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
	  }
      }

      mpfr_set_default_prec (10000);
      {
	//HPR based extensions
	for (size_t i = 0; i < nmQ; i++)
	  {
	    mpz_class mod (Pq::base.params.P[Pq::indmod0+i]);
	    mpfr_t mpfr_q_i, w;
	    mpfr_init (mpfr_q_i);
	    mpfr_set_z (mpfr_q_i, mod.get_mpz_t (), MPFR_RNDN);
	    mpfr_ui_div (mpfr_q_i, 1, mpfr_q_i, MPFR_RNDN);
	    __M_q_i_inv[i] = mpfr_get_ld (mpfr_q_i, MPFR_RNDN);

	    mpfr_clear (mpfr_q_i);
	  }

	for (size_t i = 0; i < nmQ; i++)
	  {
	    mpz_class qi (Pq::base.params.P[Pq::indmod0+i]);
	    mpz_class q_ov_qi (q/qi);
	    for (size_t j = 0; j < nmB; j++)
	      {
		mpz_class pj (Pb::base.params.P[Pb::indmod0+j]);
		mpz_class qi_pj (q_ov_qi % pj);
		__M_q_i[j][i] = mpz_get_ui (qi_pj.get_mpz_t ());
	      }
	    mpz_class msk (Psk::base.params.P[Psk::indmod0]);
	    mpz_class qi_msk (q_ov_qi % msk);
	    __M_q_i_msk[i] = mpz_get_ui (qi_msk.get_mpz_t ());
	  }

	for (size_t i = 0; i < nmB; i++)
	  {
	    mpz_class pi (Pb::base.params.P[Pb::indmod0+i]);
	    mpz_class q_pi (q % pi);
	    __M_q[i] = mpz_get_ui (q_pi.get_mpz_t ());
	  }
	mpz_class msk (Psk::base.params.P[Psk::indmod0]);
	mpz_class q_msk (q % msk);
	__M_q_msk = mpz_get_ui (q_msk.get_mpz_t ());

	for (size_t i = 0; i < nmB; i++)
	  {
	    mpz_class mod (Pb::base.params.P[Pb::indmod0+i]);
	    mpfr_t mpfr_b_i;
	    mpfr_init (mpfr_b_i);
	    mpfr_set_z (mpfr_b_i, mod.get_mpz_t (), MPFR_RNDN);
	    mpfr_ui_div (mpfr_b_i, 1, mpfr_b_i, MPFR_RNDN);
	    __M_b_i_inv[i] = mpfr_get_ld (mpfr_b_i, MPFR_RNDN);

	    mpfr_clear (mpfr_b_i);
	  }

	for (size_t i = 0; i < nmB; i++)
	  {
	    mpz_class bi (Pb::base.params.P[Pb::indmod0+i]);
	    mpz_class b_ov_bi (b/bi);
	    for (size_t j = 0; j < nmQ; j++)
	      {
		mpz_class qj (Pq::base.params.P[Pq::indmod0+j]);
		mpz_class pi_qj (b_ov_bi % qj);
		__M_b_i[j][i] = mpz_get_ui (pi_qj.get_mpz_t ());
	      }
	  }

	for (size_t i = 0; i < nmQ; i++)
	  {
	    mpz_class qi (Pq::base.params.P[Pq::indmod0+i]);
	    mpz_class p_qi (b % qi);
	    __M_b[i] = mpz_get_ui (p_qi.get_mpz_t ());
	  }

      }
    }

    /// print fcts (debug)

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::print_Mq2b()
    {
      std::cout << "printing matrix q -> b: (q/q1,...,q/qk) mod (b1,...,bl)" <<
	std::endl << std::endl;
      for(auto* itB = std::begin(_Mq2b); itB < std::end(_Mq2b); ++itB)
	{
	  for(auto* itQ = std::begin(*itB); itQ < std::end(*itB); ++itQ)
	    std::cout << *itQ << " ";
	  std::cout << std::endl;
	}
      std::cout << std::endl;
      std::cout << "end matrix q -> b" << std::endl << std::endl;
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::print_Mb2q()
    {
      std::cout << "printing matrix b -> q: (b/b1,...,b/bl) mod (q1,...,qk)" <<
	std::endl << std::endl;
      for(auto* itQ = std::begin(_Mb2q); itQ < std::end(_Mb2q); ++itQ)
	{
	  for(auto* itB = std::begin(*itQ); itB < std::end(*itQ); ++itB)
	    std::cout << *itB << " ";
	  std::cout << std::endl;
	}
      std::cout << std::endl;
      std::cout << "end matrix b -> q" << std::endl << std::endl;
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::print_Vq2sk()
    {
      std::cout << "printing matrix q -> msk: (q/q1,...,q/qk) mod msk";
      std::cout << std::endl << std::endl;
      for(auto* itQ = std::begin(_Vq2sk); itQ < std::end(_Vq2sk); ++itQ)
	std::cout << *itQ << " ";
      std::cout << std::endl << std::endl;
      std::cout << "end matrix q -> msk";
      std::cout << std::endl << std::endl;
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::print_Vb2sk()
    {
      std::cout << "printing matrix b -> msk: (b/b1,...,b/bl, -1/b) mod msk";
      std::cout << std::endl << std::endl;
      for(auto* itB = std::begin(_Vb2sk); itB < std::end(_Vb2sk); ++itB)
	std::cout << *itB << " ";
      std::cout << std::endl << std::endl;
      std::cout << "end matrix b -> msk";
      std::cout << std::endl << std::endl;
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::print_Vsk2q()
    {
      std::cout << "printing matrix msk -> q: (-b,...,-b) mod (q1,...,qk)";
      std::cout << std::endl << std::endl;
      for(auto* it = std::begin(_Vsk2q); it < std::end(_Vsk2q); ++it)
	std::cout << *it << " ";
      std::cout << std::endl << std::endl;
      std::cout << "end matrix msk -> q";
      std::cout << std::endl << std::endl;
    }


    /// base conv.

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    template<class Rm, class R1, class E, typename Mat_E_Rm, typename Vec_E_R1>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::fbe_1m_to_1m_1s(
									   Rm &rm,
									   R1 &r1,
									   const E &e,
									   const Mat_E_Rm& mat_E_Rm,
									   const Vec_E_R1& vec_E_R1) const
    {
      // accumulates the conversion to rm and r1; e.g. r1 <-- fbe(e) + r1
      constexpr size_t degree = Pq::degree;

      std::array<main_greater_type, degree> accm, acc1;
      acc1.fill(0);

      auto iter_c_E_Rm = mat_E_Rm.begin();
      auto iter_c_E_R1 = vec_E_R1.begin();
      auto iter_rm = rm.begin();
      auto iter_acc1 = std::begin(acc1);

      constexpr size_t nmoduli = E::nmoduli;
      constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					Rm::nbits);

      main_type mod_R1 = R1::base.params.P[R1::indmod0];
      main_greater_type w0_R1 = static_cast<main_greater_type>
	(R1::base.params.Pn[R1::indmod0]);

      for(uint cm = 0; cm < Rm::nmoduli; cm++)
	{
	  accm.fill(0);
	  auto iter_c_E_Rm_bis = iter_c_E_Rm->begin();
	  main_type mod_Rm = Rm::base.params.P[cm+Rm::indmod0];
	  main_greater_type w0_Rm = static_cast<main_greater_type>
	    (Rm::base.params.Pn[cm+Rm::indmod0]);

	  auto* iter_e = e.begin();
	  uint cpt = 0;
	  while(iter_e < e.end())
	    {
	      if(cm == 0)
		iter_acc1 = acc1.begin();
	      if constexpr (nmoduli >= gap_lzr)
		  {
		    ++cpt;
		  }

	      for(auto iter_accm = accm.begin(); iter_accm < accm.end();
		  ++iter_accm, ++iter_e)
		{
		  main_greater_type tmp = static_cast<main_greater_type>(*iter_e);
		  *iter_accm += (*iter_c_E_Rm_bis)*tmp;
		  if(cm == 0)
		    {
		      *iter_acc1 += (*iter_c_E_R1)*tmp;
		      ++iter_acc1;
		    }
		  if constexpr (nmoduli >= gap_lzr)
				 {
				   if(cpt == gap_lzr)
				     {
				       main_type tmp = commons::reductions::IMR<Rm>(*iter_accm, w0_Rm, mod_Rm);
				       *iter_accm = static_cast<main_greater_type>(tmp);
				       if(cm == 0)
					 {
					   --iter_acc1;
					   main_type tmp = commons::reductions::IMR<R1>(*iter_acc1, w0_R1, mod_R1);
					   *iter_acc1 = static_cast<main_greater_type>(tmp);
					   ++iter_acc1;
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

	  for(auto iter_accm = accm.begin(); iter_accm < accm.end();
	      ++iter_accm, ++iter_rm)
	    {
	      *iter_accm += static_cast<main_greater_type>(*iter_rm);
	      *iter_rm = commons::reductions::IMR_full<Rm>(*iter_accm, w0_Rm, mod_Rm);
	    }
	  if(cm == 0)
	    {
	      auto iter_r1 = r1.begin();
	      for(auto iter_acc1 = acc1.begin(); iter_acc1 < acc1.end();
		  ++iter_acc1, ++iter_r1)
		{
		  *iter_acc1 += static_cast<main_greater_type>(*iter_r1);
		  *iter_r1 = commons::reductions::IMR_full<R1>(*iter_acc1, w0_R1, mod_R1);
		}
	    }
	  ++iter_c_E_Rm;
	}
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    template<class Rm, class E, typename Mat_E_Rm>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::fbe_1m_to_1m(
									Rm &rm,
									const E &e,
									const Mat_E_Rm& mat_E_Rm) const
    {
      // accumulates the conversion to rm and r1; e.g. r1 <-- fbe(e) + r1
      constexpr size_t degree = Pq::degree;
      std::array<main_greater_type, degree> accm;
      constexpr size_t nmoduli = E::nmoduli;
      constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					Rm::nbits);


      auto iter_c_E_Rm = mat_E_Rm.begin();
      auto iter_rm = rm.begin();

      for(uint cm = 0; cm < Rm::nmoduli; cm++)
	{
	  accm.fill(0);
	  auto iter_c_E_Rm_bis = iter_c_E_Rm->begin();
	  main_type mod_Rm = Rm::base.params.P[cm+Rm::indmod0];
	  main_greater_type w0_Rm = static_cast<main_greater_type>
	    (Rm::base.params.Pn[cm+Rm::indmod0]);

	  auto* iter_e = e.begin();
	  uint cpt = 0;

	  while(iter_e < e.end())
	    {
	      if constexpr (nmoduli >= gap_lzr)
		  {
		    ++cpt;
		  }

	      for(auto iter_accm = accm.begin(); iter_accm < accm.end();
		  ++iter_accm, ++iter_e)
		{
		  main_greater_type tmp = static_cast<main_greater_type>(*iter_e);
		  *iter_accm += (*iter_c_E_Rm_bis)*tmp;
		  if constexpr (nmoduli >= gap_lzr)
				 {
				   if(cpt == gap_lzr)
				     {
				       main_type tmp = commons::reductions::IMR<Rm>(*iter_accm, w0_Rm, mod_Rm);
				       *iter_accm = static_cast<main_greater_type>(tmp);
				     }
				 }
		}
	      if constexpr (nmoduli >= gap_lzr)
			     {
			       if(cpt == gap_lzr)
				 cpt = 0;
			     }

	      ++iter_c_E_Rm_bis;
	    }

	  for(auto iter_accm = accm.begin(); iter_accm < accm.end();
	      ++iter_accm, ++iter_rm)
	    {
	      *iter_accm += static_cast<main_greater_type>(*iter_rm);
	      *iter_rm = commons::reductions::IMR_full<Rm>(*iter_accm, w0_Rm, mod_Rm);
	    }
	  ++iter_c_E_Rm;
	}
    }


    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    template<class R1, class Em, typename Vec_Em_R1>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::fbe_1m_to_1s(R1 &r,
									const Em &e, const Vec_Em_R1 &coeffs_Em_R1) const
    {
      constexpr size_t degree = Pq::degree;
      std::array<main_greater_type, degree> acc;
      acc.fill(0);

      main_type mod_R = R1::base.params.P[R1::indmod0];
      main_greater_type w0_R = static_cast<main_greater_type>
	(R1::base.params.Pn[R1::indmod0]);
      constexpr size_t nmoduli = 2*Em::nmoduli;
      constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					R1::nbits);

      uint cpt = 0;

      auto iter_c_Em_R1 = coeffs_Em_R1.begin();
      auto iter_e = e.begin();
      while(iter_e < e.end())
	{
	  if constexpr (nmoduli >= gap_lzr)
	      {
		++cpt;
	      }

	  main_greater_type tmpv = static_cast<main_greater_type>(*iter_c_Em_R1++);
	  for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc)
	    {
	      *iter_acc += tmpv*static_cast<main_greater_type>(*iter_e++);
	      if constexpr (nmoduli >= gap_lzr)
			     {
			       if(cpt == gap_lzr)
				 {
				   main_type tmp = commons::reductions::IMR<R1>(*iter_acc, w0_R, mod_R);
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
	  *iter_r = commons::reductions::IMR_full<R1>(*iter_acc, w0_R, mod_R);
	}
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::float_q_to_b (Pb &pb,
									 const Pq &pq) const
    {
      fast_1m_to_1m_float<Pb, Pq> (pb, pq, __M_q_i_inv, __M_q_i, __M_q);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::float_q_to_b_msk (Pb &pb,
									     Psk &psk, const Pq &pq) const
    {
      fast_1m_to_1m_1s_float<Pb, Psk, Pq> (pb, psk, pq, __M_q_i_inv, __M_q_i,
					   __M_q_i_msk, __M_q, __M_q_msk);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::float_b_to_q (Pq &pq,
									 const Pb &pb) const
    {
      fast_1m_to_1m_float<Pq, Pb> (pq, pb, __M_b_i_inv, __M_b_i, __M_b);
    }


    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::fbe_q_to_b_msk(
									  Pb &pb,
									  Psk &psk,
									  const Pq &pq) const
    {
      this->fbe_1m_to_1m_1s<Pb, Psk, Pq, matBEq2b, vecBEq2sk>(
							      pb,
							      psk,
							      pq,
							      _Mq2b,
							      _Vq2sk);
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::fbe_q_to_b(
								      Pb &pb,
								      const Pq &pq) const
    {
      this->fbe_1m_to_1m<Pb, Pq, matBEq2b>(pb,
					   pq,
					   _Mq2b);
    }


    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::fbe_b_to_q_msk(
									  Pq &pq,
									  Psk &psk,
									  const Pb &pb) const
    {
      this->fbe_1m_to_1m_1s<Pq, Psk, Pb, matBEb2q, vecBEb2sk>(
							      pq,
							      psk,
							      pb,
							      _Mb2q,
							      _Vb2sk);
    }


    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::fbe_q_to_msk(
									Psk &psk,
									const Pq &pq) const
    {
      this->fbe_1m_to_1s<Psk, Pq, vecBEq2sk>(
					     psk,
					     pq,
					     _Vq2sk);
    }


    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::sk(
							      Pq &pq,
							      Psk &psk,
							      const Pb &pb) const
    {
      constexpr size_t degree = Pq::degree;
      // 1st step
      this->fbe_b_to_q_msk(pq, psk, pb);

      // 2nd step
      std::array<main_signed_type, degree> acc_centeredsk;

      main_type m_sk = Psk::base.params.P[Psk::indmod0];
      main_type m_sk_ov2 = (m_sk>>1);

      auto* iter_esk = psk.begin();
      for(auto* it_sk = acc_centeredsk.begin(); it_sk < acc_centeredsk.end();
	  ++it_sk, ++iter_esk)
	*it_sk = (*iter_esk <= m_sk_ov2) ? static_cast<main_signed_type>
	  (*iter_esk) : static_cast<main_signed_type>(*iter_esk) -
	  static_cast<main_signed_type>(m_sk);

      auto* iter_vsk2q = _Vsk2q.begin();
      auto* iter_q = pq.begin();
      for(uint cm = 0; cm < nmQ; cm++, ++iter_vsk2q)
	{
	  main_type mod_q = Pq::base.params.P[cm+Pq::indmod0];
	  main_greater_type w0_q = static_cast<main_greater_type>(Pq::base.params.Pn[cm
										     +Pq::indmod0]);

	  for(auto* it_sk = acc_centeredsk.begin(); it_sk < acc_centeredsk.end(); ++it_sk)
	    {
	      main_greater_type tmp = (*it_sk >= 0) ? static_cast<main_greater_type>
		(*it_sk) : static_cast<main_greater_type>(mod_q) -
		static_cast<main_greater_type>(-(*it_sk));

	      tmp = tmp * (*iter_vsk2q) + static_cast<main_greater_type>(*iter_q);
	      *iter_q++ = commons::reductions::IMR_full<Pq>(tmp, w0_q, mod_q);
	    }
	}
    }

    template <class Pq, class Pb, class Psk, class parameters, class Pmt>
    void base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt>::SmMRq(Pb &pb, Psk &psk,
								 const Pq &pq, const std::array<main_signed_type, Pq::degree> &acc) const
    {

      constexpr size_t nmoduli = Pq::nmoduli;
      constexpr size_t gap_lzr = 1ULL<<(8 * sizeof(main_greater_type) - 2 *
					Pb::nbits);
      constexpr size_t degree = Pq::degree;

      std::array<main_greater_type, degree> accb, accsk;
      accsk.fill(0);

      auto iter_mat_b = _Mq2b_mt.begin();
      auto iter_vec_sk = _Vq2sk_mt.begin();
      auto iter_polb = pb.begin();
      auto iter_polsk = psk.begin();
      auto iter_accsk = std::begin(accsk);

      main_type mod_sk = Psk::base.params.P[Psk::indmod0];
      main_greater_type w0_sk = static_cast<main_greater_type>
	(Psk::base.params.Pn[Psk::indmod0]);

      auto iter_v = _V2bsk_mt.begin();
      auto iter_v_last = std::prev(_V2bsk_mt.end());

      for(uint cm = 0; cm < nmB; cm++)
	{
	  accb.fill(0);
	  auto iter_mat_b_bis = iter_mat_b->begin();
	  main_type mod_b = Pb::base.params.P[cm+Pb::indmod0];
	  main_greater_type w0_b = static_cast<main_greater_type>(Pb::base.params.Pn[cm
										     +Pb::indmod0]);

	  auto* iter_polq = pq.begin();
	  uint cpt = 0;
	  while(iter_polq < pq.end())
	    {
	      if(cm == 0)
		iter_accsk = accsk.begin();
	      if constexpr (nmoduli >= gap_lzr)
		  {
		    ++cpt;
		  }

	      for(auto iter_accb = accb.begin(); iter_accb < accb.end();
		  ++iter_accb, ++iter_polq)
		{
		  main_greater_type tmp = static_cast<main_greater_type>(*iter_polq);
		  *iter_accb += (*iter_mat_b_bis)*tmp;
		  if(cm == 0)
		    {
		      *iter_accsk += (*iter_vec_sk)*tmp;
		      ++iter_accsk;
		    }

		  if constexpr (nmoduli >= gap_lzr)
				 {
		      if(cpt == gap_lzr)
			{
			  main_type tmp = commons::reductions::IMR<Pb>(*iter_accb, w0_b, mod_b);
			  *iter_accb = static_cast<main_greater_type>(tmp);
			  if(cm == 0)
			    {
			      --iter_accsk;
			      main_type tmp = commons::reductions::IMR<Psk>(*iter_accsk, w0_sk, mod_sk);
			      *iter_accsk = static_cast<main_greater_type>(tmp);
			      ++iter_accsk;
			    }
			}
		    }
		}
	      if constexpr (nmoduli >= gap_lzr)
			     {
			       if(cpt == gap_lzr)
				 cpt = 0;
			     }

	      ++iter_mat_b_bis;
	      if(cm == 0)
		++iter_vec_sk;
	    }

	  ++iter_mat_b;

	  // getting content of acc

	  auto iter_acc = acc.begin();
	  for(auto iter_accb = accb.begin(); iter_accb < accb.end();
	      ++iter_accb, ++iter_polb, ++iter_acc)
	    {
	      if(*iter_acc >= 0)
		*iter_accb += (main_greater_type)(*iter_v)*(*iter_acc);
	      else
		*iter_accb -= (main_greater_type)(*iter_v)*(-*iter_acc);

	      *iter_accb += static_cast<main_greater_type>(*iter_polb);
	      *iter_polb = commons::reductions::IMR_full<Pb>(*iter_accb, w0_b, mod_b);
	    }

	  ++iter_v;

	  if(cm == 0)
	    {
	      auto iter_acc = acc.begin();
	      for(auto iter_accsk = accsk.begin(); iter_accsk < accsk.end();
		  ++iter_accsk, ++iter_polsk, ++iter_acc)
		{
		  if(*iter_acc >= 0)
		    *iter_accsk += (main_greater_type)(*iter_v_last)*(*iter_acc);
		  else
		    *iter_accsk -= (main_greater_type)(*iter_v_last)*(-*iter_acc);

		  *iter_accsk += static_cast<main_greater_type>(*iter_polsk);
		  *iter_polsk = commons::reductions::IMR_full<Psk>(*iter_accsk, w0_sk, mod_sk);
		}
	    }
	}
      /*

	this->fbe_1m_to_1m_1s<Pb, Psk, Pq, matBEq2b, vecBEq2sk>(
	pb,
	psk,
	pq,
	_Mq2b_mt,
	_Vq2sk_mt);

	// reduction is simply (p + q*acc)/mt mod (b1, ...)
	auto iter_b = pb.begin();
	auto iter_v = _V2bsk_mt.begin();
	for(uint cm = 0; cm < nmB; cm++)
	{
	main_type mod = nfl::params<main_type>::P[cm+Pb::indmod0];
	//
	main_greater_type w0 = (main_greater_type)nfl::params<main_type>::Pn[cm+Pb::indmod0];
	//
	for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_b)
	{
	main_greater_type tmp = (*iter_acc >= 0) ? (main_greater_type)(*iter_v)*(*iter_acc) : (main_greater_type)(*iter_v)*(-*iter_acc);
	main_type tmpv = reductions::IMR_full(tmp, w0, mod);
	*iter_b += (*iter_acc >= 0) ? tmpv : mod-tmpv;
	if(*iter_b >= mod)
	*iter_b -= mod;
	}
	++iter_v;
	}

	auto iter_sk = psk.begin();
	main_type mod_msk = nfl::params<main_type>::P[Psk::indmod0];
	//
	main_greater_type w0 = (main_greater_type)nfl::params<main_type>::Pn[Psk::indmod0];
	//
	for(auto iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_sk)
	{
	main_greater_type tmp = (*iter_acc >= 0) ? (main_greater_type)(*iter_v)*(*iter_acc) : (main_greater_type)(*iter_v)*(-*iter_acc);
	main_type tmpv = reductions::IMR_full(tmp, w0, mod_msk);

	*iter_sk += (*iter_acc >= 0) ? tmpv : mod_msk-tmpv;

	if(*iter_sk >= mod_msk)
	*iter_sk -= mod_msk;
	}
      */
    }

    ////////////////

    template <class R, class E, class type_m>
    void copy_paste_conversion(R &r, const E &e, type_m mod_e)
    {
      using main_type = typename R::value_type;
      using main_greater_type = typename R::greater_value_type;
      using main_signed_type = typename R::signed_value_type;
      constexpr size_t degree = R::degree;

      std::array<main_signed_type, degree> acc;
      acc.fill(0);

      auto* iter_e = e.begin();
      for(auto* iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_e)
	*iter_acc = (*iter_e <= (mod_e>>1)) ? static_cast<main_signed_type>
	  (*iter_e) : -static_cast<main_signed_type>(mod_e-*iter_e);

      auto* iter_r = r.begin();
      for(uint cm = 0; cm < R::nmoduli; cm++)
	{
	  main_type mod = r.get_modulus(cm);
	  for(auto* iter_acc = acc.begin(); iter_acc < acc.end(); ++iter_acc, ++iter_r)
	    *iter_r += (*iter_acc < 0) ? mod-static_cast<main_type>(-*iter_acc) :
	      static_cast<main_type>(*iter_acc);
	}
    }

  }
}
#endif /* fvrns_be_hpr_h */
