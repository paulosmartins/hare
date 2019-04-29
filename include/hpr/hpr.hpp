#ifndef hpr_h
#define hpr_h

#include "hpr/common_hpr.hpp"
#include "hpr/be_hpr.hpp"

namespace hpr
{
  namespace digit
  {
    // structure containing a digit with two components (c0,c1)
    // both in bases q, b and m_sk
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    struct poly_digit
    {
      typedef Pq _Pq;
      typedef Pb _Pb;
      typedef Psk _Psk;
      using main_type = typename Pq::value_type;
      using greater_main_type = typename Pq::greater_value_type;
      using main_signed_type = typename Pq::signed_value_type;

      static constexpr uint __nmQ = Pq::nmoduli;
      static constexpr uint __nmB = Pb::nmoduli;
      static constexpr uint __imB = Pb::indmod0;
      static constexpr uint __imSK = Psk::indmod0;

      Pq* __cq;
      Pb* __cb;
      Psk* __csk;

      typedef typename
      base_conversions::base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt> be_tool_set;
      const be_tool_set* __be_tool_set; // stored in context structure

      void allocate();
      void copy_from(const poly_digit &other);
      poly_digit();
      poly_digit(const poly_digit &other);
      poly_digit& operator=(const poly_digit &other);
      ~poly_digit();

      void set_mat_set_conv(const be_tool_set*);

      void setq0  ();
      void setb0  ();
      void setsk0 ();
      void set0   ();
      void set0   (int);
      void set0_float   ();
      void set0_float   (int);

      void setq   (const poly_digit*);
      void setb   (const poly_digit*);
      void setsk  (const poly_digit*);
      void set    (const poly_digit*);
      void set    (int, const poly_digit*, int);
      void set_float    (const poly_digit*);
      void set_float    (int, const poly_digit*, int);


      bool operator==(const poly_digit&) const;
      bool operator!=(const poly_digit&) const;

      void add        (int, const poly_digit*, int, const poly_digit*, int);
      void add_float        (int, const poly_digit*, int, const poly_digit*, int);
      void add        (int, const poly_digit*, int);
      void add_float        (int, const poly_digit*, int);

      void sub        (int, const poly_digit*, int, const poly_digit*, int);
      void sub        (int, const poly_digit*, int);
      void sub_float        (int, const poly_digit*, int, const poly_digit*, int);
      void sub_float        (int, const poly_digit*, int);

      void mul        (int, const poly_digit*, int, const poly_digit*, int);
      void mul        (int, const poly_digit*, int);
      void mul_float        (int, const poly_digit*, int, const poly_digit*, int);
      void mul_float        (int, const poly_digit*, int);

      void mul_acc    (int, const poly_digit*, int, const poly_digit*, int);
      void mul_acc_float    (int, const poly_digit*, int, const poly_digit*, int);

      void copy(int);

      void inv_ntt(int);
      void inv_ntt();
      void inv_ntt_float(int);
      void inv_ntt_float();
      void ntt(int);
      void ntt();
      void ntt_float(int);
      void ntt_float();
      void nttq();

      void gen_noise(int);

      void fbe(int);
      void fbe();
      void fbe_float(int);
      void fbe_float();
      void skbe(int, poly_digit*);
      void skbe(poly_digit*);
      void float_back(int, poly_digit*);
      void float_back(poly_digit*);
      void float_forward (int);
      void float_forward ();
      void float_forward_sk (int);
      void float_forward_sk ();
      void fbe_q_to_msk(int, poly_digit*, int);

      void exact_conv_gmp_q_to_bsk(int, poly_digit*, int);
      void exact_conv_gmp_q_to_bsk(int);

      void spread_smallpoly_q2bsk(int);

      void print_q(int);
      void print_b(int);
      void print_sk(int);
    };

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void
    poly_digit<Pq, Pb, Psk, parameters, Pmt>::allocate()
    {
      this->__cq = commons::mem_manag::alloc_aligned<Pq, 32>(2);
      this->__cb = commons::mem_manag::alloc_aligned<Pb, 32>(2);
      this->__csk = commons::mem_manag::alloc_aligned<Psk, 32>(2);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void
    poly_digit<Pq, Pb, Psk, parameters, Pmt>::copy_from
    (const poly_digit<Pq, Pb, Psk, parameters, Pmt> &other)
    {
      __cq[0] = other.__cq[0];
      __cq[1] = other.__cq[1];
      __cb[0] = other.__cb[0];
      __cb[1] = other.__cb[1];
      __csk[0] = other.__csk[0];
      __csk[1] = other.__csk[1];
      //shallow copy
      __be_tool_set = other.__be_tool_set;
    }

    //////////////// constructor ////////////////
    /////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    poly_digit<Pq, Pb, Psk, parameters, Pmt>::poly_digit()
    {
      allocate();
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    poly_digit<Pq, Pb, Psk, parameters, Pmt>::poly_digit
    (const poly_digit<Pq, Pb, Psk, parameters, Pmt>& other)
    {
      allocate();
      copy_from(other);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    poly_digit<Pq, Pb, Psk, parameters, Pmt>&
    poly_digit<Pq, Pb, Psk, parameters, Pmt>::operator=
    (const poly_digit<Pq, Pb, Psk, parameters, Pmt>& other)
    {
      copy_from(other);
      return *this;
    }

    //////////////// destructor ////////////////
    ////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    poly_digit<Pq, Pb, Psk, parameters, Pmt>::~poly_digit()
    {
      commons::mem_manag::free_aligned<Pq>(2, this->__cq);
      commons::mem_manag::free_aligned<Pb>(2, this->__cb);
      commons::mem_manag::free_aligned<Psk>(2, this->__csk);
    }

    //////////////// setters ////////////////
    /////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set_mat_set_conv(
								    const be_tool_set* ms)
    {
      this->__be_tool_set = ms;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::setq0()
    {
      this->__cq[0] = 0;
      this->__cq[1] = 0;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::setb0()
    {
      this->__cb[0] = 0;
      this->__cb[1] = 0;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::setsk0()
    {
      this->__csk[0] = 0;
      this->__csk[1] = 0;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set0()
    {
      this->setq0();
      this->setb0();
      this->setsk0();
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set0(int j)
    {
      this->__cq[j] = 0;
      this->__cb[j] = 0;
      this->__csk[j] = 0;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set0_float()
    {
      this->setq0();
      this->setb0();
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set0_float(int j)
    {
      this->__cq[j] = 0;
      this->__cb[j] = 0;
    }


    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::setq(const poly_digit *o)
    {
      this->__cq[0] = o->__cq[0];
      this->__cq[1] = o->__cq[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::setb(const poly_digit *o)
    {
      this->__cb[0] = o->__cb[0];
      this->__cb[1] = o->__cb[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::setsk(const poly_digit *o)
    {
      this->__csk[0] = o->__csk[0];
      this->__csk[1] = o->__csk[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set(const poly_digit *o)
    {
      this->setq(o);
      this->setb(o);
      this->setsk(o);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set(int i,
							      const poly_digit *o, int j)
    {
      this->__cq[i] = o->__cq[j];
      this->__cb[i] = o->__cb[j];
      this->__csk[i] = o->__csk[j];
    }
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set_float(const poly_digit *o)
    {
      this->setq(o);
      this->setb(o);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::set_float(int i,
								    const poly_digit *o, int j)
    {
      this->__cq[i] = o->__cq[j];
      this->__cb[i] = o->__cb[j];
    }

    //////////////// operations ////////////////
    ////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    bool poly_digit<Pq, Pb, Psk, parameters, Pmt>::operator==
    (const poly_digit &o) const
    {
      for(int i = 0; i < 2; i++)
	if ( (std::equal(this->__cq[i].begin(), this->__cq[i].begin(),
			 o.__cq[i].begin()) == false)
	     || (std::equal(this->__cb[i].begin(), this->__cb[i].begin(),
			    o.__cb[i].begin()) == false)
	     || (std::equal(this->__csk[i].begin(), this->__csk[i].begin(),
			    o.__csk[i].begin()) == false) )
	  return false;
      return true;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    bool poly_digit<Pq, Pb, Psk, parameters, Pmt>::operator!=
    (const poly_digit &o) const
    {
      for(int i = 0; i < 2; i++)
	if ( (std::equal(this->__cq[i].begin(), this->__cq[i].begin(),
			 o.__cq[i].begin()) == false)
	     || (std::equal(this->__cb[i].begin(), this->__cb[i].begin(),
			    o.__cb[i].begin()) == false)
	     || (std::equal(this->__csk[i].begin(), this->__csk[i].begin(),
			    o.__csk[i].begin()) == false) )
	  return true;
      return false;
    }

    /********** addition **********/

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::add(int i0,
							      const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 = other1.c_i1 + other2.c_i2

      this->__cq[i0] = other1->__cq[i1] + other2->__cq[i2];
      this->__cb[i0] = other1->__cb[i1] + other2->__cb[i2];
      this->__csk[i0] = other1->__csk[i1] + other2->__csk[i2];
    }
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::add_float(int i0,
								    const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 = other1.c_i1 + other2.c_i2

      this->__cq[i0] = other1->__cq[i1] + other2->__cq[i2];
      this->__cb[i0] = other1->__cb[i1] + other2->__cb[i2];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::add(int i0,
							      const poly_digit *other1, int i1)
    {
      // this.c_i0 += other1.c_i1

      this->add(i0, this, i0, other1, i1);
    }
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::add_float(int i0,
								    const poly_digit *other1, int i1)
    {
      // this.c_i0 += other1.c_i1

      this->add_float(i0, this, i0, other1, i1);
    }

    /********** subtraction **********/

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::sub(int i0,
						       const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 = other1.c_i1 - other2.c_i2

      this->__cq[i0] = other1->__cq[i1] - other2->__cq[i2];
      this->__cb[i0] = other1->__cb[i1] - other2->__cb[i2];
      this->__csk[i0] = other1->__csk[i1] - other2->__csk[i2];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::sub(int i0,
						       const poly_digit *other1, int i1)
    {
      // this.c_i0 -= other1.c_i1

      this->sub(i0, this, i0, other1, i1);
    }
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::sub_float(int i0,
							     const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 = other1.c_i1 - other2.c_i2

      this->__cq[i0] = other1->__cq[i1] - other2->__cq[i2];
      this->__cb[i0] = other1->__cb[i1] - other2->__cb[i2];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::sub_float(int i0,
							     const poly_digit *other1, int i1)
    {
      // this.c_i0 -= other1.c_i1

      this->sub_float(i0, this, i0, other1, i1);
    }

    /********** multiplication **********/

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::mul(int i0,
						       const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 = other1.c_i1 * other2.c_i2

      this->__cq[i0] = other1->__cq[i1] * other2->__cq[i2];
      this->__cb[i0] = other1->__cb[i1] * other2->__cb[i2];
      this->__csk[i0] = other1->__csk[i1] * other2->__csk[i2];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::mul(int i0,
						       const poly_digit *other1, int i1)
    {
      // this.c_i0 *= other1.c_i1

      this->mul(i0, this, i0, other1, i1);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::mul_float(int i0,
							     const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 = other1.c_i1 * other2.c_i2

      this->__cq[i0] = other1->__cq[i1] * other2->__cq[i2];
      this->__cb[i0] = other1->__cb[i1] * other2->__cb[i2];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::mul_float(int i0,
							     const poly_digit *other1, int i1)
    {
      // this.c_i0 *= other1.c_i1

      this->mul(i0, this, i0, other1, i1);
    }


    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::mul_acc(int i0,
								  const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 += other1.c_i1 * other2.c_i2
      this->__cq[i0] = this->__cq[i0] + other1->__cq[i1] * other2->__cq[i2];
      this->__cb[i0] = this->__cb[i0] + other1->__cb[i1] * other2->__cb[i2];
      this->__csk[i0] = this->__csk[i0] + other1->__csk[i1] * other2->__csk[i2];
    }
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::mul_acc_float(int i0,
									const poly_digit *other1, int i1, const poly_digit *other2, int i2)
    {
      // this.c_i0 += other1.c_i1 * other2.c_i2
      this->__cq[i0] = this->__cq[i0] + other1->__cq[i1] * other2->__cq[i2];
      this->__cb[i0] = this->__cb[i0] + other1->__cb[i1] * other2->__cb[i2];
    }

    /********** misc. **********/

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::copy(int i)
    {
      // this.c_0 = this.c_1 or this.c_1 = this.c_1
      this->__cq[1-i] = this->__cq[i];
      this->__cb[1-i] = this->__cb[i];
      this->__csk[1-i] = this->__csk[i];
    }

    /********** ntt related **********/

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::inv_ntt(int i)
    {
      // inv_ntt on this.c_i

      this->__cq[i].invntt_pow_invphi();
      this->__cb[i].invntt_pow_invphi();
      this->__csk[i].invntt_pow_invphi();
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::inv_ntt()
    {
      // inv_ntt on this.c_0 and this.c_1

      this->inv_ntt(0);
      this->inv_ntt(1);
    }
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::inv_ntt_float(int i)
    {
      // inv_ntt on this.c_i

      this->__cq[i].invntt_pow_invphi();
      this->__cb[i].invntt_pow_invphi();
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::inv_ntt_float()
    {
      // inv_ntt on this.c_0 and this.c_1

      this->inv_ntt_float(0);
      this->inv_ntt_float(1);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::ntt(int i)
    {
      // ntt on this.c_i

      this->__cq[i].ntt_pow_phi();
      this->__cb[i].ntt_pow_phi();
      this->__csk[i].ntt_pow_phi();
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::ntt_float(int i)
    {
      // ntt on this.c_i

      this->__cq[i].ntt_pow_phi();
      this->__cb[i].ntt_pow_phi();
    }


    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::ntt()
    {
      // ntt on this.c_0 and this.c_1

      this->ntt(0);
      this->ntt(1);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::ntt_float()
    {
      // ntt on this.c_0 and this.c_1

      this->ntt_float(0);
      this->ntt_float(1);
    }


    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    inline void poly_digit<Pq, Pb, Psk, parameters, Pmt>::nttq()
    {
      // ntt on this.c_0 and this.c_1 in base q only

      this->__cq[0].ntt_pow_phi();
      this->__cq[1].ntt_pow_phi();
    }

    //////////////// noise generation ////////////////
    //////////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::gen_noise(int i)
    {
      using main_type = typename Pq::value_type;
      using main_greater_type = typename Pq::greater_value_type;
      using main_signed_type = typename Pq::signed_value_type;

      // i = 0: (c_0, c_1) <- (e, 0)
      // i = 1: (c_0, c_1) <- (0, e)
      // i = 2: (c_0, c_1) <- (e, e)
      // else: (c_0, c_1) <- (e0, e1)

      if(i == 0 || i == 1)
	{
	  this->__cq[i] = nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng);
	  spread_smallpoly_q2bsk(i);
	}
      else
	{
	  this->__cq[0] = nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng);
	  if(i == 2)
	    this->__cq[1] = this->__cq[0];
	  else
	    this->__cq[1] = nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng);

	  spread_smallpoly_q2bsk(0);
	  spread_smallpoly_q2bsk(1);
	}
    }

    //////////////// base conversions ////////////////
    //////////////////////////////////////////////////

    // fast base conversion : extending c mod q toward base b, m_sk
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::fbe(int i)
    {
      // extending this.c_i mod q to Z/bZ x Z/mskZ

      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[0];
      this->__cb[i] = 0;
      this->__csk[i] = 0;
      this->__be_tool_set->fbe_q_to_b_msk(this->__cb[i], this->__csk[i],
					  this->__cq[i]);
      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::fbe_float(int i)
    {
      // extending this.c_i mod q to Z/bZ x Z/mskZ

      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[0];
      this->__cb[i] = 0;
      this->__be_tool_set->fbe_q_to_b(this->__cb[i], this->__cq[i]);
      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::fbe()
    {
      // extending this.c_0 mod q and this.c_1 to Z/bZ x Z/mskZ

      this->fbe(0);
      this->fbe(1);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::fbe_float()
    {
      // extending this.c_0 mod q and this.c_1 to Z/bZ x Z/mskZ

      this->fbe_float(0);
      this->fbe_float(1);
    }

    // fast base conversion : extending c mod q toward base b, m_sk
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::float_forward(int i)
    {
      // extending this.c_i mod q to Z/bZ x Z/mskZ

      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[0];
      this->__cb[i] = 0;
      this->__be_tool_set->float_q_to_b (this->__cb[i], this->__cq[i]);
      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::float_forward()
    {
      // extending this.c_0 mod q and this.c_1 to Z/bZ x Z/mskZ

      this->float_forward(0);
      this->float_forward(1);
    }

    // fast base conversion : extending c mod q toward base b, m_sk
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::float_forward_sk(int i)
    {
      // extending this.c_i mod q to Z/bZ x Z/mskZ

      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[0];
      this->__cb[i] = 0;
      this->__be_tool_set->float_q_to_b_msk (this->__cb[i], this->__csk[i],
					     this->__cq[i]);
      this->__cq[i] = this->__cq[i] * this->__be_tool_set->_xi_q[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::float_forward_sk()
    {
      // extending this.c_0 mod q and this.c_1 to Z/bZ x Z/mskZ

      this->float_forward_sk(0);
      this->float_forward_sk(1);
    }


    // exact base conversion : extending this mod b*msk towards base q for c
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::skbe(int i, poly_digit *other1)
    {
      // Shenoy-Kumaresan conversion of this.c_i mod b*msk
      // result stored in other1.c_i mod q

      this->__cb[i] = this->__cb[i] * __be_tool_set->_xi_b[0];
      Psk esk = this->__csk[i] * __be_tool_set->_xi_sk[0];
      this->__be_tool_set->sk(other1->__cq[i], esk, this->__cb[i]);
      this->__cb[i] = this->__cb[i] * __be_tool_set->_xi_b[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::skbe(poly_digit *other1)
    {
      // Shenoy-Kumaresan conversion of this.c_0 and this.c_1 mod b*msk
      // result stored in other1.c_0 and other1.c_1 mod q

      this->skbe(0, other1);
      this->skbe(1, other1);
    }



    // exact base conversion : extending this mod b*msk towards base q for c
    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::float_back(int i,
							      poly_digit *other1)
    {
      // Shenoy-Kumaresan conversion of this.c_i mod b*msk
      // result stored in other1.c_i mod q

      this->__cb[i] = this->__cb[i] * __be_tool_set->_xi_b[0];
      Pq tmp;
      this->__be_tool_set->float_b_to_q(tmp, this->__cb[i]);
      other1->__cq[i] = other1->__cq[i] + tmp;
      this->__cb[i] = this->__cb[i] * __be_tool_set->_xi_b[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::float_back(poly_digit *other1)
    {
      // Shenoy-Kumaresan conversion of this.c_0 and this.c_1 mod b*msk
      // result stored in other1.c_0 and other1.c_1 mod q

      this->float_back(0, other1);
      this->float_back(1, other1);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::fbe_q_to_msk(int i0,
								poly_digit *other1, int i1)
    {
      // fast conversion from Z/qZ to Z/mskZ of this.c_i0
      // result stored in other1.c_i1

      this->__cq[i0] = this->__cq[i0] * __be_tool_set->_xi_q[0];
      this->__be_tool_set->fbe_q_to_msk(other1->__csk[i1], this->__cq[i0]);
      this->__cq[i0] = this->__cq[i0] * __be_tool_set->_xi_q[1];
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::exact_conv_gmp_q_to_bsk(int i0,
									   poly_digit* other1, int i1)
    {
      // exact base conversion from q to bsk by using gmp for reduction modulo q
      // conversion of this.c_i0, result stored in other1.c_i1
      constexpr size_t degree = Pq::degree;
      mpz_class q = 1, tmp;
      std::vector<mpz_class> moduli;
      for (uint cmQ = 0; cmQ < this->__nmQ; cmQ++)
	{
	  mpz_set_ui(tmp.get_mpz_t(), Pq::base.params.P[Pq::indmod0+cmQ]);
	  moduli.push_back(tmp);
	  q = q*tmp;
	}

      this->__cq[i0] = this->__cq[i0] * __be_tool_set->_xi_q[0];

      for(size_t deg = 0; deg < degree; ++deg)
	{
	  mpz_class coeff = 0;
	  for (uint cmQ = 0; cmQ < this->__nmQ; cmQ++)
	    {
	      mpz_set_ui(tmp.get_mpz_t(), this->__cq[i0](cmQ, deg));
	      coeff = coeff + tmp * q / moduli[cmQ];
	    }

	  coeff = coeff%q;

	  for (uint cmB = 0; cmB < this->__nmB; cmB++)
	    {
	      mpz_set_ui(tmp.get_mpz_t(), Pb::base.params.P[Pb::indmod0+cmB]);
	      tmp = coeff%tmp;
	      other1->__cb[i1].data()[cmB*degree+deg] = static_cast<main_type>(mpz_get_ui(
											  tmp.get_mpz_t()));
	    }

	  mpz_set_ui(tmp.get_mpz_t(), Psk::base.params.P[Psk::indmod0]);
	  tmp = coeff%tmp;
	  other1->__csk[i1].data()[deg] = static_cast<main_type>(mpz_get_ui(
									    tmp.get_mpz_t()));
	}

      this->__cq[i0] = this->__cq[i0] * __be_tool_set->_xi_q[1];
      moduli.clear();
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::exact_conv_gmp_q_to_bsk(int i)
    {
      // exact base conversion from q to bsk by using gmp for reduction modulo q
      // conversion of this.c_i0, result stored in this.c_i0

      exact_conv_gmp_q_to_bsk(i, this, i);
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::spread_smallpoly_q2bsk(int i)
    {
      // conversion of small polynomial modulo q to base b*msk
      // just checking the sign, and copy paste
      constexpr size_t degree = Pq::degree;

      main_type modq = this->__cq[i].get_modulus(0);
      auto* it_pb = this->__cb[i].begin();
      auto* it_psk = this->__csk[i].begin();
      main_type modsk = this->__csk[i].get_modulus(0);
      for(uint cmB = 0; cmB < this->__nmB; cmB++)
	{
	  main_type modB = this->__cb[i].get_modulus(cmB);
	  auto* it_pq = this->__cq[i].begin();
	  for(size_t deg = 0; deg != degree; ++deg, ++it_pq, ++it_pb)
	    {
	      *it_pb = ((*it_pq<<1) >= modq) ? modB - (modq - *it_pq) : *it_pq;
	      if(cmB == 0)
		*it_psk++ = ((*it_pq<<1) >= modq) ? modsk - (modq - *it_pq) : *it_pq;
	    }
	}
    }


    //////////////// print fcts ////////////////
    ////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::print_q(int i)
    {
      std::cout << "cq["<<i<<"]: " << this->__cq[i] << std::endl;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::print_b(int i)
    {
      std::cout << "cb["<<i<<"]: " << this->__cb[i] << std::endl;
    }

    template<class Pq, class Pb, class Psk, class parameters, class Pmt>
    void poly_digit<Pq, Pb, Psk, parameters, Pmt>::print_sk(int i)
    {
      std::cout << "csk["<<i<<"]: " << this->__csk[i] << std::endl;
    }
  }

  // structure containing a two polynomials in base (q,b,msk), with degree at most d-1
  // represented by a vector of polynomial digits (see structure poly_digit above)
  template<class poly_digit, size_t d, typename Pmt>
  struct poly_hpr
  {
    typedef typename poly_digit::_Pq Pq;
    typedef typename poly_digit::_Pb Pb;
    typedef typename poly_digit::_Psk Psk;
    static constexpr size_t degree = Pq::degree;

    using main_type = typename Pq::value_type;
    using main_greater_type = typename Pq::greater_value_type;
    using main_signed_type = typename Pq::signed_value_type;


    typedef poly_digit _poly_digit;

    static constexpr uint __nmQ = poly_digit::__nmQ;
    static constexpr uint __nmB = poly_digit::__nmB;
    static constexpr uint __imB = poly_digit::__imB;
    static constexpr uint __imSK = poly_digit::__imSK;
    static constexpr uint __d = d;

    std::vector<poly_digit*> __digits; // vector of digits

    Pmt* __cmt;

    void allocate();
    void copy_from(const poly_hpr& other);
    poly_hpr();
    poly_hpr(const poly_hpr& other);
    poly_hpr& operator=(const poly_hpr& other);
    ~poly_hpr();

    typedef typename std::vector<poly_digit*>::iterator iterator;
    typedef typename std::vector<poly_digit*>::const_iterator const_iterator;

    iterator begin()
    {
      return __digits.begin();
    }
    iterator end()
    {
      return __digits.end();
    }
    const_iterator begin() const
    {
      return __digits.begin();
    }
    const_iterator end() const
    {
      return __digits.end();
    }

    void set0       ();
    void set0       (int);
    void set0_float       ();
    void set0_float       (int);
    void set        (int, const poly_hpr*, int);
    void set        (const poly_hpr*);
    void set_float        (int, const poly_hpr*, int);
    void set_float        (const poly_hpr*);

    typedef typename poly_digit::be_tool_set be_tool_set;
    void set_mat_set_conv(const be_tool_set*);

    bool operator==(const poly_hpr&) const;
    bool operator!=(const poly_hpr&) const;

    void add            (const poly_hpr*, int);
    void add            (int, const poly_hpr*, int, const poly_hpr*, int);
    void add_float            (const poly_hpr*, int);
    void add_float            (int, const poly_hpr*, int, const poly_hpr*, int);

    void sub            (int, const poly_hpr*, int);
    void sub            (int, const poly_hpr*, int, const poly_hpr*, int);
    void sub_float            (int, const poly_hpr*, int);
    void sub_float            (int, const poly_hpr*, int, const poly_hpr*, int);

    void mulDigit       (int, const poly_digit*, int);
    void mulDigit_float       (int, const poly_digit*, int);
    void mulDigit       (int, const poly_hpr*, int, const poly_digit*, int);
    void mulDigit_acc   (int, const poly_hpr*, int, const poly_digit*, int,
                         size_t s = d);
    void mulDigit_acc_float   (int, const poly_hpr*, int, const poly_digit*, int,
                               size_t s = d);

    void vtm_prod       (int, const poly_hpr*, int, const poly_hpr*, int);
    void vtm_prod_float       (int, const poly_hpr*, int, const poly_hpr*, int);

    void rotate         (size_t);

    void inv_ntt        (int);
    void inv_ntt        ();
    void inv_ntt_float        (int);
    void inv_ntt_float        ();
    void ntt            (int);
    void ntt            ();
    void ntt_float            (int);
    void ntt_float            ();

    void fbe            (int);
    void fbe            ();
    void fbe_float            (int);
    void fbe_float            ();

    void carry_propagation_1side      (int, size_t, size_t, bool exact = false);
    void carry_propagation_1side_float (int, size_t, size_t, bool exact = false);
    void carry_propagation      (size_t, size_t);
    void carry_propagation_float (size_t, size_t);

    void fast_modq ();
    void fast_modq_float ();

    void SmMRq();

    void print(int i = 0);
  };

  //////////////// constructor ////////////////
  /////////////////////////////////////////////

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::allocate()
  {
    for (size_t i = 0; i < d; ++i)
      this->__digits.push_back(new poly_digit());
    this->__cmt = commons::mem_manag::alloc_aligned<Pmt, 32>(2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::copy_from
  (const poly_hpr<poly_digit, d, Pmt> &other)
  {
    for (size_t i = 0; i < d; i++)
      __digits[i][0] = other.__digits[i][0];
    __cmt[0] = other.__cmt[0];
    __cmt[1] = other.__cmt[1];
  }  

  template<class poly_digit, size_t d, typename Pmt>
  poly_hpr<poly_digit, d, Pmt>::poly_hpr()
  {
    allocate();
  }

  template<class poly_digit, size_t d, typename Pmt>
  poly_hpr<poly_digit, d, Pmt>::poly_hpr
  (const poly_hpr<poly_digit, d, Pmt>& other)
  {
    allocate();
    copy_from(other);
  }

  template<class poly_digit, size_t d, typename Pmt>
  poly_hpr<poly_digit, d, Pmt>&
  poly_hpr<poly_digit, d, Pmt>::operator=
  (const poly_hpr<poly_digit, d, Pmt>& other)
  {
    copy_from(other);
    return *this;
  }

  //////////////// destructor ////////////////
  ////////////////////////////////////////////

  template<class poly_digit, size_t d, typename Pmt>
  poly_hpr<poly_digit, d, Pmt>::~poly_hpr()
  {
    for(iterator it = this->begin(); it != this->end(); ++it)
      delete *it;
    this->__digits.clear();
    commons::mem_manag::free_aligned<Pmt>(2, __cmt);
  }

  //////////////// setters ////////////////
  /////////////////////////////////////////

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set0()
  {
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->set0();
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set0(int c)
  {
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->set0(c);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set0_float()
  {
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->set0_float();
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set0_float(int c)
  {
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->set0_float(c);
  }


  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set(int i,
                                         const poly_hpr<poly_digit, d, Pmt> *other, int j)
  {
    const_iterator it_o = other->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      (*it)->set(i, *it_o, j);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set(const poly_hpr<poly_digit, d, Pmt>
                                         *other)
  {
    const_iterator it_o = other->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      (*it)->set(*it_o);
  }
  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set_float(int i,
					       const poly_hpr<poly_digit, d, Pmt> *other, int j)
  {
    const_iterator it_o = other->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      (*it)->set_float(i, *it_o, j);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set_float(const poly_hpr<poly_digit, d, Pmt>
					       *other)
  {
    const_iterator it_o = other->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      (*it)->set_float(*it_o);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::set_mat_set_conv(const be_tool_set* ts)
  {
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->set_mat_set_conv(ts);
  }

  //////////////// operations ////////////////
  ////////////////////////////////////////////

  template<class poly_digit, size_t d, typename Pmt>
  bool poly_hpr<poly_digit, d, Pmt>::operator==(const poly_hpr<poly_digit, d, Pmt>
						&other) const
  {
    const_iterator it_o = other.begin();
    for(const_iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      if (**it != **it_o)
        return false;
    return true;
  }

  template<class poly_digit, size_t d, typename Pmt>
  bool poly_hpr<poly_digit, d, Pmt>::operator!=(const poly_hpr<poly_digit, d, Pmt>
						&other) const
  {
    const_iterator it_o = other.begin();
    for(const_iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      if (**it != **it_o)
        return true;
    return false;
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::add(const poly_hpr<poly_digit, d, Pmt>
                                         *other, int i)
  {
    // addition digit by digit for component c_i

    const_iterator it_o = other->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      (*it)->add(i, *it_o, i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::add(int i0,
                                         const poly_hpr<poly_digit, d, Pmt> *other1, int i1,
                                         const poly_hpr<poly_digit, d, Pmt> *other2, int i2)
  {
    // addition digit by digit for this.c_i0 = other1.i1 + other2.c_i2

    const_iterator it_o1 = other1->begin(), it_o2 = other2->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1, ++it_o2)
      (*it)->add(i0, *it_o1, i1, *it_o2, i2);
  }
  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::add_float(const poly_hpr<poly_digit, d, Pmt>
					       *other, int i)
  {
    // addition digit by digit for component c_i

    const_iterator it_o = other->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      (*it)->add_float(i, *it_o, i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::add_float(int i0,
					       const poly_hpr<poly_digit, d, Pmt> *other1, int i1,
					       const poly_hpr<poly_digit, d, Pmt> *other2, int i2)
  {
    // addition digit by digit for this.c_i0 = other1.i1 + other2.c_i2

    const_iterator it_o1 = other1->begin(), it_o2 = other2->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1, ++it_o2)
      (*it)->add_float(i0, *it_o1, i1, *it_o2, i2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::sub(int i0,
                                         const poly_hpr<poly_digit, d, Pmt> *other1, int i1,
                                         const poly_hpr<poly_digit, d, Pmt> *other2, int i2)
  {
    // subtraction digit by digit for this.c_i0 = other1.i1 - other2.c_i2

    const_iterator it_o1 = other1->begin(), it_o2 = other2->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1, ++it_o2)
      (*it)->sub(i0, *it_o1, i1, *it_o2, i2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::sub_float(int i0,
					       const poly_hpr<poly_digit, d, Pmt> *other1, int i1,
					       const poly_hpr<poly_digit, d, Pmt> *other2, int i2)
  {
    // subtraction digit by digit for this.c_i0 = other1.i1 - other2.c_i2

    const_iterator it_o1 = other1->begin(), it_o2 = other2->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1, ++it_o2)
      (*it)->sub_float(i0, *it_o1, i1, *it_o2, i2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::sub(int i0,
                                         const poly_hpr<poly_digit, d, Pmt> *other1, int i1)
  {
    // subtraction digit by digit for this.c_i0 = this.i0 - other1.c_i1

    const_iterator it_o1 = other1->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1)
      (*it)->sub(i0, *it_o1, i1);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::sub_float(int i0,
					       const poly_hpr<poly_digit, d, Pmt> *other1, int i1)
  {
    // subtraction digit by digit for this.c_i0 = this.i0 - other1.c_i1

    const_iterator it_o1 = other1->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1)
      (*it)->sub_float(i0, *it_o1, i1);
  }


  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::mulDigit(int i0, const poly_digit *dig,
					      int i1)
  {
    // product by a single-digit element dig: this.c_i0 *= dig.c_i1

    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->mul(i0, dig, i1);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::mulDigit_float(int i0, const poly_digit *dig,
						    int i1)
  {
    // product by a single-digit element dig: this.c_i0 *= dig.c_i1

    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->mul_float(i0, dig, i1);
  }


  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::mulDigit(int i0, const poly_hpr *other1,
					      int i1, const poly_digit *dig, int i2)
  {
    // product by a single-digit element dig: this.c_i0 = other1.c_i1 * dig.c_i1
    const_iterator it_o = other1->begin();
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o)
      (*it)->mul(i0, *it_o, i1, dig, i2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::mulDigit_acc(int i0, const poly_hpr *other1,
						  int i1, const poly_digit *dig, int i2, size_t s)
  {
    // product by a single-digit element dig: this.c_i0 = other1.c_i1 * dig.c_i1
    // product on the first s digits of other1 (degree 0,..., s-1)
    // accumulates the result in this.c_i0
    const_iterator it_o1 = other1->begin();
    iterator it = this->begin();
    for(size_t is = 0; is != s; ++it, ++it_o1, ++is)
      (*it)->mul_acc(i0, *it_o1, i1, dig, i2);
  }
  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::mulDigit_acc_float(int i0,
							const poly_hpr *other1, int i1, const poly_digit *dig, int i2, size_t s)
  {
    // product by a single-digit element dig: this.c_i0 = other1.c_i1 * dig.c_i1
    // product on the first s digits of other1 (degree 0,..., s-1)
    // accumulates the result in this.c_i0
    const_iterator it_o1 = other1->begin();
    iterator it = this->begin();
    for(size_t is = 0; is != s; ++it, ++it_o1, ++is)
      (*it)->mul_acc_float(i0, *it_o1, i1, dig, i2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::vtm_prod(int i0, const poly_hpr *other1,
					      int i1, const poly_hpr *other2, int i2)
  {
    // noting this.c_i0 = (t_(0), ..., t_(d-1)) and other1.c_i1 = (o_(0), ..., o_(d-1))
    // computes the vector matrix product
    // (t_(0), ..., t_(d-1)) * [ o_(d-1)  0         .....       0 ]
    //                         [ o_(d-2)  o_(d-1)   .....       0 ]
    //                         [  ...                       ..... ]
    //                         [ o_(0)    o_(1)     ..... o_(d-1) ]
    // accumulates the result in other2.c_i2
    const_iterator it_o1_main = other1->begin();
    const_iterator it_o2_main = std::prev(other2->end());
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1_main)
      {
        const_iterator it_o1 = it_o1_main, it_o2 = it_o2_main;
        (*it)->mul(i0, *it_o1++, i1, *it_o2--, i2);
        for(; it_o1 != other1->end(); ++it_o1, --it_o2)
          (*it)->mul_acc(i0, *it_o1, i1, *it_o2, i2);
      }
  }
  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::vtm_prod_float(int i0,
						    const poly_hpr *other1, int i1, const poly_hpr *other2, int i2)
  {
    // noting this.c_i0 = (t_(0), ..., t_(d-1)) and other1.c_i1 = (o_(0), ..., o_(d-1))
    // computes the vector matrix product
    // (t_(0), ..., t_(d-1)) * [ o_(d-1)  0         .....       0 ]
    //                         [ o_(d-2)  o_(d-1)   .....       0 ]
    //                         [  ...                       ..... ]
    //                         [ o_(0)    o_(1)     ..... o_(d-1) ]
    // accumulates the result in other2.c_i2
    const_iterator it_o1_main = other1->begin();
    const_iterator it_o2_main = std::prev(other2->end());
    for(iterator it = this->begin(); it != this->end(); ++it, ++it_o1_main)
      {
        const_iterator it_o1 = it_o1_main, it_o2 = it_o2_main;
        (*it)->mul_float(i0, *it_o1++, i1, *it_o2--, i2);
        for(; it_o1 != other1->end(); ++it_o1, --it_o2)
          (*it)->mul_acc_float(i0, *it_o1, i1, *it_o2, i2);
      }
  }

  template<class poly_digit, size_t d, typename Pmt>
  inline void poly_hpr<poly_digit, d, Pmt>::rotate(size_t s)
  {
    // left-rotation by s positions: (t_(0), ..., t_(d-1)) becomes (t_(s), t_(s+1), ..., t_(s-1))
    std::rotate(this->begin(), this->begin()+s, this->end());
  }

  //////////////// ntt related ////////////////
  /////////////////////////////////////////////

  template<class poly_digit, size_t d, typename Pmt>
  inline void poly_hpr<poly_digit, d, Pmt>::inv_ntt(int i)
  {
    // digit-wise inv_ntt on component this.c_i
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->inv_ntt(i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::inv_ntt()
  {
    // digit-wise inv_ntt on components this.c_0 and this_c_1
    this->inv_ntt(0);
    this->inv_ntt(1);
  }

  template<class poly_digit, size_t d, typename Pmt>
  inline void poly_hpr<poly_digit, d, Pmt>::inv_ntt_float(int i)
  {
    // digit-wise inv_ntt on component this.c_i
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->inv_ntt_float(i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::inv_ntt_float()
  {
    // digit-wise inv_ntt on components this.c_0 and this_c_1
    this->inv_ntt_float(0);
    this->inv_ntt_float(1);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::ntt(int i)
  {
    // digit-wise ntt on component this.c_i
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->ntt(i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::ntt()
  {
    // digit-wise ntt on components this.c_0 and this_c_1
    this->ntt(0);
    this->ntt(1);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::ntt_float(int i)
  {
    // digit-wise ntt on component this.c_i
    for(iterator it = this->begin(); it != this->end(); ++it)
      (*it)->ntt_float(i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::ntt_float()
  {
    // digit-wise ntt on components this.c_0 and this_c_1
    this->ntt_float(0);
    this->ntt_float(1);
  }

  template<class poly_digit, size_t d, typename Pmt>
  inline void poly_hpr<poly_digit, d, Pmt>::fast_modq()
  {
    // sets residues modulo b*msk of last digit to 0 (useless when operating modulo q^d)
    this->__digits[d-1]->setb0();
    this->__digits[d-1]->setsk0();
  }

  template<class poly_digit, size_t d, typename Pmt>
  inline void poly_hpr<poly_digit, d, Pmt>::fast_modq_float()
  {
    // sets residues modulo b*msk of last digit to 0 (useless when operating modulo q^d)
    this->__digits[d-1]->setb0();
  }

  //////////////// base conversions ////////////////
  //////////////////////////////////////////////////

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::fbe(int i)
  {
    // digit-wise fast base extension of component this.c_i
    for(iterator it = this->begin(); it != std::prev(this->end()); ++it)
      (*it)->fbe(i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::fbe()
  {
    // digit-wise fast base extension of components this.c_0 and this.c_1
    for(iterator it = this->begin(); it != std::prev(this->end()); ++it)
      (*it)->fbe();
  }
  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::fbe_float(int i)
  {
    // digit-wise fast base extension of component this.c_i
    for(iterator it = this->begin(); it != std::prev(this->end()); ++it)
      (*it)->fbe_float(i);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::fbe_float()
  {
    // digit-wise fast base extension of components this.c_0 and this.c_1
    for(iterator it = this->begin(); it != std::prev(this->end()); ++it)
      (*it)->fbe_float();
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::SmMRq()
  {
    // fast reduction modulo q of most significant digit
    // by using a Montgomery approach with special modulus 2^BS_MT (mtilde)
    // concretely, extends this.c_0_(d-1) mod q and this.c_1_(d-1) mod q
    // to base b*msk*mtilde, and reduces (mod q) the result in base b*msk via residue mod mtilde

    std::array<main_signed_type, degree> acc;

    ////
    iterator it_last = std::prev(this->end());

    (*it_last)->__cb[0] = 0;
    (*it_last)->__csk[0] = 0;
    this->__cmt[0].set0();

    (*it_last)->__cq[0] = (*it_last)->__cq[0] *
      (*it_last)->__be_tool_set->_xi_q_mt[0];

    this->__cmt[0].template fast_base_conversion<Pq>((*it_last)->__cq[0]);
    this->__cmt[0].get_modc(acc);

    (*it_last)->__be_tool_set->SmMRq((*it_last)->__cb[0], (*it_last)->__csk[0],
                                     (*it_last)->__cq[0], acc);

    (*it_last)->__cq[0] = (*it_last)->__cq[0] *
      (*it_last)->__be_tool_set->_xi_q_mt[1];

    ////
    (*it_last)->__cb[1] = 0;
    (*it_last)->__csk[1] = 0;
    this->__cmt[1].set0();

    (*it_last)->__cq[1] = (*it_last)->__cq[1] *
      (*it_last)->__be_tool_set->_xi_q_mt[0];

    this->__cmt[1].template fast_base_conversion<Pq>((*it_last)->__cq[1]);
    this->__cmt[1].get_modc(acc);

    (*it_last)->__be_tool_set->SmMRq((*it_last)->__cb[1], (*it_last)->__csk[1],
                                     (*it_last)->__cq[1], acc);

    (*it_last)->__cq[1] = (*it_last)->__cq[1] *
      (*it_last)->__be_tool_set->_xi_q_mt[1];
  }

  //////////////// carry propagation ////////////////
  ///////////////////////////////////////////////////

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::carry_propagation_1side(int c, size_t i0,
							     size_t i1, bool exact)
  {
    // propagates the carries from digit number i0 to digit number i1-1
    // for compoent this.c_c
    // 2 steps: fast base extension (or exact using gmp if exact=true) from q to b*msk in same digit,
    // then exact conversion (shenoy-kumaresan) from b*msk to q for next digit.

    Pb* tmp_polb = commons::mem_manag::alloc_aligned<Pb, 32>(2);
    Psk* tmp_polsk = commons::mem_manag::alloc_aligned<Psk, 32>(2);
    Pb *tmp_polb2 = commons::mem_manag::alloc_aligned<Pb, 32>(1);
    Psk *tmp_polsk2 = commons::mem_manag::alloc_aligned<Psk, 32>(1);

    for(iterator it = this->begin() + i0; it != this->begin() + i1; ++it)
      {
        tmp_polb[0] = (*it)->__cb[c];
        tmp_polsk[0] = (*it)->__csk[c];

        // p -> Bsk: carry comp.
        if (exact)
          (*it)->exact_conv_gmp_q_to_bsk(c);
        else
          (*it)->fbe(c);

        tmp_polb[1] = (*it)->__cb[c];
        tmp_polsk[1] = (*it)->__csk[c];

        *tmp_polb2 = tmp_polb[0]-tmp_polb[1];
        *tmp_polsk2 = tmp_polsk[0]-tmp_polsk[1];
        (*it)->__cb[c] = (*tmp_polb2) * (*it)->__be_tool_set->_qinv_modB[0];
        (*it)->__csk[c] = (*tmp_polsk2) * (*it)->__be_tool_set->_qinv_modmsk[0];

        // Bsk -> p: carry propag.
        (*it)->skbe(c, *(it+1));

        (*(it+1))->__cb[c] = (*(it+1))->__cb[c] + (*it)->__cb[c];
        (*(it+1))->__csk[c] = (*(it+1))->__csk[c] + (*it)->__csk[c];

        //
        (*it)->__cb[c] = tmp_polb[1];
        (*it)->__csk[c] = tmp_polsk[1];
      }

    commons::mem_manag::free_aligned<Pb>(2, tmp_polb);
    commons::mem_manag::free_aligned<Psk>(2, tmp_polsk);
    commons::mem_manag::free_aligned<Pb>(1, tmp_polb2);
    commons::mem_manag::free_aligned<Psk>(1, tmp_polsk2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::carry_propagation_1side_float(int c,
								   size_t i0, size_t i1, bool exact)
  {
    // propagates the carries from digit number i0 to digit number i1-1
    // for compoent this.c_c
    // 2 steps: fast base extension (or exact using gmp if exact=true) from q to b*msk in same digit,
    // then exact conversion (shenoy-kumaresan) from b*msk to q for next digit.

    Pb* tmp_polb = commons::mem_manag::alloc_aligned<Pb, 32>(2);
    Psk* tmp_polsk = commons::mem_manag::alloc_aligned<Psk, 32>(2);
    Pb *tmp_polb2 = commons::mem_manag::alloc_aligned<Pb, 32>(1);
    Psk *tmp_polsk2 = commons::mem_manag::alloc_aligned<Psk, 32>(1);

    for(iterator it = this->begin() + i0; it != this->begin() + i1; ++it)
      {
        tmp_polb[0] = (*it)->__cb[c];
        tmp_polsk[0] = (*it)->__csk[c];

        // p -> Bsk: carry comp.
        if (exact)
          (*it)->exact_conv_gmp_q_to_bsk(c);
        else
          (*it)->fbe(c);

        tmp_polb[1] = (*it)->__cb[c];
        tmp_polsk[1] = (*it)->__csk[c];

        *tmp_polb2 = tmp_polb[0]-tmp_polb[1];
        *tmp_polsk2 = tmp_polsk[0]-tmp_polsk[1];
        (*it)->__cb[c] = (*tmp_polb2) * (*it)->__be_tool_set->_qinv_modB[0];
        (*it)->__csk[c] = (*tmp_polsk2) * (*it)->__be_tool_set->_qinv_modmsk[0];

        // Bsk -> p: carry propag.
        (*it)->float_back(c, *(it+1));

        (*(it+1))->__cb[c] = (*(it+1))->__cb[c] + (*it)->__cb[c];
        (*(it+1))->__csk[c] = (*(it+1))->__csk[c] + (*it)->__csk[c];

        //
        (*it)->__cb[c] = tmp_polb[1];
        (*it)->__csk[c] = tmp_polsk[1];
      }

    commons::mem_manag::free_aligned<Pb>(2, tmp_polb);
    commons::mem_manag::free_aligned<Psk>(2, tmp_polsk);
    commons::mem_manag::free_aligned<Pb>(1, tmp_polb2);
    commons::mem_manag::free_aligned<Psk>(1, tmp_polsk2);
  }


  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::carry_propagation(size_t i0, size_t i1)
  {
    // ditto for both components this.c_i0, this.c_i1

    Pb* tmp_polb = commons::mem_manag::alloc_aligned<Pb, 32>(4);
    Psk* tmp_polsk = commons::mem_manag::alloc_aligned<Psk, 32>(4);
    Pb *tmp_polb2 = commons::mem_manag::alloc_aligned<Pb, 32>(1);
    Psk *tmp_polsk2 = commons::mem_manag::alloc_aligned<Psk, 32>(1);

    for(iterator it = this->begin() + i0; it != this->begin() + i1; ++it)
      {
        tmp_polb[0] = (*it)->__cb[0];
        tmp_polb[1] = (*it)->__cb[1];

        tmp_polsk[0] = (*it)->__csk[0];
        tmp_polsk[1] = (*it)->__csk[1];

        // p -> Bsk: carry comp.
        (*it)->fbe();

        tmp_polb[2] = (*it)->__cb[0];
        tmp_polb[3] = (*it)->__cb[1];

        tmp_polsk[2] = (*it)->__csk[0];
        tmp_polsk[3] = (*it)->__csk[1];

        *tmp_polb2 = tmp_polb[0]-tmp_polb[2];
        (*it)->__cb[0] = (*tmp_polb2) * (*it)->__be_tool_set->_qinv_modB[0];
        *tmp_polb2 = tmp_polb[1]-tmp_polb[3];
        (*it)->__cb[1] = (*tmp_polb2) * (*it)->__be_tool_set->_qinv_modB[0];

        *tmp_polsk2 = tmp_polsk[0]-tmp_polsk[2];
        (*it)->__csk[0] = (*tmp_polsk2) * (*it)->__be_tool_set->_qinv_modmsk[0];
        *tmp_polsk2 = tmp_polsk[1]-tmp_polsk[3];
        (*it)->__csk[1] = (*tmp_polsk2) * (*it)->__be_tool_set->_qinv_modmsk[0];

        // Bsk -> p: carry propag.
        (*it)->skbe(*(it+1));

        (*(it+1))->__cb[0] = (*(it+1))->__cb[0] + (*it)->__cb[0];
        (*(it+1))->__cb[1] = (*(it+1))->__cb[1] + (*it)->__cb[1];

        (*(it+1))->__csk[0] = (*(it+1))->__csk[0] + (*it)->__csk[0];
        (*(it+1))->__csk[1] = (*(it+1))->__csk[1] + (*it)->__csk[1];

        //
        (*it)->__cb[0] = tmp_polb[2];
        (*it)->__cb[1] = tmp_polb[3];

        (*it)->__csk[0] = tmp_polsk[2];
        (*it)->__csk[1] = tmp_polsk[3];
      }

    commons::mem_manag::free_aligned<Pb>(4, tmp_polb);
    commons::mem_manag::free_aligned<Psk>(4, tmp_polsk);
    commons::mem_manag::free_aligned<Pb>(1, tmp_polb2);
    commons::mem_manag::free_aligned<Psk>(1, tmp_polsk2);
  }

  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::carry_propagation_float(size_t i0, size_t i1)
  {
    // ditto for both components this.c_i0, this.c_i1

    Pb* tmp_polb = commons::mem_manag::alloc_aligned<Pb, 32>(4);
    Psk* tmp_polsk = commons::mem_manag::alloc_aligned<Psk, 32>(4);
    Pb *tmp_polb2 = commons::mem_manag::alloc_aligned<Pb, 32>(1);
    Psk *tmp_polsk2 = commons::mem_manag::alloc_aligned<Psk, 32>(1);

    for(iterator it = this->begin() + i0; it != this->begin() + i1; ++it)
      {
        tmp_polb[0] = (*it)->__cb[0];
        tmp_polb[1] = (*it)->__cb[1];

        tmp_polsk[0] = (*it)->__csk[0];
        tmp_polsk[1] = (*it)->__csk[1];

        // p -> Bsk: carry comp.
        (*it)->fbe();

        tmp_polb[2] = (*it)->__cb[0];
        tmp_polb[3] = (*it)->__cb[1];

        tmp_polsk[2] = (*it)->__csk[0];
        tmp_polsk[3] = (*it)->__csk[1];

        *tmp_polb2 = tmp_polb[0]-tmp_polb[2];
        (*it)->__cb[0] = (*tmp_polb2) * (*it)->__be_tool_set->_qinv_modB[0];
        *tmp_polb2 = tmp_polb[1]-tmp_polb[3];
        (*it)->__cb[1] = (*tmp_polb2) * (*it)->__be_tool_set->_qinv_modB[0];

        *tmp_polsk2 = tmp_polsk[0]-tmp_polsk[2];
        (*it)->__csk[0] = (*tmp_polsk2) * (*it)->__be_tool_set->_qinv_modmsk[0];
        *tmp_polsk2 = tmp_polsk[1]-tmp_polsk[3];
        (*it)->__csk[1] = (*tmp_polsk2) * (*it)->__be_tool_set->_qinv_modmsk[0];

        // Bsk -> p: carry propag.
        //(*it)->skbe(*(it+1));
        (*it)->float_back (*(it+1));

        (*(it+1))->__cb[0] = (*(it+1))->__cb[0] + (*it)->__cb[0];
        (*(it+1))->__cb[1] = (*(it+1))->__cb[1] + (*it)->__cb[1];

        (*(it+1))->__csk[0] = (*(it+1))->__csk[0] + (*it)->__csk[0];
        (*(it+1))->__csk[1] = (*(it+1))->__csk[1] + (*it)->__csk[1];

        //
        (*it)->__cb[0] = tmp_polb[2];
        (*it)->__cb[1] = tmp_polb[3];

        (*it)->__csk[0] = tmp_polsk[2];
        (*it)->__csk[1] = tmp_polsk[3];
      }

    commons::mem_manag::free_aligned<Pb>(4, tmp_polb);
    commons::mem_manag::free_aligned<Psk>(4, tmp_polsk);
    commons::mem_manag::free_aligned<Pb>(1, tmp_polb2);
    commons::mem_manag::free_aligned<Psk>(1, tmp_polsk2);
  }


  template<class poly_digit, size_t d, typename Pmt>
  void poly_hpr<poly_digit, d, Pmt>::print(int i)
  {
    // prints component this.c_i

    int cpt = 0;
    std::cout<<"component #"<<i<<":"<<std::endl<<std::endl;
    for(const_iterator it = this->begin(); it != this->end(); ++it, cpt++)
      {
        std::cout<<"digit #"<<cpt<<":"<<std::endl<<std::endl;
        (*it)->print_q(i);
        (*it)->print_b(i);
        (*it)->print_sk(i);
        std::cout<<std::endl<<std::endl;
      }
  }

  namespace context
  {
    //========== General context - structure gathering precomputations ==========//
    /*****************************************************************************/

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    struct Context
    {
      bool __red, __redhps;

      typedef typename digit::poly_digit<Pq, Pb, Psk, parameters, Pmt> poly_digit;
      typedef typename hpr::poly_hpr<poly_digit, d, Pmt> poly_hpr;
      using main_type = typename Pq::value_type;
      using greater_main_type = typename Pq::greater_value_type;
      using main_signed_type = typename Pq::signed_value_type;
      typedef commons::raw_pol::PolT<parameters, Pq::degree> Pt;

      static constexpr uint __nmQ = Pq::nmoduli;  // number of moduli in base q
      static constexpr uint __nmB = Pb::nmoduli;  // number of moduli in base b
      static constexpr uint __imQ =
	Pq::indmod0;  // index of first modulus of base q in nfl::params<main_type>::P
      static constexpr uint __imB =
	Pb::indmod0;  // index of first modulus of base b in nfl::params<main_type>::P
      static constexpr uint __imSK =
	Psk::indmod0; // index of modulus m_sk in nfl::params<main_type>::P
      static constexpr uint __d = d;

      typedef typename
      base_conversions::base_conv_tool_set<Pq, Pb, Psk, parameters, Pmt> be_tool_set;
      be_tool_set* __be_tool_set;

      typedef commons::key::key_chain<poly_hpr, parameters> key_t;
      key_t* __keys;

      // encryption
      poly_hpr* __Delta;

      // Multiplication
      poly_digit* __t;
      poly_digit* __dig_decomp;
      poly_hpr* __coeff2;
      poly_digit* __carry;


      // noise size
      mpz_class _p, _qbsk, _q;
      std::vector<mpz_class> _pi;
      std::vector<mpz_class> _pi_ov_p;

      // Relinearisation key
      std::vector<poly_hpr*> __rlk;

      using iterator = typename std::vector<poly_hpr*>::iterator;
      using const_iterator = typename std::vector<poly_hpr*>::const_iterator;

      iterator begin()
      {
	return __rlk.begin();
      }
      iterator end()
      {
	return __rlk.end();
      }

      const_iterator begin() const
      {
	return __rlk.begin();
      }
      const_iterator end() const
      {
	return __rlk.end();
      }

      void allocate();
      void copy_from(const Context &other);

      // constructor
      Context(bool red = true, bool redhps = true);
      Context(const Context &other);
      Context& operator=(const Context &other);

      // destructor
      ~Context();

      // print fcts (debug)
      void print_baseq();
      void print_baseb();
      void print_basesk();
      void print_basegt();

      void set_precomp(key_t *keys = nullptr);
      void set_precomp_relin();

      // encryption
      void hidden_embedding(poly_hpr*, const poly_digit*, size_t, bool);
      void encrypt(poly_hpr*, const Pt*) const;
      void postproc_encrypt(poly_hpr*) const;

      // decryption
      void preproc_decrypt(poly_hpr*) const;
      void decrypt(Pt*, poly_hpr*) const;

      // addition
      void addition(poly_hpr*, const poly_hpr*, const poly_hpr*) const;
      void addition(poly_hpr*, const poly_hpr*) const;

      // multiplication
      void multiplication(poly_hpr*, poly_hpr*);
      void multiplication_float(poly_hpr*, poly_hpr*);

      // relinearisation
      void relinearisation(poly_hpr*);
      void relinearisation_float(poly_hpr*);

      // noise size measuring function
      size_t noise_size(const poly_hpr*, const Pt*, bool) const;
    };

    //////////////// constructor ////////////////
    /////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::allocate()
    {
      __be_tool_set = new be_tool_set();
      __keys = new key_t();
      __t = new poly_digit();
      __Delta = new poly_hpr();
      __dig_decomp = new poly_digit();
      __coeff2 = new poly_hpr();
      __carry = new poly_digit();

      for(uint i = 0 ; i < __d*__nmQ; ++i)
	this->__rlk.push_back(new poly_hpr());
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::copy_from
    (const Context<Pq, Pb, Psk, Pmt, parameters, d> &other)
    {
      __red = other.__red;
      __redhps = other.__redhps;
      *__be_tool_set = *other.__be_tool_set;
      *__keys = *other.__keys;
      *__Delta = *other.__Delta;
      *__t = *other.__t;
      *__dig_decomp = *other.__dig_decomp;
      *__coeff2 = *other.__coeff2;
      *__carry = *other.__carry;
      _p = other._p;
      _qbsk = other._qbsk;
      _q = other._q;
      _pi = other._pi;
      _pi_ov_p = other._pi_ov_p;
      for (uint i = 0; i < __d * __nmQ; i++)
	__rlk[i] = other.__rlk[i];
    }
    
    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    Context<Pq, Pb, Psk, Pmt, parameters, d>::Context(bool red, bool redhps)
      : __red(red), __redhps(redhps)
    {
      allocate();
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    Context<Pq, Pb, Psk, Pmt, parameters, d>::Context
    (const Context<Pq, Pb, Psk, Pmt, parameters, d>& other)
    {
      allocate();
      copy_from(other);
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    Context<Pq, Pb, Psk, Pmt, parameters, d>&
    Context<Pq, Pb, Psk, Pmt, parameters, d>::operator=
    (const Context<Pq, Pb, Psk, Pmt, parameters, d>& other)
    {
      copy_from(other);
      return *this;
    }
    
    //////////////// destructor ////////////////
    ////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    Context<Pq, Pb, Psk, Pmt, parameters, d>::~Context()
    {
      delete __be_tool_set;
      delete __keys;
      delete __t;
      delete __Delta;
      delete __dig_decomp;
      delete __coeff2;
      delete __carry;

      for(iterator it = this->begin(); it != this->end(); ++it)
	delete *it;
      this->__rlk.clear();

      _pi.clear();
      _pi_ov_p.clear();
    }

    //////////////// setters ////////////////
    /////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::set_precomp(key_t* keys)
    {
      // precomputation of useful data
      constexpr size_t degree = Pq::degree;

      this->__be_tool_set->generate();
      
      if (keys == nullptr)
	{
	  this->__keys->generate(__be_tool_set);
	}
      else
	{
	  *this->__keys = *keys;
	}

      this->__t->set_mat_set_conv(this->__be_tool_set);

      this->__Delta->set_mat_set_conv(this->__be_tool_set);
      this->__coeff2->set_mat_set_conv(this->__be_tool_set);
      this->__carry->set_mat_set_conv(this->__be_tool_set);
      _q = 1;
      {
	mpz_class tmp;
	for (uint cmQ = 0; cmQ < __nmQ; cmQ++)
	  {
	    mpz_set_ui(tmp.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]);
	    _q = _q*tmp;

	    _pi.push_back(tmp);
	  }

	mpz_pow_ui(_p.get_mpz_t(), _q.get_mpz_t(), d);

	_qbsk = _q;

	for (uint cmB = 0; cmB < __nmB; cmB++)
	  {
	    mpz_set_ui(tmp.get_mpz_t(), Pb::base.params.P[__imB+cmB]);
	    _qbsk = _qbsk*tmp;
	    _pi.push_back(tmp);
	  }

	mpz_set_ui(tmp.get_mpz_t(), Psk::base.params.P[__imSK]);
	_pi.push_back(tmp);
	_qbsk = _qbsk*tmp;

	for (uint cm = 0; cm < __nmQ+__nmB+1; cm++)
	  {
	    tmp = _qbsk / _pi[cm];
	    mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(), _pi[cm].get_mpz_t());
	    _pi_ov_p.push_back(tmp);
	  }
      }

      // carry prop. and encryption precomp.
      {
	// __Delta = digit containing floor(q/t) in ntt representation

	for(int i = 0; i < 2; i++)
	  {
	    std::fill(this->__t->__cq[i].begin(), this->__t->__cq[i].end(),
		      static_cast<main_type>(parameters::t));
	    std::fill(this->__t->__cb[i].begin(), this->__t->__cb[i].end(),
		      static_cast<main_type>(parameters::t));
	    std::fill(this->__t->__csk[i].begin(), this->__t->__csk[i].end(),
		      static_cast<main_type>(parameters::t));
	  }

	this->__Delta->set0();
	mpz_class Delta = _p;
	Delta = Delta>>parameters::_bs_t;

	for(auto it = this->__Delta->begin(); it != this->__Delta->end(); ++it)
	  {
	    mpz_class cur_delta = Delta % _q;
	    auto* itq1 = (*it)->__cq[0].begin();
	    for (uint cmQ = 0; cmQ < __nmQ; cmQ++, itq1 += degree)
	      {
		mpz_class tmp;
		mpz_set_ui(tmp.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]);
		tmp = cur_delta%tmp;
		std::fill(itq1, itq1+degree,
			  static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
	      }
	    auto* itb1 = (*it)->__cb[0].begin();
	    for (uint cmB = 0; cmB < __nmB; cmB++, itb1 += degree)
	      {
		mpz_class tmp, mod;
		mpz_set_ui(mod.get_mpz_t(), Pb::base.params.P[__imB+cmB]);
		tmp = cur_delta%mod;
		std::fill(itb1, itb1+degree,
			  static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
	      }

	    {
	      mpz_class tmp, mod;
	      mpz_set_ui(mod.get_mpz_t(), Psk::base.params.P[__imSK]);
	      tmp = cur_delta%mod;
	      auto* itsk = (*it)->__csk[0].begin();
	      std::fill(itsk, itsk+degree,
			static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
	    }

	    Delta = (Delta-cur_delta)/_q;
	  }
      } // end encryption precomp.
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::set_precomp_relin()
    {
      // Relin. precomp

      mpz_class tmp;
      poly_digit* tmp_dig = new poly_digit();

      this->__dig_decomp->set_mat_set_conv(this->__be_tool_set);
      tmp_dig->set_mat_set_conv(this->__be_tool_set);

      for(iterator it_rlk = this->begin(); it_rlk != this->end(); ++it_rlk)
	(*it_rlk)->set_mat_set_conv(this->__be_tool_set);

      for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
	{
	  // multiply s2 by q/q_slot
	  tmp_dig->set(__keys->_s2);
	  tmp_dig->inv_ntt(0);

	  mpz_div_ui(tmp.get_mpz_t(), _q.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]);
	  commons::tools::mul_mpz_class<Pq>(tmp_dig->__cq[0], tmp);
	  commons::tools::mul_mpz_class<Pb>(tmp_dig->__cb[0], tmp);
	  commons::tools::mul_mpz_class<Psk>(tmp_dig->__csk[0], tmp);

	  tmp_dig->ntt(0);

	  for(uint deg = 0; deg < __d; deg++)
	    hidden_embedding(__rlk[__nmQ*deg+cmQ], tmp_dig, deg, true);
	}
      delete tmp_dig;
    }

    //////////////// operations ////////////////
    ////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::hidden_embedding(poly_hpr *c,
								    const poly_digit *m, size_t slot, bool input_in_ntt)
    {
      // embeds m in c_1_[slot], and adds an encryption of 0 to c
      // if zeroing_msd == true, sets the digits from (slot+1) to d-1 to zero

      nfl::uniform unif;

      c->set0();
      c->__digits[0]->gen_noise(0);

      int  i=0;
      for(typename poly_hpr::iterator it = c->begin(); it != c->begin()+slot+1; ++it)
	{
	  (*it)->__cq[1] = unif;
	  if(it != c->end()-1)
	    (*it)->exact_conv_gmp_q_to_bsk(1);
	}

      c->ntt();
      c->mulDigit_acc(0, c, 1, __keys->_sk, 1); // c = (e - a*s, a)

      if(input_in_ntt == true)
	c->__digits[slot]->add(0, m, 0);

      c->inv_ntt(0);

      if(input_in_ntt == false)
	c->__digits[slot]->add(0, m, 0);

      // here, c = (e - a*s + m, a)

      c->carry_propagation_1side(0, 0, slot, true);
      c->ntt(0);

      c->fast_modq();
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::encrypt(poly_hpr *c,
							   const Pt *m) const
    {
      c->set0();

      // copying message in MSD:
      base_conversions::copy_paste_conversion<Pq, Pt, typename Pt::type_t>
	(c->__digits[d-1]->__cq[0], *m, parameters::t_t);
      base_conversions::copy_paste_conversion<Pb, Pt, typename Pt::type_t>
	(c->__digits[d-1]->__cb[0], *m, parameters::t_t);
      base_conversions::copy_paste_conversion<Psk, Pt, typename Pt::type_t>
	(c->__digits[d-1]->__csk[0], *m, parameters::t_t);

      c->__digits[d-1]->ntt(0);

      // mutiplying it to precomputed polynomial Delta:
      if constexpr (parameters::_bs_t < 16)
		     {
	  c->__digits[d-1]->mul(0, __Delta->__digits[d-1], 0, c->__digits[d-1], 0);
	  auto it_end = std::prev(c->end());
	  for(auto it_c = c->begin(); it_c != std::prev(c->end()); ++it_c)
	    (*it_c)->set(0, *it_end, 0);
	}
      else
	  {
	    c->mulDigit(0, __Delta, 0, c->__digits[d-1], 0);
	  }

      // here, c = (Delta*m, 0)

      poly_digit* noise = new poly_digit();
      noise->set_mat_set_conv(this->__be_tool_set);

      noise->gen_noise(3); // (e0, e1)
      noise->ntt();
      c->__digits[0]->add(0, noise, 0);
      c->__digits[0]->add(1, noise, 1);

      // here, c = (e0 + Delta*m, e1)

      noise->set0();
      noise->gen_noise(0); // (e2, 0)
      noise->ntt(0);

      c->mulDigit_acc(0, __keys->_pk, 0, noise, 0);
      c->mulDigit_acc(1, __keys->_pk, 1, noise, 0);

      // here, c = (e2*p0 + e0 + Delta*m, e2*p1 + e1)

      delete noise;

    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::postproc_encrypt(
								    poly_hpr *c) const
    {
      //post-processing step, should be performed by client (?)
      c->inv_ntt();
      c->carry_propagation(0, d-1);
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::preproc_decrypt(
								   poly_hpr *p) const
    {
      // preprocessing step of decryption which can be performed by server
      // p->carry_propagation(d-2, d-1);
      p->__digits[d-1]->nttq();
      p->__digits[d-1]->__cq[0] = p->__digits[d
					      -1]->__cq[0]*__be_tool_set->_xi_q_gt[0];
      p->__digits[d-1]->__cq[1] = p->__digits[d
					      -1]->__cq[1]*__be_tool_set->_xi_q_gt[0];
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::decrypt(Pt *mdec,
							   poly_hpr *p) const
    {
      // classical gamma-correction style decyrption on most significant digit
      constexpr size_t degree = Pq::degree;

      poly_digit* c = p->__digits[d-1];

      c->__cq[0] = c->__cq[0] + c->__cq[1] * __keys->_sk->__cq[0];
      c->__cq[0].invntt_pow_invphi();

      auto* it_mdec = mdec->begin();
      for(uint i = 0; i < degree; i++)
	{
	  auto* it = c->__cq[0].begin()+i;
	  typename parameters::type_gamma_t acc = parameters::gamma_half;
	  for(uint cm = 0; cm < __nmQ; cm++, it += degree)
	    acc += (~(static_cast<typename parameters::type_gamma_t>(*it)))+1;

	  if constexpr (parameters::_bs_t == 8 || parameters::_bs_t == 16)
	      {
		*it_mdec++ = static_cast<typename Pt::type_t>(acc>>parameters::_bs_g);
	      }
	      else
	      {
		*it_mdec++ = static_cast<typename Pt::type_t>(acc>>parameters::_bs_g)
		  &parameters::mask_t;
	      }
	}
    }


    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::addition(poly_hpr *s,
							    const poly_hpr *a1, const poly_hpr *a2) const
    {
      s->add(a1, a2);
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    inline void Context<Pq, Pb, Psk, Pmt, parameters, d>::addition(poly_hpr *s,
								   const poly_hpr *a) const
    {
      s->add(a, 0);
      s->add(a, 1);
    }


    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::multiplication(poly_hpr *c1,
								  poly_hpr *c2)
    {
      {
	// reducing mod q

	c1->fast_modq();
	c2->fast_modq();

	if (__red)
	  {
	    if (!__redhps)
	      {
		c1->SmMRq();
		c2->SmMRq();
	      }
	    else
	      {
		c1->__digits[d-1]->float_forward_sk ();
		c2->__digits[d-1]->float_forward_sk ();
	      }
	  }
	else
	  {
	    c1->__digits[d-1]->fbe();
	    c2->__digits[d-1]->fbe();
	  }
      }

      {
	// multiply by t

	c1->mulDigit(0, __t, 0);
	c1->mulDigit(1, __t, 1);
      }

      {
	// ntt the input

	c1->ntt();
	c2->ntt();
      }

      {
	// karatsuba pattern

	this->__coeff2->set0();

	this->__coeff2->add(0, c1, 0, c1, 1); // coeff2 = (c10+c11, ***)
	this->__coeff2->add(1, c2, 0, c2, 1); // coeff2 = (c10+c11, c20+C21)

	this->__coeff2->vtm_prod(0, this->__coeff2, 0, this->__coeff2,
				 1); // coeff2 = ( (c10+c11)*(c20+C21), ***)

	c1->vtm_prod(0, c1, 0, c2, 0); // c1 = (c10*c20, c11)
	c1->vtm_prod(1, c1, 1, c2, 1); // c1 = (c10*c20, c11*c21)

	this->__coeff2->sub(0, c1, 0); // coeff2 = ( (c10+c11)*(c20+C21) - c10*c20, ***)
	this->__coeff2->sub(0, c1,
			    1); // coeff2 = ( (c10+c11)*(c20+C21) - c10*c20 - c11*c21, ***)

	this->__coeff2->set(1, c1,
			    1); // coeff2 = ( (c10+c11)*(c20+C21) - c10*c20 - c11*c21, c11*c21)
	c1->set(1, this->__coeff2,
		0); // c1 = (c10*c20, (c10+c11)*(c20+C21) - c10*c20 - c11*c21)

      }

      {
	// ending the division+rounding

	c1->__digits[0]->inv_ntt();
	// c1->__digits[1]->inv_ntt();
	// c1->carry_propagation(0, 1); // propagates first carry
	__carry->__cb[0] = c1->__digits[0]->__cb[0];
	__carry->__cb[1] = c1->__digits[0]->__cb[1];
	__carry->__csk[0] = c1->__digits[0]->__csk[0];
	__carry->__csk[1] = c1->__digits[0]->__csk[1];

	c1->__digits[0]->fbe();

	__carry->__cb[0] = __carry->__cb[0] - c1->__digits[0]->__cb[0];
	__carry->__cb[1] = __carry->__cb[1] - c1->__digits[0]->__cb[1];
	__carry->__cb[0] = __carry->__cb[0] * __carry->__be_tool_set->_qinv_modB[0];
	__carry->__cb[1] = __carry->__cb[1] * __carry->__be_tool_set->_qinv_modB[0];
	__carry->__csk[0] = __carry->__csk[0] - c1->__digits[0]->__csk[0];
	__carry->__csk[1] = __carry->__csk[1] - c1->__digits[0]->__csk[1];
	__carry->__csk[0] = __carry->__csk[0] * __carry->__be_tool_set->_qinv_modmsk[0];
	__carry->__csk[1] = __carry->__csk[1] * __carry->__be_tool_set->_qinv_modmsk[0];

	__carry->__cq[0].set(0);
	__carry->__cq[1].set(0);

	__carry->skbe(__carry);

	c1->__digits[0]->set0();
	// c1->__digits[1]->ntt();
	c1->rotate(1);


	this->__coeff2->inv_ntt(1);

	this->__coeff2->carry_propagation_1side(1, 0, 1);  // propagates first carry
	this->__coeff2->__digits[0]->set0(1);
	this->__coeff2->rotate(1);

	this->__coeff2->carry_propagation_1side(1, 0, d-1);

      }
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::multiplication_float(
									poly_hpr *c1, poly_hpr *c2)
    {
      {
	// reducing mod q

	c1->fast_modq();
	c2->fast_modq();
	c1->__digits[d-1]->float_forward();
	c2->__digits[d-1]->float_forward();
      }

      {
	// multiply by t

	c1->mulDigit_float(0, __t, 0);
	c1->mulDigit_float(1, __t, 1);
      }

      {
	// ntt the input

	c1->ntt_float();
	c2->ntt_float();
      }

      {
	// karatsuba pattern

	this->__coeff2->set0_float();

	this->__coeff2->add_float(0, c1, 0, c1, 1); // coeff2 = (c10+c11, ***)
	this->__coeff2->add_float(1, c2, 0, c2, 1); // coeff2 = (c10+c11, c20+C21)

	this->__coeff2->vtm_prod_float(0, this->__coeff2, 0, this->__coeff2,
				       1); // coeff2 = ( (c10+c11)*(c20+C21), ***)

	c1->vtm_prod_float(0, c1, 0, c2, 0); // c1 = (c10*c20, c11)
	c1->vtm_prod_float(1, c1, 1, c2, 1); // c1 = (c10*c20, c11*c21)

	this->__coeff2->sub_float(0, c1,
				  0); // coeff2 = ( (c10+c11)*(c20+C21) - c10*c20, ***)
	this->__coeff2->sub_float(0, c1,
				  1); // coeff2 = ( (c10+c11)*(c20+C21) - c10*c20 - c11*c21, ***)

	this->__coeff2->set_float(1, c1,
				  1); // coeff2 = ( (c10+c11)*(c20+C21) - c10*c20 - c11*c21, c11*c21)
	c1->set_float(1, this->__coeff2,
		      0); // c1 = (c10*c20, (c10+c11)*(c20+C21) - c10*c20 - c11*c21)

      }

      {
	// ending the division+rounding

	c1->__digits[0]->inv_ntt();
	// c1->__digits[1]->inv_ntt();
	// c1->carry_propagation(0, 1); // propagates first carry
	__carry->__cb[0] = c1->__digits[0]->__cb[0];
	__carry->__cb[1] = c1->__digits[0]->__cb[1];

	c1->__digits[0]->fbe_float();

	__carry->__cb[0] = __carry->__cb[0] - c1->__digits[0]->__cb[0];
	__carry->__cb[1] = __carry->__cb[1] - c1->__digits[0]->__cb[1];
	__carry->__cb[0] = __carry->__cb[0] * __carry->__be_tool_set->_qinv_modB[0];
	__carry->__cb[1] = __carry->__cb[1] * __carry->__be_tool_set->_qinv_modB[0];

	__carry->__cq[0].set(0);
	__carry->__cq[1].set(0);

	__carry->float_back (__carry);

	c1->__digits[0]->set0_float();
	// c1->__digits[1]->ntt();
	c1->rotate(1);


	this->__coeff2->inv_ntt_float(1);

	this->__coeff2->carry_propagation_1side_float(1, 0,
						      1);  // propagates first carry
	this->__coeff2->__digits[0]->set0_float(1);
	this->__coeff2->rotate(1);

	this->__coeff2->carry_propagation_1side_float(1, 0, d-1);

      }
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::relinearisation(poly_hpr *c1)
    {
      constexpr size_t degree = Pq::degree;
      iterator it_rlk = this->begin(); // iterator on relinearisation key
      typename poly_hpr::iterator it_c2 = this->__coeff2->begin(); // iterator on c2
      size_t shift = 1;

      for(typename poly_hpr::iterator it_c2 = this->__coeff2->begin();
	  it_c2 != this->__coeff2->end(); ++it_c2, ++shift)
	// main loop on digits of c2
	{
	  (*it_c2)->__cq[1] = (*it_c2)->__cq[1] * this->__be_tool_set->_xi_q[0];

	  for(auto it_small_pol = (*it_c2)->__cq[1].begin();
	      it_small_pol != (*it_c2)->__cq[1].end(); ++it_rlk, it_small_pol += degree)
	    //small loop on residues in base q
	    {
	      // extracting the small polynomial from c2
	      commons::tools::set_partial<Pq>(__dig_decomp->__cq[0], it_small_pol,
					      it_small_pol+degree);
	      commons::tools::set_partial<Pb>(__dig_decomp->__cb[0], it_small_pol,
					      it_small_pol+degree);
	      commons::tools::set_partial<Psk>(__dig_decomp->__csk[0], it_small_pol,
					       it_small_pol+degree);

	      __dig_decomp->ntt(0);

	      c1->mulDigit_acc(0, *it_rlk, 0, __dig_decomp, 0, shift);
	      c1->mulDigit_acc(1, *it_rlk, 1, __dig_decomp, 0, shift);
	    }
	}

      // end
      c1->inv_ntt();
      c1->__digits[0]->add(0, __carry, 0);
      c1->__digits[0]->add(1, __carry, 1);
      c1->carry_propagation(0, d-1);
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::relinearisation_float(
									 poly_hpr *c1)
    {
      constexpr size_t degree = Pq::degree;
      iterator it_rlk = this->begin(); // iterator on relinearisation key
      typename poly_hpr::iterator it_c2 = this->__coeff2->begin(); // iterator on c2
      size_t shift = 1;

      for(typename poly_hpr::iterator it_c2 = this->__coeff2->begin();
	  it_c2 != this->__coeff2->end(); ++it_c2, ++shift)
	// main loop on digits of c2
	{
	  (*it_c2)->__cq[1] = (*it_c2)->__cq[1] * this->__be_tool_set->_xi_q[0];

	  for(auto it_small_pol = (*it_c2)->__cq[1].begin();
	      it_small_pol != (*it_c2)->__cq[1].end(); ++it_rlk, it_small_pol += degree)
	    //small loop on residues in base q
	    {
	      // extracting the small polynomial from c2
	      commons::tools::set_partial<Pq>(__dig_decomp->__cq[0], it_small_pol,
					      it_small_pol+degree);
	      commons::tools::set_partial<Pb>(__dig_decomp->__cb[0], it_small_pol,
					      it_small_pol+degree);

	      __dig_decomp->ntt_float(0);

	      c1->mulDigit_acc_float(0, *it_rlk, 0, __dig_decomp, 0, shift);
	      c1->mulDigit_acc_float(1, *it_rlk, 1, __dig_decomp, 0, shift);
	    }
	}

      // end
      c1->inv_ntt_float();
      c1->__digits[0]->add_float(0, __carry, 0);
      c1->__digits[0]->add_float(1, __carry, 1);
      c1->carry_propagation_float (0, d-1);
    }



    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    size_t Context<Pq, Pb, Psk, Pmt, parameters, d>::noise_size(const poly_hpr* pol,
								const Pt* m, bool input_in_ntt) const
    {
      // gets log2(noise) for pol encrypting message m
      // pol = (Delta*m + e*(e0 - as) + e1, e*a+e2)
      //     = (Delta*m - e*a*s + e*e0 + e1, e*a + e*e2)
      constexpr size_t degree = Pq::degree;

      poly_hpr* pol_tmp = new poly_hpr();
      pol_tmp->set_mat_set_conv(this->__be_tool_set);

      pol_tmp->set0();

      base_conversions::copy_paste_conversion<Pq, Pt, typename Pt::type_t>
	(pol_tmp->__digits[d-1]->__cq[0], *m, parameters::t_t);
      base_conversions::copy_paste_conversion<Pb, Pt, typename Pt::type_t>
	(pol_tmp->__digits[d-1]->__cb[0], *m, parameters::t_t);
      base_conversions::copy_paste_conversion<Psk, Pt, typename Pt::type_t>
	(pol_tmp->__digits[d-1]->__csk[0], *m, parameters::t_t);

      pol_tmp->__digits[d-1]->ntt(0);
      pol_tmp->mulDigit(0, __Delta, 0, pol_tmp->__digits[d-1], 0);

      pol_tmp->inv_ntt(0);
      pol_tmp->sub(0, pol, 0); // pol_tmp = (e*a*s - e*e0 - e1, 0)

      pol_tmp->set(1, pol, 1); // pol_tmp = (e*a*s - e*e0 - e1, e*a + e*e2)
      pol_tmp->ntt(1);
      pol_tmp->mulDigit(1, pol_tmp, 1, this->__keys->_sk,
			1); // pol_tmp = (e*a*s - e*e0 - e1, - e*a*s - e*e2*s)
      pol_tmp->inv_ntt(1);

      pol_tmp->add(0, pol_tmp, 1, pol_tmp,
		   0); // pol_tmp = (- e*e0 - e1 - e*e2*s, - e*a*s - e*e2*s)

      pol_tmp->carry_propagation_1side(0, 0, d-1);

      size_t res = 0;
      for(size_t deg = 0; deg < degree; ++deg)
	{
	  mpz_class coeff = 0, tmp;
	  for(long int dig = d-1; dig >= 0; --dig)
	    {
	      mpz_class tmp2 = 0;
	      for (uint cmQ = 0; cmQ < __nmQ; cmQ++)
		{
		  mpz_set_ui(tmp.get_mpz_t(), pol_tmp->__digits[dig]->__cq[0](cmQ, deg));
		  tmp = (tmp * _pi_ov_p[cmQ]) % _pi[cmQ];
		  tmp2 += tmp * _qbsk / _pi[cmQ];
		}

	      for (uint cmB = 0; cmB < __nmB; cmB++)
		{
		  mpz_set_ui(tmp.get_mpz_t(), pol_tmp->__digits[dig]->__cb[0](cmB, deg));
		  tmp = (tmp * _pi_ov_p[__nmQ+cmB]) % _pi[__nmQ+cmB];
		  tmp2 += tmp * _qbsk / _pi[__nmQ+cmB];
		}

	      mpz_set_ui(tmp.get_mpz_t(), pol_tmp->__digits[dig]->__csk[0](0, deg));
	      tmp = (tmp * _pi_ov_p[__nmQ+__nmB]) % _pi[__nmQ+__nmB];
	      tmp2 += tmp * _qbsk / _pi[__nmQ+__nmB];

	      coeff = coeff * _q + (tmp2 % _qbsk);
	    }

	  coeff = coeff % _p;

	  if(2*coeff >= _p)
	    coeff = coeff - _p;

	  coeff = abs(coeff);

	  size_t log_coeff = mpz_sizeinbase(coeff.get_mpz_t(), 2);
	  if(log_coeff > res)
	    res = log_coeff;
	}

      delete pol_tmp;
      return res;
    }



    //////////////// print fcts. ////////////////
    /////////////////////////////////////////////

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::print_baseq()
    {
      std::cout << std::endl << std::endl;
      std::cout << "printing base q";
      std::cout << std::endl << std::endl;

      for(uint i = 0; i < __nmQ; i++)
	std::cout << Pq::base.params.P[__imQ+i] << " ";

      std::cout << std::endl << std::endl;
      std::cout << "end base q";
      std::cout << std::endl << std::endl;
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::print_baseb()
    {
      std::cout << std::endl << std::endl;
      std::cout << "printing base b";
      std::cout << std::endl << std::endl;

      for(uint i = 0; i < __nmB; i++)
	std::cout << Pb::base.params.P[__imB+i] << " ";

      std::cout << std::endl << std::endl;
      std::cout << "end base b";
      std::cout << std::endl << std::endl;
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::print_basesk()
    {
      std::cout << std::endl << std::endl;
      std::cout << "printing base sk";
      std::cout << std::endl << std::endl;

      std::cout << Psk::base.params.P[__imSK];

      std::cout << std::endl << std::endl;
      std::cout << "end base sk";
      std::cout << std::endl << std::endl;
    }

    template<class Pq, class Pb, class Psk, class Pmt, class parameters, size_t d>
    void Context<Pq, Pb, Psk, Pmt, parameters, d>::print_basegt()
    {
      std::cout << std::endl << std::endl;
      std::cout << "printing base gamma*t";
      std::cout << std::endl << std::endl;

      std::cout << parameters::gamma << " * " << parameters::t << " = " << (((
									      uint)parameters::gamma)*parameters::t);

      std::cout << std::endl << std::endl;
      std::cout << "end base gamma*t";
      std::cout << std::endl << std::endl;
    }
  }
}


#endif /* hpr_h */
