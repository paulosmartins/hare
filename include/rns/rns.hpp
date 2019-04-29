#ifndef rns_h
#define rns_h

#include "rns/common_rns.hpp"
#include <gmp.h>
#include <mpfr.h>

namespace rns
{
// contains ct=(c0,c1) in bases q, b and m_sk (and m_tilde if RED)
  template<class Pq, class Pb, class Psk, class Pmt>
  struct Ciphertext
  {
    using main_type = typename Pq::value_type;
    using main_greater_type = typename Pq::greater_value_type;
    using main_signed_type = typename Pq::signed_value_type;
    using main_double = typename Pq::main_double;
    using main_greater_double = typename Pq::main_greater_double;


    static constexpr uint __nmQ = Pq::nmoduli;
    static constexpr uint __nmB = Pb::nmoduli;
    static constexpr uint __imQ = Pq::indmod0;
    static constexpr uint __imB = Pb::indmod0;
    static constexpr uint __imSK = Psk::indmod0;

    Pq* __cq;
    Pb* __cb;
    Psk* __csk;
    // typedef commons::fast_pol::PolMT<main_signed_type> Pmt;
    Pmt* __ct;
    bool __red;

    Ciphertext();
    ~Ciphertext();
    Ciphertext(const Ciphertext&);
    Ciphertext& operator=(const Ciphertext&);

    void set(const Ciphertext&);
    void set0();

    bool operator==(const Ciphertext &other) const;

    void add(const Ciphertext &c);
    void add(const Ciphertext &c1, const Ciphertext &c2);
    void addE(const Ciphertext &c1, const Ciphertext &c2);
    void sub(const Ciphertext &c);
    void subE(const Ciphertext &c1, const Ciphertext &c2);

    inline void inv_ntt0();
    inline void inv_ntt1();
    inline void ntt();
    inline void inv_ntt0_float();
    inline void inv_ntt1_float();
    inline void ntt_float();

    void mul(Ciphertext &c2, const Ciphertext &c, const Pq &aq, const Pb &aB,
             const Psk &ask);
    void mul_float(Ciphertext &c2, const Ciphertext &c, const Pq &aq);

    void SmMRq(const std::array<main_type, __nmB> &vB,
               const std::array<main_type, 1> &vsk);
    void fast_base_conversion_red(const
                                  std::array<std::array<main_type, __nmQ>, __nmB> &Mb,
                                  const std::array<main_type, __nmQ> &Msk,
                                  const std::array<typename Pmt::type_mtilde_, __nmQ> &Mt);
    void fast_base_conversion(const std::array<std::array<main_type, __nmQ>, __nmB>
                              &Mb, const std::array<main_type, __nmQ> &Msk);
    void base_conversion_float(const std::array<main_double, __nmQ> &q_i_inv,
                               const std::array<std::array<main_type, __nmQ>, __nmB> &q_i,
                               const std::array<main_type, __nmB> &q);
    void base_conversion_back_float(const std::array<main_double, __nmB> &b_i_inv,
                                    const std::array<std::array<main_type, __nmB>, __nmQ> &b_i,
                                    const std::array<main_type, __nmQ> &b);
  };

  template<class Pq, class Pb, class Psk, class Pmt>
  Ciphertext<Pq, Pb, Psk, Pmt>::Ciphertext()
  {
    __cq = commons::mem_manag::alloc_aligned<Pq, 32>(2);
    __cb = commons::mem_manag::alloc_aligned<Pb, 32>(2);
    __csk = commons::mem_manag::alloc_aligned<Psk, 32>(2);
    __ct = commons::mem_manag::alloc_aligned<Pmt, 32>(2);
    __red = false;
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  Ciphertext<Pq, Pb, Psk, Pmt>::Ciphertext(const Ciphertext<Pq, Pb, Psk, Pmt> &other)
  {
    __cq = commons::mem_manag::alloc_aligned<Pq, 32>(2);
    __cb = commons::mem_manag::alloc_aligned<Pb, 32>(2);
    __csk = commons::mem_manag::alloc_aligned<Psk, 32>(2);
    __ct = commons::mem_manag::alloc_aligned<Pmt, 32>(2);

    __cq[0] = other.__cq[0];
    __cq[1] = other.__cq[1];
    __cb[0] = other.__cb[1];
    __csk[0] = other.__csk[0];
    __csk[1] = other.__csk[1];
    __ct[0] = other.__ct[0];
    __ct[1] = other.__ct[1];
    __red = other.__red;
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  Ciphertext<Pq, Pb, Psk, Pmt>&
  Ciphertext<Pq, Pb, Psk, Pmt>::operator=
  (const Ciphertext<Pq, Pb, Psk, Pmt> &other)
  {
    __cq[0] = other.__cq[0];
    __cq[1] = other.__cq[1];
    __cb[0] = other.__cb[1];
    __csk[0] = other.__csk[0];
    __csk[1] = other.__csk[1];
    __ct[0] = other.__ct[0];
    __ct[1] = other.__ct[1];
    __red = other.__red;
    return *this;
  }
  
  
  template<class Pq, class Pb, class Psk, class Pmt>
  Ciphertext<Pq, Pb, Psk, Pmt>::~Ciphertext()
  {
    commons::mem_manag::free_aligned<Pq>(2, __cq);
    commons::mem_manag::free_aligned<Pb>(2, __cb);
    commons::mem_manag::free_aligned<Psk>(2, __csk);
    commons::mem_manag::free_aligned<Pmt>(2, __ct);
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::set(const Ciphertext &other)
  {
    __cq[0] = other.__cq[0];
    __cq[1] = other.__cq[1];
    __cb[0] = other.__cb[0];
    __cb[1] = other.__cb[1];
    __csk[0] = other.__csk[0];
    __csk[1] = other.__csk[1];
    __ct[0] = other.__ct[0];
    __ct[1] = other.__ct[1];
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::set0()
  {
    __cq[0] = 0;
    __cq[1] = 0;
    __cb[0] = 0;
    __cb[1] = 0;
    __csk[0] = 0;
    __csk[1] = 0;
    __ct[0].set0();
    __ct[1].set0();
    __red = false;
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  bool Ciphertext<Pq, Pb, Psk, Pmt>::operator==(const Ciphertext &other) const
  {
    return (__cq[0] == other.__cq[0]) && (__cb[0] == other.__cb[0])
           && (__ct[0] == other.__ct[0]) && (__cq[1] == other.__cq[1])
           && (__cb[1] == other.__cb[1]) && (__ct[1] == other.__ct[1])
           && (__csk[0] == other.__csk[0]) && (__csk[1] == other.__csk[1]);
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::add(const Ciphertext &c)
  {
    // only add the Pq components
    __cq[0] = __cq[0] + c.__cq[0];
    __cq[1] = __cq[1] + c.__cq[1];
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::add(const Ciphertext &c1,
                                         const Ciphertext &c2)
  {
    // only add the Pq components
    __cq[0] = __cq[0] + c1.__cq[0] + c2.__cq[0];
    __cq[1] = __cq[1] + c1.__cq[1] + c2.__cq[1];
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::addE(const Ciphertext &c1,
                                          const Ciphertext &c2)
  {
    // only add the Pq components
    __cq[0] = c1.__cq[0] + c2.__cq[0];
    __cq[1] = c1.__cq[1] + c2.__cq[1];
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::sub(const Ciphertext &c)
  {
    // only sub the Pq components
    __cq[0] = __cq[0] - c.__cq[0];
    __cq[1] = __cq[1] - c.__cq[1];
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::subE(const Ciphertext &c1,
                                          const Ciphertext &c2)
  {
    // only sub the Pq components
    __cq[0] = c1.__cq[0] - c2.__cq[0];
    __cq[1] = c1.__cq[1] - c2.__cq[1];
  }


  template<class Pq, class Pb, class Psk, class Pmt>
  inline void Ciphertext<Pq, Pb, Psk, Pmt>::inv_ntt0()
  {
    __cq[0].invntt_pow_invphi();
    __cb[0].invntt_pow_invphi();
    __csk[0].invntt_pow_invphi();
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  inline void Ciphertext<Pq, Pb, Psk, Pmt>::inv_ntt1()
  {
    __cq[1].invntt_pow_invphi();
    __cb[1].invntt_pow_invphi();
    __csk[1].invntt_pow_invphi();
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  inline void Ciphertext<Pq, Pb, Psk, Pmt>::ntt()
  {
    __cq[0].ntt_pow_phi();
    __cb[0].ntt_pow_phi();
    __csk[0].ntt_pow_phi();
    __cq[1].ntt_pow_phi();
    __cb[1].ntt_pow_phi();
    __csk[1].ntt_pow_phi();
  }
  template<class Pq, class Pb, class Psk, class Pmt>
  inline void Ciphertext<Pq, Pb, Psk, Pmt>::inv_ntt0_float()
  {
    __cq[0].invntt_pow_invphi();
    __cb[0].invntt_pow_invphi();
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  inline void Ciphertext<Pq, Pb, Psk, Pmt>::inv_ntt1_float()
  {
    __cq[1].invntt_pow_invphi();
    __cb[1].invntt_pow_invphi();
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  inline void Ciphertext<Pq, Pb, Psk, Pmt>::ntt_float()
  {
    __cq[0].ntt_pow_phi();
    __cb[0].ntt_pow_phi();
    __cq[1].ntt_pow_phi();
    __cb[1].ntt_pow_phi();
  }


// multiply 2 ciphertexts
// and multiply by some precomputed data
  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::mul(Ciphertext &c2, const Ciphertext &c,
                                         const Pq &aq, const Pb &aB, const Psk &ask)
  {
    __cq[0] = __cq[0]*aq;
    __cq[1] = __cq[1]*aq;

    __cb[0] = __cb[0]*aB;
    __cb[1] = __cb[1]*aB;

    __csk[0] = __csk[0]*ask;
    __csk[1] = __csk[1]*ask;

    c2.__cq[0] = __cq[0] - __cq[1]; // c0-c1
    c2.__cb[0] = __cb[0] - __cb[1];
    c2.__csk[0] = __csk[0] - __csk[1];

    c2.__cq[1] = c.__cq[1] - c.__cq[0]; // c'1-c'0
    c2.__cb[1] = c.__cb[1] - c.__cb[0];
    c2.__csk[1] = c.__csk[1] - c.__csk[0];

    c2.__cq[1] = c2.__cq[0]*c2.__cq[1]; // (c0-c1)(c'1-c'0)
    c2.__cb[1] = c2.__cb[0]*c2.__cb[1];
    c2.__csk[1] = c2.__csk[0]*c2.__csk[1];

    c2.__cq[0] = __cq[1]*c.__cq[1]; // c1c'1
    c2.__cb[0] = __cb[1]*c.__cb[1];
    c2.__csk[0] = __csk[1]*c.__csk[1];

    __cq[0] = __cq[0]*c.__cq[0]; // c0c'0
    __cb[0] = __cb[0]*c.__cb[0];
    __csk[0] = __csk[0]*c.__csk[0];

    __cq[1] = c2.__cq[0] + c2.__cq[1] + __cq[0];
    __cb[1] = c2.__cb[0] + c2.__cb[1] + __cb[0];
    __csk[1] = c2.__csk[0] + c2.__csk[1] + __csk[0];
    __red = false;
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::mul_float(Ciphertext &c2,
					       const Ciphertext &c,
					       const Pq &aq)
  {
    __cq[0] = __cq[0] * aq;
    __cq[1] = __cq[1] * aq;
    
    c2.__cq[0] = __cq[0] - __cq[1]; // c0-c1
    c2.__cb[0] = __cb[0] - __cb[1];

    c2.__cq[1] = c.__cq[1] - c.__cq[0]; // c'1-c'0
    c2.__cb[1] = c.__cb[1] - c.__cb[0];

    c2.__cq[1] = c2.__cq[0]*c2.__cq[1]; // (c0-c1)(c'1-c'0)
    c2.__cb[1] = c2.__cb[0]*c2.__cb[1];

    c2.__cq[0] = __cq[1]*c.__cq[1]; // c1c'1
    c2.__cb[0] = __cb[1]*c.__cb[1];

    __cq[0] = __cq[0]*c.__cq[0]; // c0c'0
    __cb[0] = __cb[0]*c.__cb[0];

    __cq[1] = c2.__cq[0] + c2.__cq[1] + __cq[0];
    __cb[1] = c2.__cb[0] + c2.__cb[1] + __cb[0];
  }


// Montgomery reduction modulo q (using m_tilde)
// Algo. 2 in paper
  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::SmMRq(const std::array<main_type, __nmB> &vB,
      const std::array<main_type, 1> &vsk)
  {
    if(not __red)
      {
        std::array<main_signed_type, Pq::degree> acc;
        __ct[0].get_modc(acc);
        commons::reductions::SmMRq<Pb>(__cb[0], acc, vB);
        commons::reductions::SmMRq<Psk>(__csk[0], acc, vsk);
        __ct[1].get_modc(acc);
        commons::reductions::SmMRq<Pb>(__cb[1], acc, vB);
        commons::reductions::SmMRq<Psk>(__csk[1], acc, vsk);
        __red = true;
      }
  }

// base conversion : extending c mod q toward base b, m_sk and m_tilde
// eq. (2) in paper
  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::fast_base_conversion_red(
    const std::array<std::array<main_type, __nmQ>,__nmB> &Mb,
    const std::array<main_type, __nmQ> &Msk
    ,
    const std::array<typename Pmt::type_mtilde_, __nmQ> &Mt
  )
  {
    // extends cq to other bases
    __cb[0] = 0;
    __csk[0] = 0;
    commons::base_conversions::fast_1m_to_1m_11<Pb, Psk, Pq>(__cb[0], __csk[0],
        __cq[0], Mb, Msk);
    __ct[0].set0();
    __ct[0].template fast_base_conversion<Pq>(__cq[0], Mt);

    __cb[1] = 0;
    __csk[1] = 0;
    commons::base_conversions::fast_1m_to_1m_11<Pb, Psk, Pq>(__cb[1], __csk[1],
        __cq[1], Mb, Msk);
    __ct[1].set0();
    __ct[1].template fast_base_conversion<Pq>(__cq[1], Mt);
    __red = false;
  }

// base conversion : extending c mod q toward base b, m_sk and m_tilde
// eq. (2) in paper
  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::fast_base_conversion(
    const std::array<std::array<main_type, __nmQ>,__nmB> &Mb,
    const std::array<main_type, __nmQ> &Msk
  )
  {
    // extends cq to other bases
    __cb[0] = 0;
    __csk[0] = 0;
    commons::base_conversions::fast_1m_to_1m_11<Pb, Psk, Pq>(__cb[0], __csk[0],
        __cq[0], Mb, Msk);

    __cb[1] = 0;
    __csk[1] = 0;
    commons::base_conversions::fast_1m_to_1m_11<Pb, Psk, Pq>(__cb[1], __csk[1],
        __cq[1], Mb, Msk);
  }


  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::base_conversion_float(
    const std::array<typename Pq::main_double, __nmQ> &q_i_inv,
    const std::array<std::array<main_type, __nmQ>, __nmB> &q_i,
    const std::array<main_type, __nmB> &q)
  {
    commons::base_conversions::fast_1m_to_1m_float<Pb, Pq>(__cb[0], __cq[0], q_i_inv,
        q_i, q);
    commons::base_conversions::fast_1m_to_1m_float<Pb, Pq>(__cb[1], __cq[1], q_i_inv,
        q_i, q);
  }

  template<class Pq, class Pb, class Psk, class Pmt>
  void Ciphertext<Pq, Pb, Psk, Pmt>::base_conversion_back_float(
    const std::array<typename Pq::main_double, __nmB> &b_i_inv,
    const std::array<std::array<main_type, __nmB>, __nmQ> &b_i,
    const std::array<main_type, __nmQ> &b)
  {
    commons::base_conversions::fast_1m_to_1m_float<Pq, Pb>(__cq[0], __cb[0], b_i_inv,
        b_i, b);
    commons::base_conversions::fast_1m_to_1m_float<Pq, Pb>(__cq[1], __cb[1], b_i_inv,
        b_i, b);
  }


//========== General context - structure gathering precomputations ==========//

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  struct Context
  {

    typedef commons::fast_pol::PolT<parameters, Pq::degree> PT; // message space
    using main_type = typename Pq::value_type;
    using main_greater_type = typename Pq::greater_value_type;
    using main_signed_type = typename Pq::signed_value_type;
    using main_double = typename Pq::main_double;
    using main_greater_double = typename Pq::main_greater_double;

    typedef Ciphertext<Pq, Pb, Psk, Pmt> _Ciphertext;
    bool __red;

    static constexpr uint __nmQ = Pq::nmoduli;  // number of moduli in base q
    static constexpr uint __nmB = Pb::nmoduli;  // number of moduli in base b
    static constexpr uint __imQ =
      Pq::indmod0;  // index of first modulus of base q in nfl::params<main_type>::P
    static constexpr uint __imB =
      Pb::indmod0;  // index of first modulus of base b in nfl::params<main_type>::P
    static constexpr uint __imSK =
      Psk::indmod0; // index of modulus m_sk in nfl::params<main_type>::P
    typedef commons::fast_pol::PolT<parameters, Pq::degree> Pt;
    // useful for computing noise size
    // not critical: only for test purpose (can use gmp)
    mpz_class __q, __b;
    Pq* __pol_noise; // 0: xi_q(Delta); 1: 1/mtilde in NTT; 2: s/mtilde in NTT
    std::vector<mpz_class> __q_ov_qi; // (q/qi)_i=1..k
    std::vector<mpz_class> __b_ov_bi; // (q/qi)_i=1..k

    // precomputations for encryption
    Pq* __pol_enc; // 0: a_q; 1: b_q; 2: Delta_q; 3: mtilde*qi/q
    Pq* __pol_enc_float; // 0: a_q; 1: b_q; 2: Delta_q;
    std::array<std::array<main_type, __nmQ>, __nmB> __M_q_enc_b;
    std::array<main_type, __nmQ> __M_q_enc_sk;
    std::array<typename Pmt::type_mtilde_, __nmQ> __M_q_enc_tilde;

    // decryption
    Pq* __pol_dec;
    Pq* __pol_dec_float;
    Pq* __sk;
    std::array<typename parameters::type_gamma_t, __nmQ> __M_q_dec_gamma_t;
    main_greater_double __t_ov_qi[__nmQ];

    // SMRmq
    std::array<main_type, __nmB> __delta_b_tilde_mul;
    std::array<main_type, 1> __delta_sk_tilde_mul;

    // Multiplication
    Pq* __evka, *__evkb;
    Pq* __evka_float, *__evkb_float;
    Pq* __decomp;
    _Ciphertext* __coeff2;
    // b
    Pb* __polb_mul;
    std::array<std::array<main_type, __nmQ>, __nmB> __M_q_mul_b_1;
    // sk
    Psk* __polsk_mul;
    std::array<main_type, __nmQ> __M_q_mul_sk_1;
    std::array<main_type, __nmB> __M_b_mul_sk;
    // q
    Pq* __polq_mul;
    Pq* __polq_mul_1;
    std::array<main_type,__nmQ> __delta_q_sk_mul;
    std::array<std::array<main_type, __nmB>, __nmQ> __M_b_mul_q;
    std::array<main_type, __nmQ> __delta_q_sk_mul2;
    std::array<std::array<main_type, __nmB>, __nmQ> __M_b_mul_q2;

    //HPS-based extensions
    std::array<main_double, __nmQ> __M_q_i_inv;
    std::array<std::array<main_type, __nmQ>, __nmB> __M_q_i;
    std::array<main_type, __nmB> __M_q;
    std::array<std::array<main_type, __nmQ>, __nmB> __M_q_i_1;
    std::array<main_type, __nmB> __M_q_1;

    std::array<main_double, __nmQ> thetas;
    std::array<main_type, __nmB> lambdas;
    std::array<std::array<main_type, __nmQ>, __nmB> omegas;

    std::array<main_double, __nmB> __M_b_i_inv;
    std::array<std::array<main_type, __nmB>, __nmQ> __M_b_i;
    std::array<std::array<main_type, __nmB>, __nmQ> __M_b_i_1;
    std::array<main_type, __nmQ> __M_b;
    std::array<main_type, __nmQ> __M_b_1;

    Pq* __qi_ov_q;
    Pb* __bi_ov_b;

    std::array<main_greater_double, __nmQ> dec_thetas;

    // constructor
    Context(bool red = true);
    Context(const Context &other);
    Context& operator=(const Context &other);
    void allocate();
    void copy_from(const Context &other);

    // destructor
    ~Context();

    void set_precomp(const Pq* keys);

    size_t size_noise(const _Ciphertext &c, const Pt &m, bool input_in_ntt) const;
    inline void extend(_Ciphertext &c);
    inline void extend_float(_Ciphertext &c);
    inline void extend_float_1(_Ciphertext &c);
    inline void extend_back_float(_Ciphertext &c);
    void encrypt(_Ciphertext &c, const Pq &m, bool output_in_ntt = false);
    void encrypt_float(_Ciphertext &c, const Pq &m, bool output_in_ntt = false);
    void decrypt(Pt &mdec, _Ciphertext &c);
    void decrypt_not_ntt(Pt &mdec, _Ciphertext &c);
    void decrypt_float_old(Pt &mdec, _Ciphertext &c);
    void decrypt_float(Pt &mdec, _Ciphertext &c);
    void decrypt_float_not_ntt(Pt &mdec, _Ciphertext &c);
    void mult(_Ciphertext &c1, _Ciphertext &c2, bool output_in_ntt = false);
    void mult_float_step1(_Ciphertext &c1, _Ciphertext &c2,
                          bool output_in_ntt = false);
    void mult_float_step2(_Ciphertext &c1, _Ciphertext &c2,
                          bool output_in_ntt = false);
    void mult_float_step3(_Ciphertext &c1, _Ciphertext &c2,
                          bool output_in_ntt = false);
    void mult_float_step4(_Ciphertext &c1, _Ciphertext &c2,
                          bool output_in_ntt = false);
    void mult_float(_Ciphertext &c1, _Ciphertext &c2, bool output_in_ntt = false);
    void relinearisation(_Ciphertext &c1, bool output_in_ntt = false);
    void relinearisation_float(_Ciphertext &c1, bool output_in_ntt = false);
  };

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::allocate()
  {
    __pol_noise = commons::mem_manag::alloc_aligned<Pq, 32>(3);
    __pol_enc = commons::mem_manag::alloc_aligned<Pq, 32>(4);
    __pol_enc_float = commons::mem_manag::alloc_aligned<Pq, 32>(4);
    __polb_mul = commons::mem_manag::alloc_aligned<Pb, 32>(1);
    __sk = commons::mem_manag::alloc_aligned<Pq, 32>(1);
    __pol_dec = commons::mem_manag::alloc_aligned<Pq, 32>(2);
    __pol_dec_float = commons::mem_manag::alloc_aligned<Pq, 32>(1);
    __polq_mul = commons::mem_manag::alloc_aligned<Pq, 32>(1);
    __polq_mul_1 = commons::mem_manag::alloc_aligned<Pq, 32>(1);
    __polsk_mul = commons::mem_manag::alloc_aligned<Psk, 32>(1);
    __evka = commons::mem_manag::alloc_aligned<Pq, 32>(__nmQ);
    __evkb = commons::mem_manag::alloc_aligned<Pq, 32>(__nmQ);
    __evka_float = commons::mem_manag::alloc_aligned<Pq, 32>(__nmQ);
    __evkb_float = commons::mem_manag::alloc_aligned<Pq, 32>(__nmQ);
    __decomp = commons::mem_manag::alloc_aligned<Pq, 32>(3);
    __coeff2 = commons::mem_manag::alloc_aligned<_Ciphertext, 32>(1);
    __qi_ov_q = commons::mem_manag::alloc_aligned<Pq, 32>(1);
    __bi_ov_b = commons::mem_manag::alloc_aligned<Pb, 32>(1);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::copy_from
  (const Context<Pq, Pb, Psk, Pmt, parameters>& other)
  {
    __red = other.__red;
    __q = other.__q;
    __b = other.__b;
    __pol_noise[0] = other.__pol_noise[0];
    __pol_noise[1] = other.__pol_noise[1];
    __pol_noise[2] = other.__pol_noise[2];
    __q_ov_qi = other.__q_ov_qi;
    __b_ov_bi = other.__b_ov_bi;
    __pol_enc[0] = other.__pol_enc[0];
    __pol_enc[1] = other.__pol_enc[1];
    __pol_enc[2] = other.__pol_enc[2];
    __pol_enc[3] = other.__pol_enc[3];
    __pol_enc_float[0] = other.__pol_enc_float[0];
    __pol_enc_float[1] = other.__pol_enc_float[1];
    __pol_enc_float[2] = other.__pol_enc_float[2];
    __pol_enc_float[3] = other.__pol_enc_float[3];
    __M_q_enc_b = other.__M_q_enc_b;
    __M_q_enc_sk = other.__M_q_enc_sk;
    __M_q_enc_tilde = other.__M_q_enc_tilde;
    __pol_dec[0] = other.__pol_dec[0];
    __pol_dec[1] = other.__pol_dec[1];
    __pol_dec_float[0] = other.__pol_dec_float[0];
    __sk[0] = other.__sk[0];
    __M_q_dec_gamma_t = other.__M_q_dec_gamma_t;
    std::copy(&other.__t_ov_qi[0], &other.__t_ov_qi[0]+__nmQ, &__t_ov_qi[0]);
    __delta_b_tilde_mul = other.__delta_b_tilde_mul;
    __delta_sk_tilde_mul = other.__delta_sk_tilde_mul;
    for (size_t i = 0; i < __nmQ; i++)
      {
	__evka[i] = other.__evka[i];
	__evkb[i] = other.__evkb[i];
	__evka_float[i] = other.__evka_float[i];
	__evkb_float[i] = other.__evkb_float[i];
      }
    __decomp[0] = other.__decomp[0];
    __decomp[1] = other.__decomp[1];
    __decomp[2] = other.__decomp[2];
    __coeff2[0] = other.__coeff2[0];
    __polb_mul[0] = other.__polb_mul[0];
    __M_q_mul_b_1 = other.__M_q_mul_b_1;
    __polsk_mul[0] = other.__polsk_mul[0];
    __M_q_mul_sk_1 = other.__M_q_mul_sk_1;
    __M_b_mul_sk = other.__M_b_mul_sk;
    __polq_mul[0] = other.__polq_mul[0];
    __polq_mul_1[0] = other.__polq_mul_1[0];
    __delta_q_sk_mul = other.__delta_q_sk_mul;
    __M_b_mul_q = other.__M_b_mul_q;
    __delta_q_sk_mul2 = other.__delta_q_sk_mul2;
    __M_b_mul_q2 = other.__M_b_mul_q2;
    __M_q_i_inv = other.__M_q_i_inv;
    __M_q_i = other.__M_q_i;
    __M_q = other.__M_q;
    __M_q_i_1 = other.__M_q_i_1;
    __M_q_1 = other.__M_q_1;
    thetas = other.thetas;
    lambdas = other.lambdas;
    omegas = other.omegas;
    __M_b_i_inv = other.__M_b_i_inv;
    __M_b_i = other.__M_b_i;
    __M_b_i_1 = other.__M_b_i_1;
    __M_b = other.__M_b;
    __M_b_1 = other.__M_b_1;
    __qi_ov_q[0] = other.__qi_ov_q[0];
    __bi_ov_b[0] = other.__bi_ov_b[0];
    dec_thetas = other.dec_thetas;
  }
  
  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  Context<Pq, Pb, Psk, Pmt, parameters>::Context(bool red)
    : __red(red)
  {
    allocate();
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  Context<Pq, Pb, Psk, Pmt, parameters>::Context
  (const Context<Pq, Pb, Psk, Pmt, parameters>& other)
  {
    allocate();
    copy_from(other);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  Context<Pq, Pb, Psk, Pmt, parameters>&
  Context<Pq, Pb, Psk, Pmt, parameters>::operator=
  (const Context<Pq, Pb, Psk, Pmt, parameters>& other)
  {
    copy_from(other);
    return *this;
  }
  
  
  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  Context<Pq, Pb, Psk, Pmt, parameters>::~Context()
  {
    commons::mem_manag::free_aligned<Pq>(1, __sk);
    commons::mem_manag::free_aligned<Pq>(2, __pol_dec);
    commons::mem_manag::free_aligned<Pq>(1, __pol_dec_float);
    commons::mem_manag::free_aligned<Pq>(1, __polq_mul);
    commons::mem_manag::free_aligned<Pq>(1, __polq_mul_1);
    commons::mem_manag::free_aligned<Psk>(1, __polsk_mul);
    commons::mem_manag::free_aligned<Pq>(__nmQ, __evka);
    commons::mem_manag::free_aligned<Pq>(__nmQ, __evkb);
    commons::mem_manag::free_aligned<Pq>(__nmQ, __evka_float);
    commons::mem_manag::free_aligned<Pq>(__nmQ, __evkb_float);
    commons::mem_manag::free_aligned<Pq>(3, __decomp);
    commons::mem_manag::free_aligned<Pq>(3, __pol_noise);
    commons::mem_manag::free_aligned<Pq>(4, __pol_enc);
    commons::mem_manag::free_aligned<Pq>(4, __pol_enc_float);
    commons::mem_manag::free_aligned<Pb>(1, __polb_mul);
    commons::mem_manag::free_aligned<_Ciphertext>(1, __coeff2);
    commons::mem_manag::free_aligned<Pq>(1, __qi_ov_q);
    commons::mem_manag::free_aligned<Pb>(1, __bi_ov_b);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::set_precomp(const Pq* keys)
  {
    constexpr size_t degree = Pq::degree;
    const Pq *a = &keys[0];
    const Pq *b = &keys[1];
    const Pq *s = &keys[2];

    mpz_class tmpg, tmpzmt = 1, tmpt, modsk;
    mpz_set_ui(modsk.get_mpz_t(), Psk::base.params.P[__imSK]); // msk
    mpz_set_ui(tmpzmt.get_mpz_t(), Pmt::mtilde);
    mpz_set_ui(tmpg.get_mpz_t(), parameters::gamma);
    mpz_set_ui(tmpt.get_mpz_t(), parameters::t);

    // some useful quantities
    std::array<main_type, __nmQ> invQi;
    {
      mpz_class tmp, mod;
      __q = 1;
      for (uint cmQ = 0; cmQ < __nmQ; cmQ++)
        {
          mpz_set_ui(mod.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]);
          __q *= mod;
        }
      __b = 1;
      for (uint cmB = 0; cmB < __nmB; cmB++)
        {
          mpz_set_ui(mod.get_mpz_t(), Pb::base.params.P[__imB+cmB]);
          __b *= mod;
        }

      for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
        {
          mpz_set_ui(mod.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]);
          tmp = __q/mod;
          __q_ov_qi.push_back(tmp);
          mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(), mod.get_mpz_t());
          invQi[cmQ] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
        }

      for (uint cmB = 0; cmB < __nmB; cmB++)
        {
          mpz_set_ui (mod.get_mpz_t (),
                      Pb::base.params.P[__imB+cmB]);
          tmp = __b / mod;
          __b_ov_bi.push_back (tmp);
        }

    }

    mpz_class Delta = __q>>parameters::_bs_t;

    // for noise analysis
    {
      mpz_class mod, tmp;

      __pol_noise[0] = 0;
      __pol_noise[1] = 0;

      auto* it0 = __pol_noise[0].begin(), *it1 = __pol_noise[1].begin();

      for(uint cmQ = 0; cmQ < __nmQ; cmQ++, it0 += degree, it1 += degree)
        {
          mpz_set_ui(mod.get_mpz_t(), Pq::base.params.P[cmQ+__imQ]);
          mpz_set_ui(tmp.get_mpz_t(), invQi[cmQ]);
          tmp = (Delta*tmp) % mod;
          std::fill(it0, it0+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
          mpz_invert(tmp.get_mpz_t(), tmpzmt.get_mpz_t(), mod.get_mpz_t());
          std::fill(it1, it1+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
        }
      __pol_noise[2] = __pol_noise[1]**s;
    }

    //Encryption
    {
      __pol_enc[2] = 0;
      __pol_enc[3] = 0;

      mpz_class modq, tmp;

      auto* it2 = __pol_enc[2].begin(), *it3 = __pol_enc[3].begin();

      for(uint cmQ = 0; cmQ < __nmQ; cmQ++, it2 += degree, it3 += degree)
        {
          // q
          mpz_set_ui(tmp.get_mpz_t(), invQi[cmQ]);
          mpz_set_ui(modq.get_mpz_t(), Pq::base.params.P[cmQ+__imQ]);
          tmp = (tmp*tmpzmt) % modq; // mtilde*qi/q mod qi

          std::fill(it3, it3+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

          tmp = __q;
          if(__nmQ > 1)
            tmp = __q/modq;

          // b
          mpz_class modb, tmp2;
          for(uint cmB = 0; cmB < __nmB; cmB++)
            {
              mpz_set_ui(modb.get_mpz_t(), Pb::base.params.P[cmB+__imB]);
              mpz_invert(tmp2.get_mpz_t(), tmpzmt.get_mpz_t(), modb.get_mpz_t());
              tmp2 = (tmp2*tmp) % modb;
              __M_q_enc_b[cmB][cmQ] = static_cast<main_type>(mpz_get_ui(
                                        tmp2.get_mpz_t())); // (q/mtilde/q_cmQ) mod m_cmB
            }

          // msk
          mpz_invert(tmp2.get_mpz_t(), tmpzmt.get_mpz_t(),
                     modsk.get_mpz_t()); // 1/mtilde mod msk
          tmp2 *= __q_ov_qi[cmQ]; // q/q_cmQ*((1/mtilde) mod msk)
          tmp2 = tmp2 % modsk; // (q/q_cmQ/mtilde) mod msk
          __M_q_enc_sk[cmQ] = static_cast<main_type>(mpz_get_ui(
                                tmp2.get_mpz_t())); // (q/mtilde/q_cmQ) mod m_SK
          // mtilde
          tmp *= -modq; // -q_cmQ
          mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(),
                     tmpzmt.get_mpz_t()); // (-1/q_cmQ) mod mtilde
          __M_q_enc_tilde[cmQ] = static_cast<typename Pmt::type_mtilde_>(mpz_get_ui(
                                   tmp.get_mpz_t())); // (-1/q_cmQ) mod m_tilde
        }
      commons::tools::set_1coeff_mpz<Pq>(__pol_enc[2], Delta, 0);
      __pol_enc[2].ntt_pow_phi();
      __pol_enc[2] = __pol_enc[2]*__pol_enc[3]; // Delta*mtilde*qi/q
      __pol_enc[0] = *a*__pol_enc[3];      // a*mtilde*qi/q
      __pol_enc[1] = *b*__pol_enc[3];      // b*mtilde*qi/q
    }

    // Decryption
    {
      mpz_class tmp = tmpg*tmpt, tmp2;
      for(uint cm = 0; cm < __nmQ; cm++)
        {
          // gamma
          mpz_set_ui(tmp2.get_mpz_t(), Pq::base.params.P[cm+__imQ]);
          tmp2 = -tmp2; // -qi
          mpz_invert(tmp2.get_mpz_t(), tmp2.get_mpz_t(),
                     tmp.get_mpz_t()); // (-1/qi) mod gamma*t
          __M_q_dec_gamma_t[cm] = static_cast<typename parameters::type_gamma_t>
                                  (mpz_get_ui(tmp2.get_mpz_t()));
        }
      __pol_dec[0] = 0;
      __pol_dec[1] = 0;
      mpz_invert(tmp2.get_mpz_t(), tmpzmt.get_mpz_t(), __q.get_mpz_t());
      tmp2 = (tmp2*tmp)%__q;
      commons::tools::set_1coeff_mpz<Pq>(__pol_dec[0], tmp2, 0);
      commons::tools::set_1coeff_mpz<Pq>(__pol_dec[1], tmp, 0);
      __pol_dec[0].ntt_pow_phi();  // gamma*t/mtilde
      __pol_dec[1].ntt_pow_phi();
      __sk[0] = *s;  // s

      __pol_dec_float[0] = 0;
      mpz_invert(tmp2.get_mpz_t(), tmpzmt.get_mpz_t(), __q.get_mpz_t());
      commons::tools::set_1coeff_mpz<Pq>(__pol_dec_float[0], tmp2, 0);
      __pol_dec_float[0].ntt_pow_phi(); // 1/mtilde
      for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
        __t_ov_qi[cmQ] = static_cast<main_greater_double>(1.0L)
                         *static_cast<main_greater_double>
                         (parameters::t)/static_cast<main_greater_double>(Pq::base.params.P[cmQ+__imQ]);
    }

    // SMRmq
    {
      mpz_class tmp, modb;
      for(uint cmB = 0; cmB < __nmB; cmB++)
        {
          mpz_set_ui(modb.get_mpz_t(), Pb::base.params.P[__imB+cmB]); // bi
          mpz_invert(tmp.get_mpz_t(), tmpzmt.get_mpz_t(),
                     modb.get_mpz_t()); // 1/mtilde mod bi
          tmp = (__q*tmp) % modb; // (q/mtilde) mod bi
          __delta_b_tilde_mul[cmB] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
        }
      mpz_invert(tmp.get_mpz_t(), tmpzmt.get_mpz_t(),
                 modsk.get_mpz_t()); // 1/mtilde mod msk
      tmp = (__q*tmp) % modsk;
      __delta_sk_tilde_mul[0] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
    }

    // Mult
    {
      mpz_class tmp, tmp2, modb, modq;
      __polb_mul[0] = 0;
      auto* itb0 = __polb_mul[0].begin();
      for(uint cmB = 0; cmB < __nmB; cmB++, itb0 += degree)
        {
          mpz_set_ui(modb.get_mpz_t(), Pb::base.params.P[__imB+cmB]); // bi

          tmp = -__b/modb; // -b/bi
          for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
            {
              mpz_set_ui(modq.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]);
              tmp2 = tmp*modq; // -(b/bi)*qj
              mpz_invert(tmp2.get_mpz_t(), tmp2.get_mpz_t(), modb.get_mpz_t());
              __M_q_mul_b_1[cmB][cmQ] = static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t()));
            }
          mpz_invert(tmp2.get_mpz_t(), tmpzmt.get_mpz_t(),
                     modb.get_mpz_t()); // 1/mtilde mod bi
          tmp = (__q*tmp2)%modb; // (q/mtilde) mod bi
          __delta_b_tilde_mul[cmB] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
          tmp = __b*__q/modb; // q*b/bi
          mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(), modb.get_mpz_t());
          tmp = (tmp*tmpt)%modb; // (t*bi/b/q) mod bi

          //*(__polb_mul[0].begin()+cmB*Pb::degree) = mpz_get_ui(tmpz1.get_mpz_t())%__polb_mul[0].get_modulus(cmB);
          std::fill(itb0, itb0+degree,
                    static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
        }
      //__polb_mul[0].ntt_pow_phi();

      // sk
      tmp2 = -__q*__b;
      mpz_invert(tmp2.get_mpz_t(), tmp2.get_mpz_t(),
                 modsk.get_mpz_t()); // -1/QB mod msk
      tmp = (tmp2*tmpt)%modsk; // -t/QB mod msk
      __polsk_mul[0] = 0;

      auto* itsk0 = __polsk_mul[0].begin();
      std::fill(itsk0, itsk0+degree,
                static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

      mpz_invert(tmp2.get_mpz_t(), tmpzmt.get_mpz_t(),
                 modsk.get_mpz_t()); // 1/mtilde mod msk
      tmp = (__q*tmp2)%modsk;
      __delta_sk_tilde_mul[0] = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
      tmp = __b%modsk;

      for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
        {
          mpz_set_ui(modq.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]);
          tmp2 = tmp*modq;
          mpz_invert(tmp2.get_mpz_t(), tmp2.get_mpz_t(), modsk.get_mpz_t());
          __M_q_mul_sk_1[cmQ] = static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t()));
        }
      for(uint cmB = 0; cmB < __nmB; cmB++)
        {
          mpz_set_ui(modb.get_mpz_t(), Pb::base.params.P[__imB+cmB]);
          mpz_invert(tmp2.get_mpz_t(), modb.get_mpz_t(), modsk.get_mpz_t());
          __M_b_mul_sk[cmB] = static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t()));
        }

      // q
      __polq_mul[0] = 0;
      auto* itq0 = __polq_mul[0].begin();
      for(uint cmQ = 0; cmQ < __nmQ; cmQ++, itq0 += degree)
        {
          tmp = tmpzmt*tmpzmt;
          mpz_set_ui(modq.get_mpz_t(), Pq::base.params.P[__imQ+cmQ]); // qi
          mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(),
                     modq.get_mpz_t()); // 1/mtilde**2 mod qi
          tmp = (tmpt*tmp*__q/modq)%modq; // t*Q/qi/mtilde**2 mod qi
          std::fill(itq0, itq0+degree,
                    static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));

          tmp = __q_ov_qi[cmQ]; // Q/qi
          mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(), modq.get_mpz_t()); // qi/Q
          tmp2 = modq-(__b*tmp)%modq; // -M qi/Q
          __delta_q_sk_mul2[cmQ] = static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t()));
          tmp2 = (tmp2*tmpzmt)%modq; // -M*mtilde*qi/Q
          __delta_q_sk_mul[cmQ] = static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t()));

          for(uint cmB = 0; cmB < __nmB; cmB++)
            {
              mpz_set_ui(modb.get_mpz_t(), Pb::base.params.P[__imB+cmB]);
              tmp2 = (tmp*__b/modb)%modq; // Mj qi/Q
              __M_b_mul_q2[cmQ][cmB] = static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t()));
              tmp2 = (tmp2*tmpzmt)%modq; // Mj*mtilde*qi/Q
              __M_b_mul_q[cmQ][cmB] = static_cast<main_type>(mpz_get_ui(tmp2.get_mpz_t()));
            }
        }
    }

    // evk
    {
      mpz_class modq, tmp;
      Pq* es = commons::mem_manag::alloc_aligned<Pq, 32>(2);
      nfl::uniform unif;
      for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
        {
          es[1] = *s**s;
          commons::tools::mul_mpz_class<Pq>(es[1], __q_ov_qi[cmQ]); // s*s*Qi
          __evka[cmQ] = unif;
          __evka[cmQ].ntt_pow_phi();
          es[0] = nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng);
          es[0].ntt_pow_phi();
          __evkb[cmQ] = es[0]+__evka[cmQ]**s;
          __evkb[cmQ] = es[1]-__evkb[cmQ]; // s*s*Qi-(e+a*s)
          for(uint cmQ2 = 0; cmQ2 < __nmQ; cmQ2++)
            {
              mpz_set_ui(modq.get_mpz_t(), Pq::base.params.P[__imQ+cmQ2]);
              tmp = __q/modq;
              mpz_invert(tmp.get_mpz_t(), tmp.get_mpz_t(), modq.get_mpz_t()); // qj/Q mod qj
              tmp = (tmp*tmpzmt)%modq;  // mtilde/Qj mod qj
              main_type val = static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t()));
              for(uint d = 0; d < degree; d++)
                {
                  commons::tools::mul_val<Pq>(__evka[cmQ], val, cmQ2, d);
                  commons::tools::mul_val<Pq>(__evkb[cmQ], val, cmQ2, d);
                }
            }
        }
      commons::mem_manag::free_aligned<Pq>(2, es);
    }

    mpfr_set_default_prec (10000);
    // HPS-based extensions
    {
      for (size_t i = 0; i < __nmQ; i++)
        {
          mpz_class mod (Pq::base.params.P[__imQ+i]);
          mpfr_t mpfr_q_i, w;
          mpfr_init (mpfr_q_i);
          mpfr_set_z (mpfr_q_i, mod.get_mpz_t (), MPFR_RNDN);
          mpfr_ui_div (mpfr_q_i, 1, mpfr_q_i, MPFR_RNDN);
          __M_q_i_inv[i] = mpfr_get_ld (mpfr_q_i, MPFR_RNDN);

          mpfr_clear (mpfr_q_i);
        }

      for (size_t i = 0; i < __nmQ; i++)
        {
          for (size_t j = 0; j < __nmB; j++)
            {
              mpz_class pj (Pb::base.params.P[__imB+j]);
              mpz_class qi_pj (__q_ov_qi[i] % pj);
              __M_q_i[j][i] = mpz_get_ui (qi_pj.get_mpz_t ());
            }
        }

      for (size_t i = 0; i < __nmQ; i++)
	{
	  for (size_t j = 0; j < __nmB; j++)
	    {
	      mpz_class pj(Pb::base.params.P[__imB+j]);
	      mpz_class Qj_tilde(__q * __b / pj);
	      mpz_invert(Qj_tilde.get_mpz_t(), Qj_tilde.get_mpz_t(),
			 pj.get_mpz_t());
	      mpz_class Qj_tilde_t = Qj_tilde * tmpt;
	      mpz_class new_M_q_i = (mpz_class(__M_q_i[j][i]) * Qj_tilde_t) % pj;
	      __M_q_i_1[j][i] = mpz_get_ui(new_M_q_i.get_mpz_t());
	    }
	}

      for (size_t i = 0; i < __nmB; i++)
        {
          mpz_class pi (Pb::base.params.P[__imB+i]);
          mpz_class q_pi (__q % pi);
          __M_q[i] = mpz_get_ui (q_pi.get_mpz_t ());
        }

      for (size_t i = 0; i < __nmB; i++)
	{
	  mpz_class pi(Pb::base.params.P[__imB+i]);
	  mpz_class Qi_tilde(__q * __b / pi);
	  mpz_invert(Qi_tilde.get_mpz_t(), Qi_tilde.get_mpz_t(),
		     pi.get_mpz_t());
	  mpz_class Qi_tilde_t = Qi_tilde * tmpt;
	  mpz_class new_M_q = (mpz_class(__M_q[i]) * Qi_tilde_t) % pi;
	  __M_q_1[i] = mpz_get_ui(new_M_q.get_mpz_t());
	}
	  
      for (size_t i = 0; i < __nmB; i++)
        {
          mpz_class mod (Pb::base.params.P[__imB+i]);
          mpfr_t mpfr_b_i;
          mpfr_init (mpfr_b_i);
          mpfr_set_z (mpfr_b_i, mod.get_mpz_t (), MPFR_RNDN);
          mpfr_ui_div (mpfr_b_i, 1, mpfr_b_i, MPFR_RNDN);
          __M_b_i_inv[i] = mpfr_get_ld (mpfr_b_i, MPFR_RNDN);

          mpfr_clear (mpfr_b_i);
        }

      for (size_t i = 0; i < __nmB; i++)
        {
          for (size_t j = 0; j < __nmQ; j++)
            {
              mpz_class qj (Pq::base.params.P[__imQ+j]);
              mpz_class pi_qj (__b_ov_bi[i] % qj);
              __M_b_i[j][i] = mpz_get_ui (pi_qj.get_mpz_t ());
	      
	      mpz_class qj_ov_q(__q_ov_qi[j] % qj);
	      mpz_invert(qj_ov_q.get_mpz_t(),
			 qj_ov_q.get_mpz_t(),
			 qj.get_mpz_t());
	      mpz_class new_M_b_i = (__M_b_i[j][i] * qj_ov_q) % qj;
	      __M_b_i_1[j][i] = mpz_get_ui(new_M_b_i.get_mpz_t());
            }
        }

      for (size_t i = 0; i < __nmQ; i++)
        {
          mpz_class qi (Pq::base.params.P[__imQ+i]);
          mpz_class p_qi (__b % qi);
          __M_b[i] = mpz_get_ui (p_qi.get_mpz_t ());
	  
	  mpz_class qi_ov_q(__q_ov_qi[i] % qi);
	  mpz_invert(qi_ov_q.get_mpz_t(),
		     qi_ov_q.get_mpz_t(),
		     qi.get_mpz_t());
	  mpz_class new_M_b = (__M_b[i] * qi_ov_q) % qi;
	  __M_b_1[i] = mpz_get_ui(new_M_b.get_mpz_t());
        }

      for (size_t i = 0; i < __nmQ; i++)
	{
	  mpz_class qi(Pq::base.params.P[__imQ+i]);
	  mpz_class q_ov_qi_mod_qi = __q_ov_qi[i] % qi;
	  
	  main_type ul_q_ov_qi_mod_qi = mpz_get_ui
	    (q_ov_qi_mod_qi.get_mpz_t());
	  
	  std::fill(&__polq_mul_1[0](i, 0),
		    &__polq_mul_1[0](i, 0) + degree,
		    ul_q_ov_qi_mod_qi);
	}

      for (size_t i = 0; i < __nmQ; i++)
        {
          mpz_class inv_b (__b);
          inv_b %= mpz_class (Pq::base.params.P[__imQ+i]);
          mpz_invert (inv_b.get_mpz_t (),
                      inv_b.get_mpz_t (),
                      mpz_class (Pq::base.params.P[__imQ+i]).get_mpz_t ());

          mpz_class num;
          num = inv_b * tmpt * __b;
          mpz_class den (Pq::base.params.P[__imQ+i]);
          mpz_class omegai;
          mpfr_t z, y;
          mpfr_init (z);
          mpfr_init (y);
          mpfr_set_z (z, num.get_mpz_t (), MPFR_RNDN);
          mpfr_set_z (y, den.get_mpz_t (), MPFR_RNDN);
          mpfr_div (z, z, y, MPFR_RNDN); // num / den
          mpfr_get_z (omegai.get_mpz_t (), z, MPFR_RNDN);
          mpfr_set_z (y, omegai.get_mpz_t (), MPFR_RNDN);
          mpfr_sub (y, z, y, MPFR_RNDN); // num / den - round (num / den)

          thetas[i] = mpfr_get_ld (y, MPFR_RNDN);
          // std::cout << "thetas[" << i << "] = " << thetas[i] << "\n";
          // std::cout << "omegas[" << i << "] = " << omegai << "\n";

          for (size_t j = 0; j < __nmB; j++)
            {
              mpz_class omegaij;
	      mpz_class pj(Pb::base.params.P[j+__imB]);
	      mpz_class tilde_p_j(__b / pj);
	      tilde_p_j %= pj;
	      mpz_invert(tilde_p_j.get_mpz_t(), tilde_p_j.get_mpz_t(),
			 pj.get_mpz_t());
              omegaij = (omegai*tilde_p_j) % pj;
              omegas[j][i] = static_cast<main_type>(mpz_get_ui(omegaij.get_mpz_t()));
            }

          mpfr_clear (z);
          mpfr_clear (y);
        }

      for (size_t i = 0; i < __nmB; i++)
        {
	  mpz_class pi(Pb::base.params.P[i+__imB]);
	  mpz_class tilde_p_i(__b / pi);
	  tilde_p_i %= pi;
	  mpz_invert(tilde_p_i.get_mpz_t(),
		     tilde_p_i.get_mpz_t(),
		     pi.get_mpz_t());
          lambdas[i] = static_cast<main_type>(mpz_get_ui(tilde_p_i.get_mpz_t ()));
          // std::cout << "lambdas[" << i << "] = " << lambdas [i] << "\n";
        }
    }

    //Encryption for mul_float
    {
      __pol_enc_float[2] = 0;
      __pol_enc_float[3] = 0;

      mpz_class modq, tmp;

      auto* it2 = __pol_enc_float[2].begin(), *it3 = __pol_enc_float[3].begin();

      for(uint cmQ = 0; cmQ < __nmQ; cmQ++, it2 += degree, it3 += degree)
        {
          // q
          mpz_set_ui(tmp.get_mpz_t(), invQi[cmQ]);
          mpz_set_ui(modq.get_mpz_t(), Pq::base.params.P[cmQ+__imQ]);
          tmp = (tmp) % modq; // qi/q mod qi

          std::fill(it3, it3+degree, static_cast<main_type>(mpz_get_ui(tmp.get_mpz_t())));
        }
      commons::tools::set_1coeff_mpz<Pq>(__pol_enc_float[2], Delta, 0);
      __pol_enc_float[2].ntt_pow_phi();
      __pol_enc_float[2] = __pol_enc_float[2] * __pol_enc_float[3];
      __pol_enc_float[0] = *a * __pol_enc_float[3];      // a
      __pol_enc_float[1] = *b * __pol_enc_float[3];      // b
      __qi_ov_q[0] = __pol_enc_float[3];
      // std::cout << "__qi_ov_q = " << __qi_ov_q << "\n";

      for (size_t i = 0; i < __nmB; i++)
        {
          mpz_class mod (Pb::base.params.P[i+__imB]);
          tmp = (__b / mod) % mod;
          mpz_invert (tmp.get_mpz_t (),
                      tmp.get_mpz_t (),
                      mod.get_mpz_t ());

          for (size_t j = 0; j < degree; j++)
            {
              __bi_ov_b[0] (i, j) = static_cast<main_type> (mpz_get_ui(tmp.get_mpz_t ()));
            }
        }
    }

    {
      // Decryption for HPS
      for (size_t i = 0; i < __nmQ; i++)
        {
          mpz_class mod (Pq::base.params.P[i+__imQ]);
          mpz_class num = tmpt;
          mpz_class den = mod;

          mpz_class omegai;
          mpfr_t z, y;
          mpfr_init(z);
          mpfr_init(y);
          mpfr_set_z(z, num.get_mpz_t(), MPFR_RNDN);
          mpfr_set_z(y, den.get_mpz_t(), MPFR_RNDN);
          mpfr_div(z, z, y, MPFR_RNDN);

          dec_thetas[i] = mpfr_get_ld(z, MPFR_RNDN);

          mpfr_clear(z);
          mpfr_clear(y);
        }
    }

    {
      mpz_class modq, tmp;
      Pq* es = commons::mem_manag::alloc_aligned<Pq, 32>(2);
      nfl::uniform unif;
      for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
        {
          es[1] = *s * *s * __qi_ov_q[0];
          commons::tools::mul_mpz_class<Pq>(es[1], __q_ov_qi[cmQ]); // s*s*Qi
          __evka_float[cmQ] = unif;
          __evka_float[cmQ].ntt_pow_phi();
          es[0] = nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng);
          es[0].ntt_pow_phi();
          __evkb_float[cmQ] = es[0]+__evka_float[cmQ]**s;
          __evkb_float[cmQ] = es[1]-__evkb_float[cmQ]; // s*s*Qi*Qi^(-1)-(e+a*s)

	  __evka_float[cmQ] = __evka_float[cmQ] * __qi_ov_q[0];
	  __evkb_float[cmQ] = __evkb_float[cmQ] * __qi_ov_q[0];
        }
      commons::mem_manag::free_aligned<Pq>(2, es);

    }
  }


// computes the infinity norm of the noise of c
  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  size_t Context<Pq, Pb, Psk, Pmt, parameters>::size_noise(const _Ciphertext &c,
      const Pt &m, bool input_in_ntt) const
  {
    constexpr size_t degree = Pq::degree;
    size_t res = 0;
    mpz_class coeff, tmpz;
    Pq* tmp = commons::mem_manag::alloc_aligned<Pq, 32>(3);
    m.template get<Pq>(tmp[0]);

    tmp[1] = c.__cq[0];
    tmp[2] = c.__cq[1];
    if(input_in_ntt)
      {
        tmp[1].invntt_pow_invphi();
        tmp[2].invntt_pow_invphi();
      }
    tmp[0] = tmp[0]*__pol_noise[0]; // xi_q(m*Delta)
    if (__red)
      {
        tmp[1] = tmp[1]*__pol_noise[1]; // xi_q(c0/mtilde)
      }
    tmp[2].ntt_pow_phi();
    tmp[2] = tmp[2]*__pol_noise[2]; // xi_q(c1*s/mtilde)
    tmp[2].invntt_pow_invphi();
    tmp[0] = tmp[1]+tmp[2]-tmp[0];  // c0+c1*s-Delta*m

    for(uint d = 0; d < degree; d++)
      {
        coeff = 0;
        for(uint cmQ = 0; cmQ < __nmQ; cmQ++)
          {
            mpz_set_ui(tmpz.get_mpz_t(), tmp[0](cmQ, d));
            coeff += tmpz*__q_ov_qi[cmQ];
          }
        coeff = coeff%__q;
        if(coeff > (__q-1)/2)
          coeff -= __q;
        coeff = abs(coeff);
        size_t log_coeff = mpz_sizeinbase(coeff.get_mpz_t(), 2);
        if(log_coeff > res)
          res = log_coeff;
      }
    commons::mem_manag::free_aligned<Pq>(3, tmp);
    return res;
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  inline void Context<Pq, Pb, Psk, Pmt, parameters>::extend(_Ciphertext &c)
  {
    // extends cq to other bases
    if (__red)
      {
        c.fast_base_conversion_red(
          __M_q_enc_b,
          __M_q_enc_sk
          , __M_q_enc_tilde
        );
      }
    else
      {
        c.fast_base_conversion(
          __M_q_enc_b,
          __M_q_enc_sk
        );
      }
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  inline void Context<Pq, Pb, Psk, Pmt, parameters>::extend_float(_Ciphertext &c)
  {
    // extends cq to other bases
    c.base_conversion_float(__M_q_i_inv,
                            __M_q_i,
                            __M_q);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  inline void Context<Pq, Pb, Psk, Pmt, parameters>::extend_float_1(_Ciphertext &c)
  {
    // extends cq to other bases
    c.base_conversion_float(__M_q_i_inv,
                            __M_q_i_1,
                            __M_q_1);
  }  

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  inline void Context<Pq, Pb, Psk, Pmt, parameters>::extend_back_float(
    _Ciphertext &c)
  {
    // extends cb to other bases
    c.base_conversion_back_float(__M_b_i_inv,
                                 __M_b_i_1,
                                 __M_b_1);
  }


  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::encrypt(_Ciphertext &c, const Pq &m,
      bool output_in_ntt)
  {
    Pq* e = commons::mem_manag::alloc_aligned<Pq, 32, nfl::gaussian<uint8_t, main_type, 2>>
            (2, nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng));
    Pq* u = commons::mem_manag::alloc_aligned<Pq, 32, nfl::gaussian<uint8_t, main_type, 2>>
            (1, nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng));
    u[0].ntt_pow_phi();

    c.__cq[0] = u[0]*__pol_enc[1];
    c.__cq[1] = m*__pol_enc[2]+e[0]*__pol_enc[3];
    c.__cq[1].ntt_pow_phi();
    c.__cq[0] = c.__cq[0]+c.__cq[1];

    c.__cq[1] = u[0]*__pol_enc[0];
    e[1] = e[1]*__pol_enc[3];
    e[1].ntt_pow_phi();
    c.__cq[1] = c.__cq[1]+e[1];

    if(not output_in_ntt)
      {
        c.__cq[0].invntt_pow_invphi();
        c.__cq[1].invntt_pow_invphi();
      }
    commons::mem_manag::free_aligned<Pq>(2, e);
    commons::mem_manag::free_aligned<Pq>(1, u);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::encrypt_float(_Ciphertext &c,
      const Pq &m, bool output_in_ntt)
  {
    Pq* e = commons::mem_manag::alloc_aligned<Pq, 32, nfl::gaussian<uint8_t, main_type, 2>>
            (2, nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng));
    Pq* u = commons::mem_manag::alloc_aligned<Pq, 32, nfl::gaussian<uint8_t, main_type, 2>>
            (1, nfl::gaussian<uint8_t, main_type, 2>(&parameters::g_prng));
    u[0].ntt_pow_phi();

    c.__cq[0] = u[0]*__pol_enc_float[1];
    c.__cq[1] = m*__pol_enc_float[2]+e[0]*__pol_enc_float[3];
    c.__cq[1].ntt_pow_phi();
    c.__cq[0] = c.__cq[0]+c.__cq[1];

    c.__cq[1] = u[0]*__pol_enc_float[0];
    e[1] = e[1]*__pol_enc_float[3];
    e[1].ntt_pow_phi();
    c.__cq[1] = c.__cq[1]+e[1];

    if(not output_in_ntt)
      {
        c.__cq[0].invntt_pow_invphi();
        c.__cq[1].invntt_pow_invphi();
      }
    commons::mem_manag::free_aligned<Pq>(2, e);
    commons::mem_manag::free_aligned<Pq>(1, u);
  }


// floating point variant of full RNS decryption (cf. Section 3.5 in paper)
  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::decrypt_float_old(Pt &mdec,
      _Ciphertext &c)
  {
    constexpr size_t degree = Pq::degree;
    c.__cq[0] = c.__cq[0]+c.__cq[1]*__sk[0];
    if (__red)
      {
        c.__cq[0] = c.__cq[0]*__pol_dec_float[0];
      }
    c.__cq[0].invntt_pow_invphi();

    for(uint d = 0; d < degree; d++)
      {
        main_greater_double acc = 0;
        for(uint cm = 0; cm < __nmQ; cm++)
          acc += static_cast<main_greater_double>(c.__cq[0](cm, d))*__t_ov_qi[cm];

        mdec.__data[d] = static_cast<typename Pt::type_t>(round(
                           acc))&parameters::mask_t;
      }
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::decrypt_float(Pt &mdec,
      _Ciphertext &c)
  {
    constexpr size_t degree = Pq::degree;
    c.__cq[0] = c.__cq[0]+c.__cq[1]*__sk[0];
    c.__cq[0].invntt_pow_invphi();

    for (size_t i = 0; i < degree; i++)
      {
        main_greater_double v = 0;

        for (size_t j = 0; j < __nmQ; j++)
          {
            main_greater_double cij = c.__cq[0](j, i);

            if (c.__cq[0](j, i) > Pq::base.params.P[j+__imQ]/2)
              cij -= Pq::base.params.P[j+__imQ];

            v += cij * dec_thetas[j];
          }

        mdec.__data[i] = static_cast<typename Pt::type_t>
	  (static_cast<main_signed_type>(std::roundl (v)) & parameters::mask_t);
      }
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::decrypt_float_not_ntt(Pt &mdec,
      _Ciphertext &c)
  {
    constexpr size_t degree = Pq::degree;
    c.__cq[1].ntt_pow_phi();
    c.__cq[1] = c.__cq[1]*__sk[0];
    c.__cq[1].invntt_pow_invphi();
    c.__cq[0] = c.__cq[0]+c.__cq[1];

    for (size_t i = 0; i < degree; i++)
      {
        main_greater_double v = 0;

        for (size_t j = 0; j < __nmQ; j++)
          {
            main_greater_double cij = c.__cq[0](j, i);

            if (c.__cq[0](j, i) > Pq::base.params.P[j+__imQ]/2)
              cij -= Pq::base.params.P[j+__imQ];

            v += cij * dec_thetas[j];
          }

        mdec.__data[i] = static_cast<typename Pt::type_t>
	  (static_cast<main_signed_type>(std::roundl (v)) & parameters::mask_t);
      }
  }
  

// RNS decryption; Algo. 1, variant with gamma = power of 2 (cf. Section 3.5)
// computes the fast base conversion from base q to base with one modulus gamma*t
// then, recovers the error (centered residue modulo gamma) and corrects the residue modulo t
  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::decrypt(Pt &mdec, _Ciphertext &c)
  {
    constexpr size_t degree = Pq::degree;
    c.__cq[0] = c.__cq[0]+c.__cq[1]*__sk[0];
    c.__cq[0] = c.__cq[0]* (__red ?
                            __pol_dec[0] :
                            __pol_dec[1]);
    c.__cq[0].invntt_pow_invphi();

    auto* it_mdec = mdec.begin();
    for(uint d = 0; d < degree; d++)
      {
        typename parameters::type_gamma_t acc = parameters::gamma_half;
        auto* it = c.__cq[0].begin()+d;
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

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::decrypt_not_ntt(Pt &mdec, _Ciphertext &c)
  {
    constexpr size_t degree = Pq::degree;
    c.__cq[1].ntt_pow_phi();
    c.__cq[1] = c.__cq[1]*__sk[0];
    c.__cq[1].invntt_pow_invphi();
    c.__cq[0] = c.__cq[0]+c.__cq[1];
    c.__cq[0] = c.__cq[0]* (__red ?
                            __pol_dec[0] :
                            __pol_dec[1]);

    auto* it_mdec = mdec.begin();
    for(uint d = 0; d < degree; d++)
      {
        typename parameters::type_gamma_t acc = parameters::gamma_half;
        auto* it = c.__cq[0].begin()+d;
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
  
// cf. Algo. 3 in paper
  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::mult(_Ciphertext &c1, _Ciphertext &c2,
      bool output_in_ntt)
  {
    // step 0
    this->extend(c1);
    this->extend(c2);
    if (__red)
      {
        // step 1
        c1.SmMRq(__delta_b_tilde_mul, __delta_sk_tilde_mul);
        c2.SmMRq(__delta_b_tilde_mul, __delta_sk_tilde_mul);
      }
    // step 2
    c1.ntt();
    c2.ntt();
    c1.mul(__coeff2[0], c2, __polq_mul[0], __polb_mul[0], __polsk_mul[0]);
    c1.inv_ntt0();
    c1.inv_ntt1();
    __coeff2->inv_ntt0();

    // steps 3, 4

    // q -> b
    commons::base_conversions::fast_1m_to_1m<Pb, Pq>(c1.__cb[0], c1.__cq[0],
        __M_q_mul_b_1);
    commons::base_conversions::fast_1m_to_1m<Pb, Pq>(c1.__cb[1], c1.__cq[1],
        __M_q_mul_b_1);
    commons::base_conversions::fast_1m_to_1m<Pb, Pq>(__coeff2->__cb[0],
        __coeff2->__cq[0], __M_q_mul_b_1);

    // b,q -> m_sk
    commons::base_conversions::fast_2m_to_11<Psk, Pq, Pb>(c1.__csk[0], c1.__cq[0],
        c1.__cb[0], __M_q_mul_sk_1, __M_b_mul_sk);
    commons::base_conversions::fast_2m_to_11<Psk, Pq, Pb>(c1.__csk[1], c1.__cq[1],
        c1.__cb[1], __M_q_mul_sk_1, __M_b_mul_sk);
    commons::base_conversions::fast_2m_to_11<Psk, Pq, Pb>(__coeff2->__csk[0],
        __coeff2->__cq[0], __coeff2->__cb[0], __M_q_mul_sk_1, __M_b_mul_sk);


    // here, the approximate rounding is in b+m_sk
    // it just remains to convert it exactly in base q
    c1.__cq[0] = 0;
    c1.__cq[1] = 0;
    __coeff2->__cq[0] = 0;

    // b, m_sk -> q
    commons::base_conversions::sk<Pq, Pb, Psk>(c1.__cq[0], c1.__cb[0], c1.__csk[0],
        __M_b_mul_q, __delta_q_sk_mul);
    commons::base_conversions::sk<Pq, Pb, Psk>(c1.__cq[1], c1.__cb[1], c1.__csk[1],
        __M_b_mul_q, __delta_q_sk_mul);
    commons::base_conversions::sk<Pq, Pb, Psk>(__coeff2->__cq[0], __coeff2->__cb[0],
        __coeff2->__csk[0], __M_b_mul_q2, __delta_q_sk_mul2);

    if(output_in_ntt)
      {
        c1.__cq[0].ntt_pow_phi();
        c1.__cq[1].ntt_pow_phi();
      }
  }



  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::mult_float_step1(_Ciphertext &c1,
      _Ciphertext &c2, bool output_in_ntt)
  {
    // step 0
    this->extend_float(c1);
    this->extend_float_1(c2);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::mult_float_step2(_Ciphertext &c1,
      _Ciphertext &c2, bool output_in_ntt)
  {
    // step 2
    c1.ntt_float();
    c2.ntt_float();
    c1.mul_float(__coeff2[0], c2, __polq_mul_1[0]);
    c1.inv_ntt0_float();
    c1.inv_ntt1_float();
    __coeff2->inv_ntt0_float();
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::mult_float_step3(_Ciphertext &c1,
      _Ciphertext &c2, bool output_in_ntt)
  {
    // steps 3, 4
    commons::base_conversions::fast_div_round<Pb, Pq>
    (c1.__cb[0], c1.__cq[0],
     thetas, lambdas, omegas);
    commons::base_conversions::fast_div_round<Pb, Pq>
    (c1.__cb[1], c1.__cq[1],
     thetas, lambdas, omegas);
    commons::base_conversions::fast_div_round<Pb, Pq>
    (__coeff2->__cb[0], __coeff2->__cq[0],
     thetas, lambdas, omegas);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::mult_float_step4(_Ciphertext &c1,
      _Ciphertext &c2, bool output_in_ntt)
  {
    // b -> q
    this->extend_back_float(c1);
    commons::base_conversions::fast_1m_to_1m_float<Pq, Pb>
    (__coeff2->__cq[0], __coeff2->__cb[0], this->__M_b_i_inv,
     this->__M_b_i, this->__M_b);
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::mult_float(_Ciphertext &c1,
      _Ciphertext &c2, bool output_in_ntt)
  {
    mult_float_step1(c1, c2, output_in_ntt);
    mult_float_step2(c1, c2, output_in_ntt);
    mult_float_step3(c1, c2, output_in_ntt);
    mult_float_step4(c1, c2, output_in_ntt);

    if(output_in_ntt)
      {
        c1.__cq[0].ntt_pow_phi();
        c1.__cq[1].ntt_pow_phi();
      }
  }


  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::relinearisation(_Ciphertext &c1,
      bool output_in_ntt)
  {
    constexpr size_t degree = Pq::degree;
    __decomp[1] = 0;
    __decomp[2] = 0;
    auto* it = __coeff2->__cq[0].begin();
    for(uint cmQ = 0; cmQ < __nmQ; cmQ++, it+=degree)
      {
        __decomp[0] = 0;
        commons::tools::set_partial(__decomp[0], it, it+degree);
        __decomp[0].ntt_pow_phi();
        __decomp[1] = __decomp[1]+__decomp[0]*__evkb[cmQ];
        __decomp[2] = __decomp[2]+__decomp[0]*__evka[cmQ];
      }

    if(not output_in_ntt)
      {
        __decomp[1].invntt_pow_invphi();
        __decomp[2].invntt_pow_invphi();
      }
    c1.__cq[0] = c1.__cq[0]+__decomp[1];
    c1.__cq[1] = c1.__cq[1]+__decomp[2];
  }

  template<class Pq, class Pb, class Psk, class Pmt, class parameters>
  void Context<Pq, Pb, Psk, Pmt, parameters>::relinearisation_float(
    _Ciphertext &c1, bool output_in_ntt)
  {
    constexpr size_t degree = Pq::degree;
    __decomp[1] = 0;
    __decomp[2] = 0;
    auto* it = __coeff2->__cq[0].begin();
    for(uint cmQ = 0; cmQ < __nmQ; cmQ++, it+=degree)
      {
        __decomp[0] = 0;
        commons::tools::set_partial(__decomp[0], it, it+degree);
        __decomp[0].ntt_pow_phi();
        __decomp[1] = __decomp[1]+__decomp[0]*__evkb_float[cmQ];
        __decomp[2] = __decomp[2]+__decomp[0]*__evka_float[cmQ];
      }

    if(not output_in_ntt)
      {
        __decomp[1].invntt_pow_invphi();
        __decomp[2].invntt_pow_invphi();
      }
    c1.__cq[0] = c1.__cq[0]+__decomp[1];
    c1.__cq[1] = c1.__cq[1]+__decomp[2];
  }


}


#endif /* rns_h */
