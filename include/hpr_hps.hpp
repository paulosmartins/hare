#ifndef __HPR_HPS_HPP__
#define __HPR_HPS_HPP__

template<typename nfl_params_t, size_t degree, size_t nmoduli,
	 size_t ndigits, typename gamma_t_parameters_t, size_t bs_mt>
struct context_hpr_hps_t
{
  typedef nfl::poly<nfl_params_t, degree, nmoduli, 0> Pq;
  typedef nfl::poly<nfl_params_t, degree, nmoduli+1, nmoduli> Pb;
  //Psk is unused. We need this definition because it's the same
  //class that implements HPR and HPR-HPS, but with different methods
  typedef nfl::poly<nfl_params_t, degree, 1, 2*nmoduli+1> Psk;
  
  using main_type = typename Pq::value_type;
  using main_greater_type = typename Pq::greater_value_type;
  using main_signed_type = typename Pq::signed_value_type;

  typedef typename std::conditional<bs_mt == 8,
				    uint8_t,
				    typename
				    std::conditional<bs_mt == 16,
						     uint16_t,
						     std::nullptr_t>::type>::type
  type_mtilde;
  
  typedef typename std::conditional<bs_mt == 8,
				    int8_t,
				    typename
				    std::conditional<bs_mt == 16,
						     int16_t,
						     std::nullptr_t>::type>::type
  signed_type_mtilde;
  
  typedef hpr::commons::raw_pol::PolMT<main_type, main_signed_type,
				       type_mtilde, signed_type_mtilde,
				       bs_mt, degree> Pmt;

  typedef hpr::commons::raw_pol::PolT<gamma_t_parameters_t, degree> Pt;

  typedef hpr::context::Context<Pq, Pb, Psk, Pmt, gamma_t_parameters_t, ndigits> context_hpr_t;

  typedef hpr::digit::poly_digit<Pq, Pb, Psk, gamma_t_parameters_t, Pmt> poly_digit;
  typedef hpr::poly_hpr<poly_digit, ndigits, Pmt> poly_hpr;

  typedef hpr::commons::key::key_chain<poly_hpr, gamma_t_parameters_t> key_hpr_t;

  context_hpr_t *cont_HPR;

  context_hpr_hps_t()
    : cont_HPR(new context_hpr_t())
  {
  }

  context_hpr_hps_t(const context_hpr_hps_t& other)
    : cont_HPR(new context_hpr_t())
  {
    *cont_HPR = *other.cont_HPR;
  }

  context_hpr_hps_t& operator=(const context_hpr_hps_t& other)
  {
    *cont_HPR = *other.cont_HPR;
  }

  ~context_hpr_hps_t()
  {
    delete cont_HPR;
  }

  struct plaintext_t
  {
    Pt *_p;

    plaintext_t()
      : _p(new Pt())
    {
    }

    plaintext_t(const plaintext_t& other)
      : _p(new Pt())
    {
      *_p = *other._p;
    }

    plaintext_t &operator=(const plaintext_t& other)
    {
      *_p = *other._p;
      return *this;
    }

    ~plaintext_t()
    {
      delete _p;
    }
  };

  struct ciphertext_t
  {
    poly_hpr *_c;

    ciphertext_t()
      : _c(new poly_hpr())
    {
    }

    ciphertext_t(const ciphertext_t &other)
      : _c(new poly_hpr())
    {
      *_c = *other._c;
    }

    ciphertext_t &operator=(const ciphertext_t &other)
    {
      *_c = *other._c;
      return *this;
    }

    ~ciphertext_t()
    {
      delete _c;
    }
  };

  struct keys_t
  {
    key_hpr_t *_keys;

    keys_t()
      : _keys(new key_hpr_t())
    {
    }

    keys_t(const keys_t &other)
      : _keys(new key_hpr_t())
    {
      *_keys = *other._keys;
    }

    keys_t& operator=(const keys_t &other)
    {
      *_keys = *other._keys;
      return *this;
    }

    ~keys_t()
    {
      delete _keys;
    }
  };

  template<typename iterator_t>
  void import_it(plaintext_t &p, iterator_t it)
  {
    for (size_t i = 0; i < degree; i++)
      {
	p._p->__data[i] = *it++;
      }
  }

  template<typename iterator_t>
  void export_it(iterator_t it, plaintext_t &p)
  {
    for (size_t i = 0; i < degree; i++)
      {
	*it++ = p._p->__data[i];
      }
  }

  void generate_keys(keys_t &k)
  {
    cont_HPR->__be_tool_set->generate();
    k._keys->generate(cont_HPR->__be_tool_set);
  }

  void do_precomputations(keys_t &k)
  {
    cont_HPR->set_precomp(k._keys);
    cont_HPR->set_precomp_relin();
  }

  void encrypt(ciphertext_t &c, plaintext_t &p)
  {
    c._c->set_mat_set_conv(cont_HPR->__be_tool_set);
    cont_HPR->encrypt(c._c, p._p);
    cont_HPR->postproc_encrypt(c._c);
  }

  void decrypt(plaintext_t &p, ciphertext_t &c)
  {
    cont_HPR->preproc_decrypt(c._c);
    cont_HPR->decrypt(p._p, c._c);
  }

  void multiply(ciphertext_t &c1, ciphertext_t &c2)
  {
    cont_HPR->multiplication_float(c1._c, c2._c);
  }

  void relinearise(ciphertext_t &c1)
  {
    cont_HPR->relinearisation_float(c1._c);
  }

  void add(ciphertext_t &c1, ciphertext_t &c2)
  {
    cont_HPR->addition(c1._c, c2._c);
  }

  size_t get_degree()
  {
    return degree;
  }
};

#endif
