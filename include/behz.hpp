#ifndef __BEHZ_HPP__
#define __BEHZ_HPP__

template
<typename nfl_params_t, size_t degree, size_t nmoduli,
 typename gamma_t_parameters_t, size_t bs_mt>
struct context_behz_t
{
  typedef nfl::poly<nfl_params_t, degree, nmoduli, 0> Pq;
  typedef nfl::poly<nfl_params_t, degree, nmoduli, nmoduli> Pb;
  typedef nfl::poly<nfl_params_t, degree, 1, 2*nmoduli> Psk;
  
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
  
  typedef rns::commons::fast_pol::PolMT<main_type, main_signed_type,
					type_mtilde, signed_type_mtilde,
					bs_mt, degree> Pmt;

  typedef rns::commons::fast_pol::PolT<gamma_t_parameters_t, degree> P_t;

  typedef rns::Ciphertext<Pq, Pb, Psk, Pmt> ciphertext_rns_t;

  typedef rns::Context<Pq, Pb, Psk, Pmt, gamma_t_parameters_t> context_rns_t;

  context_rns_t *cont_RNS;

  context_behz_t(bool red = true)
    : cont_RNS(new context_rns_t(red))
  {
  }

  context_behz_t(const context_behz_t &other)
    : cont_RNS(new context_rns_t())
  {
    *cont_RNS = *other.cont_RNS;
  }
  
  context_behz_t &operator=(const context_behz_t &other)
  {
    *cont_RNS = *other.cont_RNS;
    return *this;
  }
  
  ~context_behz_t()
  {
    delete cont_RNS;
  }
  
  struct plaintext_t
  {
    P_t *_p; //output
    Pq *_m; //input
    
    plaintext_t()
    {
      _p = rns::commons::mem_manag::alloc_aligned<P_t, 32>(1);
      _m = rns::commons::mem_manag::alloc_aligned<Pq, 32>(1);
    }
    
    plaintext_t(const plaintext_t &other)
    {
      _p = rns::commons::mem_manag::alloc_aligned<P_t, 32>(1);
      _m = rns::commons::mem_manag::alloc_aligned<Pq, 32>(1);
      *_p = *other._p;
      *_m = *other._m;
    }

    plaintext_t& operator=(const plaintext_t &other)
    {
      *_p = *other._p;
      *_m = *other._m;
      return *this;
    }
    
    ~plaintext_t()
    {
      rns::commons::mem_manag::free_aligned<P_t>(1, _p);
      rns::commons::mem_manag::free_aligned<Pq>(1, _m);
    }
  };

  struct ciphertext_t
  {
    ciphertext_rns_t *_c;
    
    ciphertext_t()
    {
      _c = rns::commons::mem_manag::alloc_aligned<ciphertext_rns_t, 32>(1);
    }

    ciphertext_t(const ciphertext_t &other)
    {
      _c = rns::commons::mem_manag::alloc_aligned<ciphertext_rns_t, 32>(1);
      *_c = *other._c;
    }

    ciphertext_t& operator=(const ciphertext_t &other)
    {
      *_c = *other._c;
      return *this;
    }
    
    ~ciphertext_t()
    {
      rns::commons::mem_manag::free_aligned<ciphertext_rns_t>(1, _c);
    }
  };

  struct keys_t
  {
    Pq *_keys;

    keys_t()
    {
      _keys = rns::commons::mem_manag::alloc_aligned<Pq, 32>(3);
    }

    keys_t(const keys_t &other)
    {
      _keys = rns::commons::mem_manag::alloc_aligned<Pq, 32>(3);
      _keys[0] = other._keys[0];
      _keys[1] = other._keys[1];
      _keys[2] = other._keys[2];
    }

    keys_t& operator=(const keys_t &other)
    {
      _keys[0] = other._keys[0];
      _keys[1] = other._keys[1];
      _keys[2] = other._keys[2];
      return *this;
    }
    
    ~keys_t()
    {
      rns::commons::mem_manag::free_aligned<Pq>(3, _keys);
    }
  };

  template<typename iterator_t>
  void import_it(plaintext_t &p, iterator_t it)
  {
    for (size_t i = 0; i < Pq::nmoduli; i++)
      {
	if (i == 0)
	  {
	    for (size_t j = 0; j < degree; j++)
	      {
		(*p._m)(i, j) = *it++;
	      }
	  }
	else
	  {
	    for (size_t j = 0; j < degree; j++)
	      {
		(*p._m)(i, j) = (*p._m)(i-1, j);
	      }
	  }
      }
  }

  template<typename iterator_t>
  void export_it(iterator_t it, plaintext_t &p)
  {
    typename P_t::iterator it1 = p._p->begin();
    for (size_t i = 0; i < degree; i++)
      {
	 *it++ = *it1++;
      }
  }

  void generate_keys(keys_t &k)
  {
    rns::commons::key_gen::key_gen<Pq, gamma_t_parameters_t>(k._keys);
  }

  void do_precomputations(keys_t &k)
  {
    cont_RNS->set_precomp(k._keys);
  }

  void encrypt(ciphertext_t &c, plaintext_t &p)
  {
    cont_RNS->encrypt(*c._c, *p._m);
  }

  void decrypt(plaintext_t &p, ciphertext_t &c)
  {
    cont_RNS->decrypt_not_ntt(*p._p, *c._c);
  }

  void multiply(ciphertext_t &c1, ciphertext_t &c2)
  {
    cont_RNS->mult(*c1._c, *c2._c);
  }

  void relinearise(ciphertext_t &c1)
  {
    cont_RNS->relinearisation(*c1._c);
  }

  void add(ciphertext_t &c1, ciphertext_t &c2)
  {
    c1._c->__cq[0] = c1._c->__cq[0] + c2._c->__cq[0];
    c1._c->__cq[1] = c1._c->__cq[1] + c2._c->__cq[1];
  }

  size_t get_degree()
  {
    return degree;
  }
};

#endif
