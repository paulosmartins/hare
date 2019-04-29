#ifndef __FACTORY_HPP__
#define __FACTORY_HPP__

template<size_t i>
struct get_nfl_params_t { };

#define SPECIALISE_GET_NFL_PARAMS(i)		\
  template<>					\
  struct get_nfl_params_t<i>			\
  {						\
    typedef nfl::params_ ## i type;		\
  }

SPECIALISE_GET_NFL_PARAMS(46);
SPECIALISE_GET_NFL_PARAMS(47);
SPECIALISE_GET_NFL_PARAMS(48);
SPECIALISE_GET_NFL_PARAMS(49);
SPECIALISE_GET_NFL_PARAMS(50);
SPECIALISE_GET_NFL_PARAMS(51);
SPECIALISE_GET_NFL_PARAMS(52);
SPECIALISE_GET_NFL_PARAMS(53);
SPECIALISE_GET_NFL_PARAMS(54);
SPECIALISE_GET_NFL_PARAMS(55);
SPECIALISE_GET_NFL_PARAMS(56);
SPECIALISE_GET_NFL_PARAMS(57);
SPECIALISE_GET_NFL_PARAMS(58);
SPECIALISE_GET_NFL_PARAMS(59);
SPECIALISE_GET_NFL_PARAMS(60);
SPECIALISE_GET_NFL_PARAMS(61);

#undef SPECIALISE_GET_NFL_PARAMS

template<size_t bs_deg>
struct get_mt
{
  typedef typename std::conditional<
    bs_deg <= 14,
	      std::integral_constant<size_t, 8>,
	      std::integral_constant<size_t, 16>>::type type;
  static constexpr size_t value = type::value;
};

template<size_t primesize, size_t bs_g, size_t bs_t>
struct get_parameters_gamma_t
{
  typedef typename get_nfl_params_t<primesize>::type nfl_params_t;
  typedef typename nfl_params_t::value_type value_type;
  typedef parameters_gamma_t<value_type, bs_g, bs_t> type;
};

template<size_t primesize, size_t bs_deg, size_t nmoduli, size_t bs_g, size_t bs_t>
auto make_behz(bool red = true)
  -> decltype (context_behz_t<
	       typename get_nfl_params_t<primesize>::type,
	       (1<<bs_deg),
	       nmoduli,
	       typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
	       get_mt<bs_deg>::value> (red))
{
  return context_behz_t<
    typename get_nfl_params_t<primesize>::type,
    (1<<bs_deg),
      nmoduli,
      typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
      get_mt<bs_deg>::value> (red);
}

template<size_t primesize, size_t bs_deg, size_t nmoduli, size_t bs_g, size_t bs_t>
auto make_hps()
  -> decltype (context_hps_t<
	       typename get_nfl_params_t<primesize>::type,
	       (1<<bs_deg),
	       nmoduli,
	       typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
	       get_mt<bs_deg>::value> ())
{
  return context_hps_t<
    typename get_nfl_params_t<primesize>::type,
    (1<<bs_deg),
      nmoduli,
      typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
      get_mt<bs_deg>::value> ();
}

template<size_t primesize, size_t bs_deg, size_t nmoduli, size_t ndigits,
	 size_t bs_g, size_t bs_t>
auto make_hpr_rns(bool red, bool redhps)
  -> decltype (context_hpr_rns_t<
	       typename get_nfl_params_t<primesize>::type,
	       (1<<bs_deg),
	       nmoduli,
	       ndigits,
	       typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
	       get_mt<bs_deg>::value> (red, redhps))
{
  return context_hpr_rns_t<
    typename get_nfl_params_t<primesize>::type,
    (1<<bs_deg),
      nmoduli,
      ndigits,
      typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
      get_mt<bs_deg>::value> (red, redhps);
}

template<size_t primesize, size_t bs_deg, size_t nmoduli, size_t ndigits,
	 size_t bs_g, size_t bs_t>
auto make_hpr_hps()
  -> decltype (context_hpr_hps_t<
	       typename get_nfl_params_t<primesize>::type,
	       (1<<bs_deg),
	       nmoduli,
	       ndigits,
	       typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
	       get_mt<bs_deg>::value> ())
{
  return context_hpr_hps_t<
    typename get_nfl_params_t<primesize>::type,
    (1<<bs_deg),
      nmoduli,
      ndigits,
      typename get_parameters_gamma_t<primesize, bs_g, bs_t>::type,
      get_mt<bs_deg>::value> ();
}



#endif
