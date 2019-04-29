#ifndef __PARAMS_GAMMA_T_HPP__
#define __PARAMS_GAMMA_T_HPP__

template<typename main_type, size_t bs_g, size_t bs_t>
struct parameters_gamma_t
{
  typedef typename std::conditional< bs_g <= 8,
          uint_fast8_t,
          std::conditional< bs_g <= 16,
          uint_fast16_t,
          uint_fast32_t >>::type
          type_gamma;

  typedef typename std::conditional< bs_t + bs_g <= 16,
          uint_fast16_t,
          uint_fast32_t >::type
          type_gamma_t;

  typedef typename std::conditional< bs_t <= 8,
          uint_fast8_t,
          std::conditional< bs_t <= 16,
          uint_fast16_t,
          uint_fast32_t >>::type
          type_t;

  static constexpr main_type gamma = static_cast<main_type>(1)<<bs_g;
  static constexpr type_gamma gamma_half = static_cast<type_gamma>(gamma>>1);
  static constexpr type_gamma_t mask_gt = (static_cast<type_gamma_t>(1)<<
                                          (bs_t+bs_g))-1;

  // t
  static constexpr main_type t = static_cast<main_type>(1)<<bs_t;
  static constexpr type_t t_t = static_cast<type_t>(1)<<bs_t;
  static constexpr main_type mmask_t = t-1;
  static constexpr type_t mask_t = static_cast<type_t>(t-1);

  // mtilde

  // error distribution
  static constexpr double sigma_err = 8.0;
  static constexpr int B_err = 6*sigma_err + 0.5;
  static constexpr uint64_t B_key = 1;
  //static nfl::FastGaussianNoise<uint8_t, main_type, 2> g_prng (sigma_err, 128, 1UL << 14, (double)0, false);
  static nfl::PalisadeGaussianGenerator<main_type> g_prng;

  static constexpr size_t _bs_g = bs_g;
  static constexpr size_t _bs_t = bs_t;
};

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr main_type parameters_gamma_t<main_type, bs_g, bs_t>::gamma;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr typename parameters_gamma_t<main_type, bs_g, bs_t>::type_gamma
parameters_gamma_t<main_type, bs_g, bs_t>::gamma_half;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr typename parameters_gamma_t<main_type, bs_g, bs_t>::type_gamma_t
parameters_gamma_t<main_type, bs_g, bs_t>::mask_gt;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr main_type parameters_gamma_t<main_type, bs_g, bs_t>::t;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr typename parameters_gamma_t<main_type, bs_g, bs_t>::type_t
parameters_gamma_t<main_type, bs_g, bs_t>::t_t;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr main_type parameters_gamma_t<main_type, bs_g, bs_t>::mmask_t;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr typename parameters_gamma_t<main_type, bs_g, bs_t>::type_t
parameters_gamma_t<main_type, bs_g, bs_t>::mask_t;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr double parameters_gamma_t<main_type, bs_g, bs_t>::sigma_err;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr int parameters_gamma_t<main_type, bs_g, bs_t>::B_err;

template<typename main_type, size_t bs_g, size_t bs_t>
constexpr uint64_t parameters_gamma_t<main_type, bs_g, bs_t>::B_key;

template<typename main_type, size_t bs_g, size_t bs_t>
nfl::PalisadeGaussianGenerator<main_type>
parameters_gamma_t<main_type, bs_g, bs_t>::g_prng{parameters_gamma_t<main_type, bs_g, bs_t>::sigma_err};

#endif
