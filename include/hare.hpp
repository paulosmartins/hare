#ifndef __HARE_HPP__
#define __HARE_HPP__

/*

  Common interface for all schemes

  struct context_t
  {

  typedef (...) keys_t;
  typedef (...) ciphertext_t;
  typedef (...) plaintext_t;

  template<typename iterator_t>
  void import_it(plaintext_t&, iterator_t);
  template<typename iterator_t>
  void export_it(iterator_t, plaintext_t&);

  void generate_keys(keys_t&);
  void do_precomputations(keys_t&);

  void encrypt(ciphertext_t&, plaintext_t&);
  void decrypt(plaintext_t&, ciphertext_t&);

  void multiply(ciphertext_t&, ciphertext_t&);
  void relinearise(ciphertext_t&);
  void add(ciphertext_t&, ciphertext_t&);

  size_t get_degree();

  };

*/

/* Helpful macro to create functions based on the
   generic contexts */

#define BEGINMETHOD(ret, name, ...)				\
  template<typename context_t>					\
  struct name ## _t {						\
  context_t &ctx;						\
  typedef typename context_t::keys_t keys_t;			\
  typedef typename context_t::ciphertext_t ciphertext_t;	\
  typedef typename context_t::plaintext_t plaintext_t;		\
  name ## _t(context_t &_ctx) : ctx(_ctx) { }			\
  ret operator() (__VA_ARGS__)					\

#define ENDMETHOD()				\
  };
 
#define CALLMETHOD(ctx, name, ...)		\
  name ## _t<decltype(ctx)>{ctx}(__VA_ARGS__)

#include "rns/rns.hpp"
#include "hpr/hpr.hpp"

#include "behz.hpp"
#include "hps.hpp"
#include "hpr_rns.hpp"
#include "hpr_hps.hpp"
  
#include "params_gamma_t.hpp"
#include "factory.hpp"

#endif
