# Homomorphic Advanced RNS-based Encryption (HARE) library

This software was used to obtain the experimental results in "An HPR variant of the FV scheme: Computationally Cheaper, Asymptotically Faster". We are releasing it to allow for reproducible results. The library is limited, focusing on the arithmetic operations, while not providing any guarantees of security for real use-cases. We also do not provide any guarantees in terms of thread or reentrancy safety. What's more, most of the API functions modify their inputs, making them unviable for further calls, but allowing for a lower memory consumption and improved efficiency.

The file "hare.hpp" provides a common interface for all the supported cryptographic schemes (BEHZ, HPS, and HPR). More concretely, each of these schemes is associated with a class with the following typedefs/methods:

```
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
```

Thus, one might use C++ templates to write uniform code that targets all schemes. Nevertheless, to facilitate writing code for this interface, the macros BEGINMETHOD, ENDMETHOD and CALLMETHOD were defined that help with defining and calling functions that use homomorphic operations. For instance, the definition of a function returning void, named test_multiplication, accepting a parameter size_t num_mults could be defined as follows:

```
  BEGINMETHOD(void, //return type
	      test_multiplication, //name
	      size_t num_mults) //arguments
  {
    /* Function body */
  }
  ENDMETHOD()
```

Inside the function a ctx object is available, featuring the above-described functions. The types keys_t, ciphertext_t and plaintext_t are available both inside the functions and when defining the function arguments. The function might then be called with, for instance

```
  auto ctx = make_behz<primesize, bs_deg, nmoduli, bs_g, bs_t>(red);
  CALLMETHOD(ctx, test_multiplication, 10);
```

The file "factory.hpp", implicitly included by "hare.hpp" defines the functions `make_behz`, `make_hps` and `make_hpr_rns` that facilitate the creation of contexts.

## Modified NFLlib

This library relies on a heavily modified version of NFLlib. For more details check the "nfl" folder.

## Install

The code has the GMP and MPFR libraries as depedencies. The following commands might be used to generate makefiles on Linux

```
  cmake -DCMAKE_BUILD_TYPE=Release -H. -B_builds -G "Unix Makefiles"
```

and build the code

```
  cmake --build _builds/
```

Variations to build the code on other systems should be available by consulting the manpages of cmake and changing the "-G" flag accordingly.

## License

The code is licensed under GPLv3

## Contributors

The code has been developed by Julien Eynard (CEA), Paulo Martins (INESC-ID) and Vincent Zucca (University of Wollongong) (in alphabetical order).
