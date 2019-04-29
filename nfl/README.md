# Modified NFLlib

This is an heavily modified version of [NFLlib](https://github.com/quarkslab/NFLlib/). We have changed NFLlib so that it would support primes of different sizes; polynomials could be defined with respect to RNS bases defined as a range of primes starting from an arbitrary index of the arrays defined in the nfl::params_x::P structures, where x is a value between 46 and 62; and it would also produce numbers according to PALISADE Gaussian number generator.

## License

The original and the modified code are licensed under GPLv3

## Contributors

The original library was an extension/evolution of the NTTTools module from [XPIR](https://github.com/XPIR-team/XPIR) done by members of [CryptoExperts](https://www.cryptoexperts.com), [INP ENSEEIHT](http://www.enseeiht.com), [Quarkslab](http://www.quarkslab.com) (in alphabetical order).

Now, it has been modified by Julien Eynard (CEA), Paulo Martins (INESC-ID) and Vincent Zucca (University of Wollongong) (in alphabetical order).
