namespace nfl {

  namespace ops {

    // Fused Multiplications-Additions
    // FMA with division for modular reduction (expensive)
    template<class type, class tag> struct muladd;

    template<class T>
    struct muladd<T, simd::serial> {
      using poly_type = T;
      using simd_mode = simd::serial;
      using value_type = typename poly_type::value_type;
      using greater_value_type = typename poly_type::greater_value_type;
      static_assert(std::is_same<value_type, uint64_t>(), "this code has only been tested for uint64_t");
      value_type operator()(value_type rop, value_type x, value_type y, size_t cm) const
      {
	auto const p = poly_type::base.params.P[cm];
	auto const shift = poly_type::base.params.kModulusRepresentationBitsize;
	auto const s0 = shift-poly_type::base.params.kModulusBitsize;;
	ASSERT_STRICTMOD((x<p) && (y<p));
	greater_value_type res = (greater_value_type)x * y;
	ASSERT_STRICTMOD(res < (std::numeric_limits<greater_value_type>::max()-rop));
	greater_value_type q = ((greater_value_type)poly_type::base.params.Pn[cm] * (res >> shift)) + (res<<s0) ;
	value_type r  = res - (q>>shift) * p;
	if (r >= p) { r -= p; }
	r += rop;
	if (r >= p) { r -= p; }
	ASSERT_STRICTMOD(r == ((((greater_value_type)(x) * y) % p) + rop) % p);
	return r;
      }
    };

    // FMA using Shoup's precomputed approach again
    // Works only if yshoup = ((uint128_t) y << 64) / p (for T=uint64_t)
    // has been pre-computed AND p is small enough (2 bits less than limb).
    // OUTPUT: x * y mod p
    template<class T, class tag> struct muladd_shoup;

    template<class T>
    struct muladd_shoup<T, simd::serial> {
      using poly_type = T;
      using simd_mode = simd::serial;
      using value_type = typename poly_type::value_type;
      using greater_value_type = typename poly_type::greater_value_type;
      value_type operator()(value_type rop, value_type x, value_type y, value_type yprime, size_t cm) const
      {
	auto const p = poly_type::base.params.P[cm];
	ASSERT_STRICTMOD(rop<p && x<p && y<p);
	value_type q = ((greater_value_type) x * yprime) >> poly_type::base.params.kModulusRepresentationBitsize;
#ifdef ENFORCE_STRICTMOD
	value_type res;
	res = x * y - q * p;
	res -= ((res>=p) ? p : 0);
	rop += res;
#else
	rop += x * y - q * p;
#endif
	ASSERT_STRICTMOD(rop - ((rop>=p) ? p : 0)<p);
	return rop - ((rop>=p) ? p : 0);
      }
    };



  } // ops

} // nfl
