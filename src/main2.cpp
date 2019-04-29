#include "hare.hpp"
#include <random>
#include <chrono>
#include <map>

std::vector<int> multiply(const std::vector<int> &a,
			  const std::vector<int> &b)
{
  size_t degree = a.size();
  std::vector<int> c(degree, 0);

  for (size_t i = 0; i < degree; i++)
    {
      for (size_t j = 0; j < degree; j++)
	{
	  if (j <= i)
	    c[i] ^= a[j] & b[i-j];
	  else
	    c[i] ^= a[j] & b[degree+i-j];
	}
    }

  return c;
}

/*BEGINMETHOD makes ctx available, featuring the
  functions described in hare.hpp, as well as the types
  keys_t, ciphertext_t and plaintext_t */
BEGINMETHOD(void, //return type
	    test_multiplication, //name
	    size_t num_mults) //arguments
{
  std::default_random_engine gen;
  std::uniform_int_distribution<int> dist(0,1);

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed;
  double enc_time = 0.0;
  double mult_time = 0.0;
  double relin_time = 0.0;
  double dec_time = 0.0;
  
  for (size_t i = 0; i < num_mults; i++)
    {
      std::vector<int> msg1(ctx.get_degree());
      std::vector<int> msg2(ctx.get_degree());

      for (size_t j = 0; j < ctx.get_degree(); j++)
	{
	  msg1[j] = dist(gen);
	  msg2[j] = dist(gen);
	}

      plaintext_t pt1;
      ctx.import_it(pt1, std::begin(msg1));
      plaintext_t pt2;
      ctx.import_it(pt2, std::begin(msg2));

      ciphertext_t ct1;
      ctx.encrypt(ct1, pt1);
      ciphertext_t ct2;
      start = std::chrono::high_resolution_clock::now();
      ctx.encrypt(ct2, pt2);
      end = std::chrono::high_resolution_clock::now();
      elapsed = end - start;

      enc_time += elapsed.count();

      start = std::chrono::high_resolution_clock::now();
      ctx.multiply(ct1, ct2);
      end = std::chrono::high_resolution_clock::now();
      elapsed = end - start;

      mult_time += elapsed.count();

      start = std::chrono::high_resolution_clock::now();
      ctx.relinearise(ct1);
      end = std::chrono::high_resolution_clock::now();
      elapsed = end - start;

      relin_time += elapsed.count();

      std::vector<int> msg3_expected = multiply(msg1, msg2);
      std::vector<int> msg3(ctx.get_degree());

      plaintext_t pt3;
      
      start = std::chrono::high_resolution_clock::now();
      ctx.decrypt(pt3, ct1);
      end = std::chrono::high_resolution_clock::now();
      elapsed = end - start;

      dec_time += elapsed.count();
      
      ctx.export_it(std::begin(msg3), pt3);

      for (size_t j = 0; j < ctx.get_degree(); j++)
	{
	  assert(msg3[j] == msg3_expected[j]);
	}
    }

  enc_time /= num_mults;
  mult_time /= num_mults;
  relin_time /= num_mults;
  dec_time /= num_mults;

  std::cout << "enc_time = " << enc_time << "s\n";
  std::cout << "mult_time = " << mult_time << "s\n";
  std::cout << "relin_time = " << relin_time << "s\n";
  std::cout << "dec_time = " << dec_time << "s\n";
}
ENDMETHOD()

BEGINMETHOD(void,
	    test_depth,
	    size_t num_tries)
{
  std::vector<size_t> depth(num_tries);
  std::default_random_engine gen;
  std::uniform_int_distribution<int> dist(0,1);
  
  for (size_t i = 0; i < num_tries; i++)
    {
      size_t num_mults = 0;
      std::vector<int> msg1(ctx.get_degree());

      for (size_t j = 0; j < ctx.get_degree(); j++)
	{
	  msg1[j] = dist(gen);
	}

      plaintext_t pt1;
      ctx.import_it(pt1, std::begin(msg1));
      ciphertext_t ct1;
      ctx.encrypt(ct1, pt1);

      bool failed_decrypt = false;

      while (!failed_decrypt)
	{
	  ciphertext_t tmp = ct1;
	  ctx.multiply(ct1, tmp);
	  ctx.relinearise(ct1);

	  tmp = ct1;
	  std::vector<int> msg2(ctx.get_degree());
	  plaintext_t pt2;
	  ctx.decrypt(pt2, tmp);
	  ctx.export_it(std::begin(msg2), pt2);

	  msg1 = multiply(msg1, msg1);

	  for (size_t j = 0; j < ctx.get_degree(); j++)
	    {
	      if (msg1[j] != msg2[j])
		failed_decrypt = true;
	    }

	  if (!failed_decrypt)
	    num_mults++;
	}

      depth[i] = num_mults;
    }

  double avg_depth = 0.0;
  for (size_t i = 0; i < num_tries; i++)
    {
      avg_depth += depth[i];
    }
  avg_depth /= num_tries;

  double std_dev = 0.0;
  for (size_t i = 0; i < num_tries; i++)
    {
      std_dev += (depth[i] - avg_depth) * (depth[i] - avg_depth);
    }
  std_dev /= (num_tries - 1); //unbiased estimate
  std_dev = sqrt(std_dev);

  size_t min_depth = 1000000;
  for (size_t i = 0; i < num_tries; i++)
    {
      min_depth = std::min(min_depth, depth[i]);
    }

  size_t max_depth = 0;
  for (size_t i = 0; i < num_tries; i++)
    {
      max_depth = std::max(max_depth, depth[i]);
    }

  std::cout << "avg_depth = " << avg_depth << "\n";
  std::cout << "std_dev = " << std_dev << "\n";
  std::cout << "min_depth = " << min_depth << "\n";
  std::cout << "max_depth = " << max_depth << "\n";
}
ENDMETHOD()

BEGINMETHOD(void,
	    initialise)
{
  keys_t keys;
  ctx.generate_keys(keys);
  ctx.do_precomputations(keys);
}
ENDMETHOD()

void test_behz_13()
{
  static constexpr size_t primesize = 55;
  static constexpr size_t bs_deg = 13;
  static constexpr size_t nmoduli = 5;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;
  static constexpr bool red = true;

  {
    std::cout << "####################\n";
  
    std::cout << "BEHZ\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n"
	      << "red=\t" << red << "\n";
    auto ctx = make_behz<primesize, bs_deg, nmoduli, bs_g, bs_t>(red);
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_behz_14()
{
  static constexpr size_t primesize = 61;
  static constexpr size_t bs_deg = 14;
  static constexpr size_t nmoduli = 9;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;
  static constexpr bool red = true;

  {
    std::cout << "####################\n";
  
    std::cout << "BEHZ\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n"
	      << "red=\t" << red << "\n";
    auto ctx = make_behz<primesize, bs_deg, nmoduli, bs_g, bs_t>(red);
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_behz_15()
{
  static constexpr size_t primesize = 61;
  static constexpr size_t bs_deg = 15;
  static constexpr size_t nmoduli = 18;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;
  static constexpr bool red = true;

  {
    std::cout << "####################\n";
  
    std::cout << "BEHZ\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n"
	      << "red=\t" << red << "\n";
    auto ctx = make_behz<primesize, bs_deg, nmoduli, bs_g, bs_t>(red);
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_behz_16()
{
  static constexpr size_t primesize = 61;
  static constexpr size_t bs_deg = 16;
  static constexpr size_t nmoduli = 36;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;
  static constexpr bool red = true;

  {
    std::cout << "####################\n";
  
    std::cout << "BEHZ\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n"
	      << "red=\t" << red << "\n";
    auto ctx = make_behz<primesize, bs_deg, nmoduli, bs_g, bs_t>(red);
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_hps_13()
{
  static constexpr size_t primesize = 55;
  static constexpr size_t bs_deg = 13;
  static constexpr size_t nmoduli = 5;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;

  {
    std::cout << "####################\n";
  
    std::cout << "HPS\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n";
    auto ctx = make_hps<primesize, bs_deg, nmoduli, bs_g, bs_t>();
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_hps_14()
{
  static constexpr size_t primesize = 61;
  static constexpr size_t bs_deg = 14;
  static constexpr size_t nmoduli = 9;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;

  {
    std::cout << "####################\n";
  
    std::cout << "HPS\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n";
    auto ctx = make_hps<primesize, bs_deg, nmoduli, bs_g, bs_t>();
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_hps_15()
{
  static constexpr size_t primesize = 61;
  static constexpr size_t bs_deg = 15;
  static constexpr size_t nmoduli = 18;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;

  {
    std::cout << "####################\n";
  
    std::cout << "HPS\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n";
    auto ctx = make_hps<primesize, bs_deg, nmoduli, bs_g, bs_t>();
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_hps_16()
{
  static constexpr size_t primesize = 61;
  static constexpr size_t bs_deg = 16;
  static constexpr size_t nmoduli = 36;
  static constexpr size_t bs_g = 8;
  static constexpr size_t bs_t = 1;
  static constexpr bool red = true;

  {
    std::cout << "####################\n";
  
    std::cout << "HPS\n"
	      << "primesize=\t" << primesize << "\n"
	      << "bs_deg=\t" << bs_deg << "\n"
	      << "nmoduli=\t" << nmoduli << "\n"
	      << "bs_g=\t" << bs_g << "\n"
	      << "bs_t=\t" << bs_t << "\n";
    auto ctx = make_hps<primesize, bs_deg, nmoduli, bs_g, bs_t>();
    CALLMETHOD(ctx, initialise);
    CALLMETHOD(ctx, test_multiplication, 100);
    CALLMETHOD(ctx, test_depth, 30);
  }
}

void test_hpr_13()
{
  static constexpr size_t hpr_primesize = 55;
  static constexpr size_t hpr_bs_deg = 13;
  static constexpr size_t hpr_nmoduli = 1;
  static constexpr size_t hpr_ndigits = 5;
  static constexpr size_t hpr_bs_g = 8;
  static constexpr size_t hpr_bs_t = 1;
  static constexpr bool hpr_red = true;
  static constexpr bool hpr_redhps = true;
  
  {
    std::cout << "####################\n";
  
    std::cout << "HPR-RNS\n"
	      << "primesize=\t" << hpr_primesize << "\n"
	      << "bs_deg=\t" << hpr_bs_deg << "\n"
	      << "nmoduli=\t" << hpr_nmoduli << "\n"
	      << "ndigits=\t" << hpr_ndigits << "\n"
	      << "bs_g=\t" << hpr_bs_g << "\n"
	      << "bs_t=\t" << hpr_bs_t << "\n"
	      << "red=\t" << hpr_red << "\n"
	      << "redhps=\t" << hpr_redhps << "\n";
  
    auto ctx3 = make_hpr_rns<hpr_primesize,
			     hpr_bs_deg,
			     hpr_nmoduli,
			     hpr_ndigits,
			     hpr_bs_g,
			     hpr_bs_t>(hpr_red, hpr_redhps);
    CALLMETHOD(ctx3, initialise);
    CALLMETHOD(ctx3, test_multiplication, 100);
    CALLMETHOD(ctx3, test_depth, 30);
  }
}

void test_hpr_14()
{
  static constexpr size_t hpr_primesize = 61;
  static constexpr size_t hpr_bs_deg = 14;
  static constexpr size_t hpr_nmoduli = 3;
  static constexpr size_t hpr_ndigits = 3;
  static constexpr size_t hpr_bs_g = 8;
  static constexpr size_t hpr_bs_t = 1;
  static constexpr bool hpr_red = true;
  static constexpr bool hpr_redhps = true;
  
  {
    std::cout << "####################\n";
  
    std::cout << "HPR-RNS\n"
	      << "primesize=\t" << hpr_primesize << "\n"
	      << "bs_deg=\t" << hpr_bs_deg << "\n"
	      << "nmoduli=\t" << hpr_nmoduli << "\n"
	      << "ndigits=\t" << hpr_ndigits << "\n"
	      << "bs_g=\t" << hpr_bs_g << "\n"
	      << "bs_t=\t" << hpr_bs_t << "\n"
	      << "red=\t" << hpr_red << "\n"
	      << "redhps=\t" << hpr_redhps << "\n";
  
    auto ctx3 = make_hpr_rns<hpr_primesize,
			     hpr_bs_deg,
			     hpr_nmoduli,
			     hpr_ndigits,
			     hpr_bs_g,
			     hpr_bs_t>(hpr_red, hpr_redhps);
    CALLMETHOD(ctx3, initialise);
    CALLMETHOD(ctx3, test_multiplication, 100);
    CALLMETHOD(ctx3, test_depth, 30);
  }
}

void test_hpr_15()
{
  static constexpr size_t hpr_primesize = 61;
  static constexpr size_t hpr_bs_deg = 15;
  static constexpr size_t hpr_nmoduli = 3;
  static constexpr size_t hpr_ndigits = 6;
  static constexpr size_t hpr_bs_g = 8;
  static constexpr size_t hpr_bs_t = 1;
  static constexpr bool hpr_red = true;
  static constexpr bool hpr_redhps = true;
  
  {
    std::cout << "####################\n";
  
    std::cout << "HPR-RNS\n"
	      << "primesize=\t" << hpr_primesize << "\n"
	      << "bs_deg=\t" << hpr_bs_deg << "\n"
	      << "nmoduli=\t" << hpr_nmoduli << "\n"
	      << "ndigits=\t" << hpr_ndigits << "\n"
	      << "bs_g=\t" << hpr_bs_g << "\n"
	      << "bs_t=\t" << hpr_bs_t << "\n"
	      << "red=\t" << hpr_red << "\n"
	      << "redhps=\t" << hpr_redhps << "\n";
  
    auto ctx3 = make_hpr_rns<hpr_primesize,
			     hpr_bs_deg,
			     hpr_nmoduli,
			     hpr_ndigits,
			     hpr_bs_g,
			     hpr_bs_t>(hpr_red, hpr_redhps);
    CALLMETHOD(ctx3, initialise);
    CALLMETHOD(ctx3, test_multiplication, 100);
    CALLMETHOD(ctx3, test_depth, 30);
  }
}

void test_hpr_16()
{
  static constexpr size_t hpr_primesize = 61;
  static constexpr size_t hpr_bs_deg = 16;
  static constexpr size_t hpr_nmoduli = 3;
  static constexpr size_t hpr_ndigits = 12;
  static constexpr size_t hpr_bs_g = 8;
  static constexpr size_t hpr_bs_t = 1;
  static constexpr bool hpr_red = true;
  static constexpr bool hpr_redhps = true;
  
  {
    std::cout << "####################\n";
  
    std::cout << "HPR-RNS\n"
	      << "primesize=\t" << hpr_primesize << "\n"
	      << "bs_deg=\t" << hpr_bs_deg << "\n"
	      << "nmoduli=\t" << hpr_nmoduli << "\n"
	      << "ndigits=\t" << hpr_ndigits << "\n"
	      << "bs_g=\t" << hpr_bs_g << "\n"
	      << "bs_t=\t" << hpr_bs_t << "\n"
	      << "red=\t" << hpr_red << "\n"
	      << "redhps=\t" << hpr_redhps << "\n";
  
    auto ctx3 = make_hpr_rns<hpr_primesize,
			     hpr_bs_deg,
			     hpr_nmoduli,
			     hpr_ndigits,
			     hpr_bs_g,
			     hpr_bs_t>(hpr_red, hpr_redhps);
    CALLMETHOD(ctx3, initialise);
    CALLMETHOD(ctx3, test_multiplication, 100);
    CALLMETHOD(ctx3, test_depth, 30);
  }
}

const std::map<std::string, std::function<void()>> functests =
  {{"test_behz_13", test_behz_13},
   {"test_behz_14", test_behz_14},
   {"test_behz_15", test_behz_15},
   {"test_behz_16", test_behz_16},
   {"test_hps_13", test_hps_13},
   {"test_hps_14", test_hps_14},
   {"test_hps_15", test_hps_15},
   {"test_hps_16", test_hps_16},
   {"test_hpr_13", test_hpr_13},
   {"test_hpr_14", test_hpr_14},
   {"test_hpr_15", test_hpr_15},
   {"test_hpr_16", test_hpr_16}};


int main(int argc, char *argv[])
{
  if (argc != 2 || functests.find(argv[1]) == std::end(functests))
    {
      std::cerr << "usage: ./main2 funcname\n";
      std::cerr << "\twith funcname = test_(behz/hps/hpr)_(13/14/15/16)\n";
      return 1;
    }

  functests.at(argv[1])();
  
  return 0;
}
