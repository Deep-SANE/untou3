// test_input.cpp - a test driver for untou3::UNtoU3 class.
//
// Program implements the U(N) to U(3) reduction for the input irrep [f]
// specified by the HO level n, N=(n+1)*(n+2)/2, and the elements of its
// representation R = [r_4,r_3,r_2,r_1,r_0] read from the standard input. For
// instance, for the input irrep [f] = [4,2,2,2,2,0], the user should provide
// the following input: 2 1 0 4 0 1.
//
// The program performs the U(N) to U(3) reduction, calculates the sum of the
// dimensions of resulting U(3) irrpes multiplied by their multiplicities,
// and prints it to the standard output. For instance, for the input irrep
// specified above, the output should read: U(3) irreps total dim = 405
//
// This sum should be equal to dim[f], which can be calculated analytically with
// the support of rational numbers. The program performs this calculcation as
// well if the Boost library is available and uses its Boost.Rational
// sublibrary. Availabitliy of Boost is indicated by users by definition of
// HAVE_BOOST preprocessor symbol. For the input irrep [f] specified above, the
// program should first print out: U(N) irrep dim = 405

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

#ifdef MEASURE_TIME
#include <chrono>
#endif

#ifdef HAVE_BOOST
#include <boost/rational.hpp>
#endif

#include <UNtoU3.h>

#ifdef HAVE_BOOST
// Impelments analytical formula for calculation of a dimension of a generic
// U(N) irrep [f]. [f] is specified by its labels passed as an array (generaly
// any indexed data structure) argument irrep.
template <typename T>
unsigned long dim(const T& irrep)
{
  const auto N = irrep.size();
  boost::rational<unsigned long> result{1};
  for (uint32_t l = 2; l <= N; l++)
    for (uint32_t k = 1; k <= l - 1; k++)
      result *= {irrep[k - 1] - irrep[l - 1] + l - k, l - k};

  assert(result.denominator() == 1);
  return result.numerator();
}
#endif

int main()
{
  // HO level
  unsigned long n;

  // specification of intput U(N) irrep
  unsigned int n4, n3, n2, n1, n0;

  std::cin >> n >> n4 >> n3 >> n2 >> n1 >> n0;

  if (n4 + n3 + n2 + n1 + n0 != (n + 1) * (n + 2) / 2)
    throw std::invalid_argument("Arguments mismatch!");

#ifdef HAVE_BOOST
  // analytical calculation of dim([f])
  std::vector<unsigned long> f(n4, 4);
  std::fill_n(std::back_inserter(f), n3, 3);
  std::fill_n(std::back_inserter(f), n2, 2);
  std::fill_n(std::back_inserter(f), n1, 1);
  std::fill_n(std::back_inserter(f), n0, 0);
  std::cout << "U(N) irrep dim = " << dim(f) << std::endl;
#endif

  // maximum U(N) irrep label = 4
  untou3::UNtoU3<4> gen;

  // set HO level
  gen.set_shell(n);

#ifdef MEASURE_TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  // U(N) to U(3) reduction
  gen.reduce({static_cast<int16_t>(n4), static_cast<int16_t>(n3),
              static_cast<int16_t>(n2), static_cast<int16_t>(n1),
              static_cast<int16_t>(n0)});

#ifdef MEASURE_TIME
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "U3 weights generation time: "
            << std::chrono::duration<double>(end - start).count() << " [s]"
            << std::endl;
#endif

  // resulting U(3) irreps
  auto res = gen.irreps();

  long total_dim = 0;
  for (const auto& [i, m] : res)
    total_dim += static_cast<long>(m) * i.dim();

  std::cout << "U(3) irreps total dim = " << total_dim << std::endl;
}
