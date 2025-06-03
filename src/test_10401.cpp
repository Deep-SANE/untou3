// test_10401.cpp - a test driver for untou3::UNtoU3 class.
//
// Program implements the U(N) to U(3) reduction where n=2 (N=6) for the input
// irrep [4,2,2,2,2,0], which is represented by R = [1,0,4,0,1]. After the
// reduction, the program iterates over resulting U(3) irrps, prints their
// multiplicities, and verify the total dimension.

#include <iostream>

#include <UNtoU3.h>

int main()
{
  // maximum U(N) irrep label = 4
  untou3::UNtoU3<4> gen;

  // n = 2 - a given HO level
  gen.set_shell(2);

  // U(N) to U(3) reduction
  gen.reduce({1, 0, 4, 0, 1});

  // resulting U(3) irreps
  auto res = gen.irreps();

  std::cout << "Resulting U(3) irrpes:" << std::endl;

  int total_dim = 0;
  for (const auto& [i, m] : res)
  {
    std::cout << i << " : mult = " << m << ", dim = " << i.dim() << std::endl;
    total_dim += m * i.dim();
  }

  std::cout << std::endl;
  std::cout << "Multiplicity of U(3) irrep [10,8,6] = "
            << gen.multiplicity({10, 8, 6}) << std::endl;

  std::cout << std::endl;
  std::cout << "Dimension of input U(N=6) irrep: 405" << std::endl;
  std::cout << "Dimension of reduced U(3) irreps: " << total_dim << std::endl;
}
