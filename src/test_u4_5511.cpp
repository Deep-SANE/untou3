// test_u4_5511.cpp - a test driver for untou3::UNtoU3 class.
//
// Program implements the conjugated U(N) to U(3) reduction where n=2 (N=6) for
// the input U(4) irrep [5,5,1,1], for which the conjugated U(6) irrep is
// [4,2,2,2,2,0] represented by R = [1,0,4,0,1]. After the reduction, the
// program iterates over resulting U(3) irrps and prints their multiplicities.

#include <iostream>

#include <UNtoU3.h>

int main()
{
  // maximum U(N) irrep label = 4
  untou3::UNtoU3<4> gen;

  // n = 2 - a given HO level
  gen.set_shell(2);

  // U(4) to U(3) conjugated reduction
  untou3::U4Irrep i = {5, 5, 1, 1};
  gen.reduceU4(i);

  // resulting U(3) irreps
  auto res = gen.irreps();

  std::cout << "Resulting U(3) irrpes:" << std::endl;

  int total_dim = 0;
  for (const auto& [i, m] : res)
  {
    std::cout << i << " : mult = " << m << ", dim = " << i.dim() << std::endl;
    total_dim += m * i.dim();
  }
}
