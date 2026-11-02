// test_custom_irreps.cpp - a test driver for untou3::UNtoU3 class.
//
// Program demonstrates the use of untou3::UNtoU3 with user-defined types
// representing U(3) and U(4) irreps, which are plugged in by specializing
// untou3::IrrepTraits. It performs the U(N) to U(3) reduction where n=2 (N=6)
// for the input irrep [4,2,2,2,2,0], which is represented by R = [1,0,4,0,1],
// and the conjugated reduction for the U(4) irrep [5,5,1,1], with both the
// user-defined and the default irrep types, and verifies that the results are
// the same.

#include <array>
#include <iostream>
#include <vector>

#include <UNtoU3.h>

// user-defined U(3) irrep type with labels accessed by member functions:
class MyU3
{
public:
  MyU3(int a, int b, int c) : l_{a, b, c} { }

  int label(int k) const { return l_[k]; }

private:
  std::array<int, 3> l_;
};

// user-defined U(4) irrep type with labels stored in an array:
struct MyU4
{
  std::array<int, 4> f;
};

template <>
struct untou3::IrrepTraits<MyU3>
{
  static std::array<int, 3> labels(const MyU3& i)
  {
    return {i.label(0), i.label(1), i.label(2)};
  }

  static MyU3 make(int e1, int e2, int e3) { return {e1, e2, e3}; }
};

template <>
struct untou3::IrrepTraits<MyU4>
{
  static std::array<int, 4> labels(const MyU4& i) { return i.f; }
};

// resulting U(3) irreps as labels and multiplicities:
template <typename Gen>
std::vector<std::array<int, 4>> results(const Gen& gen)
{
  using U3 = typename Gen::u3_irrep;

  std::vector<std::array<int, 4>> res;
  for (const auto& [i, m] : gen.irreps())
  {
    auto [e1, e2, e3] = untou3::IrrepTraits<U3>::labels(i);
    res.push_back({e1, e2, e3, m});
  }
  return res;
}

int main()
{
  // default irrep types:
  untou3::UNtoU3<4> gen;
  // user-defined irrep types:
  untou3::UNtoU3<4, MyU3, MyU4> my_gen;

  // n = 2 - a given HO level
  gen.set_shell(2);
  my_gen.set_shell(2);

  int errors = 0;

  // U(N) to U(3) reduction:
  gen.reduce({1, 0, 4, 0, 1});
  my_gen.reduce({1, 0, 4, 0, 1});
  if (results(gen) != results(my_gen))
    errors++;

  if (my_gen.multiplicity({10, 8, 6}) != 3)
    errors++;
  if (my_gen.multiplicity({20, 0, 0}) != 0)
    errors++;

  // U(4) to U(3) conjugated reduction:
  auto reduced = results(gen);
  gen.reduceU4({5, 5, 1, 1});
  my_gen.reduceU4(MyU4{{5, 5, 1, 1}});
  if ((results(gen) != results(my_gen)) || (results(gen) != reduced))
    errors++;

  std::cout << "Resulting U(3) irreps (user-defined types):" << std::endl;
  for (const auto& [i, m] : my_gen.irreps())
    std::cout << "[" << i.label(0) << "," << i.label(1) << "," << i.label(2)
              << "] : mult = " << m << std::endl;

  std::cout << std::endl;
  std::cout << (errors == 0 ? "OK" : "FAILED") << std::endl;

  return errors == 0 ? 0 : 1;
}
