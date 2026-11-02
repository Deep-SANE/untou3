# UNtoU3

UNtoU3 is a header-only C++20 library for the reduction of representations of the
unitary group U(N) into representations of U(3), the symmetry group of the
three-dimensional harmonic oscillator (HO). For a given U(N) irreducible
representation (irrep), it enumerates all U(3) irreps contained in it together
with their multiplicities (*level dimensionalities*).

The library implements the **UNtoU3-4** algorithm, which is designed for U(N)
irreps whose labels go up to Λ = 4. This is the case of isospin-adapted nuclear
structure calculations, where protons and neutrons are treated as the same kind
of particle and each HO eigenstate can be occupied by up to 4 nucleons with
different combinations of spin and isospin. The U(N) → U(3) reduction is a key
ingredient of the symmetry-adapted no-core shell model (SA-NCSM), an *ab initio*
approach to the structure of light and medium-mass nuclei.

The algorithm and its performance are described in the
[PPAM 2024 paper](#algorithm-and-publication).

## Background

Consider the η-th HO level, η = 0, 1, 2, …, with its $N = (\eta+1)(\eta+2)/2$
degenerate eigenstates. A U(N) irrep associated with a distribution of nucleons
among these eigenstates is labeled by

$$[f] = [f_1, \ldots, f_N], \qquad f_i \in \lbrace 0,1,2,3,4 \rbrace, \qquad f_i \ge f_{i+1}.$$

Since the labels can only take five values, UNtoU3 represents [f] by the counts
of each value:

$$R = [r_4, r_3, r_2, r_1, r_0],$$

where $r_k$ is the number of labels equal to $k$ and $r_4 + \cdots + r_0 = N$.
For instance, the U(6) irrep [4,2,2,2,2,0] (η = 2) is represented by
R = [1,0,4,0,1].

The result of the reduction is the set of U(3) irreps $[e] = [e_1, e_2, e_3]$,
$e_1 \ge e_2 \ge e_3$, contained in [f], each with its level dimensionality
$D_l[e]$, i.e., the number of times [e] occurs in [f]. The result can be checked
via

$$\sum_{[e]} D_l[e] \cdot \dim[e] = \dim[f].$$

## Requirements

- A C++20 compiler (the library uses `std::span`); tested with GCC 13.
- No other dependencies. [Boost](https://www.boost.org/) (Boost.Rational) is
  optionally used by the `test_input` program only, to compute dim[f]
  analytically.

## Usage

Being header-only, the library needs no build step: add `include/` to the include
path and include `UNtoU3.h`. The following program reduces the U(6) irrep
[4,2,2,2,2,0]:

```cpp
#include <iostream>

#include <UNtoU3.h>

int main()
{
  // maximum U(N) label: 4
  untou3::UNtoU3<4> gen;

  // HO level eta = 2, i.e., N = (eta+1)(eta+2)/2 = 6
  gen.set_shell(2);

  // input U(6) irrep [4,2,2,2,2,0] given by R = [r4,r3,r2,r1,r0] = [1,0,4,0,1]
  gen.reduce({1, 0, 4, 0, 1});

  long total_dim = 0;

  // resulting U(3) irreps and their level dimensionalities
  for (const auto& [i, m] : gen.irreps())
  {
    std::cout << i << " : mult = " << m << ", dim = " << i.dim() << "\n";
    total_dim += static_cast<long>(m) * i.dim();
  }

  std::cout << "total dim = " << total_dim << "\n";

  // level dimensionality of a single U(3) irrep
  std::cout << "mult [10,8,6] = " << gen.multiplicity({10, 8, 6}) << "\n";
}
```

Compiled with `g++ -std=c++20 -O2 -I include example.cpp`, it prints:

```text
[8,8,8] : mult = 1, dim = 1
[9,8,7] : mult = 1, dim = 8
[10,8,6] : mult = 3, dim = 27
[10,9,5] : mult = 1, dim = 35
[10,10,4] : mult = 1, dim = 28
[11,7,6] : mult = 1, dim = 35
[11,8,5] : mult = 1, dim = 64
[12,6,6] : mult = 1, dim = 28
[12,8,4] : mult = 1, dim = 125
total dim = 405
mult [10,8,6] = 3
```

which agrees with dim[4,2,2,2,2,0] = 405.

The conjugated reduction takes a U(4) irrep and reduces the conjugate U(N)
irrep (with the transposed Young diagram). For the U(4) irrep [5,5,1,1] and
N = 6, the conjugate U(6) irrep is [4,2,2,2,2,0], so the following gives the
same U(3) irreps as above (see `src/test_u4_5511.cpp`):

```cpp
untou3::UNtoU3<4> gen;
gen.set_shell(2);
gen.reduceU4({5, 5, 1, 1});
auto res = gen.irreps();
```

### API

The reduction is provided by the class template
`untou3::UNtoU3<L, U3 = U3Irrep, U4 = U4Irrep>`, where `L` is the maximum U(N)
label and `U3` and `U4` are the types representing U(3) and U(4) irreps (see
[Using other irrep types](#using-other-irrep-types)):

| Member | Description |
|---|---|
| `void set_shell(int n)` | Selects the `n`-th HO level, N = (n+1)(n+2)/2. Must be called first. |
| `int shell() const` | The HO level selected by `set_shell()` (-1 if none). |
| `void reduce(const GPR_rep& R)` | Reduces the U(N) irrep given by `R = {r_L, ..., r_0}` (`GPR_rep` is `std::array<int16_t, L + 1>`); the counts must sum to N. Replaces the results of any previous reduction. |
| `void reduceU4(const U4& i)` | Reduces the U(N) irrep conjugate to the U(4) irrep `i`. Available for `L = 4`. |
| `void reduceU4(const U4& i, int N)` | The same, with `N` given explicitly; it must match the selected HO level. |
| `std::vector<std::pair<U3, int>> irreps() const` | The resulting U(3) irreps with their level dimensionalities, sorted by the irrep labels. |
| `int multiplicity(const U3& i) const` | Level dimensionality of a single U(3) irrep `i`; 0 if `i` is not contained in the reduced irrep. |

The member types `u3_irrep` and `u4_irrep` are `U3` and `U4`. By default, irreps
are represented by two simple types:

| Type | Description |
|---|---|
| `untou3::U3Irrep` | U(3) irrep [e1,e2,e3] with members `e1`, `e2`, `e3`, `dim()`, comparison operators, and `operator<<` printing `[e1,e2,e3]`. |
| `untou3::U4Irrep` | U(4) irrep [f1,f2,f3,f4] with members `f1`, ..., `f4`, comparison operators, and `operator<<`. |

One `UNtoU3` object can be reused for any number of reductions, also for
different HO levels. `L = 4` is the case described in the paper and covered by
the test programs; the template parameter allows other values, but these are
untested.

### Using other irrep types

`UNtoU3` can work directly with irrep types of another library. It accesses
irreps only through `untou3::IrrepTraits`, which is specialized for `U3Irrep`
and `U4Irrep`. For other types, provide specializations with these static
member functions:

- for a U(3) irrep type `T`: `std::array<int, 3> labels(const T&)`, returning
  {e1, e2, e3}, and `T make(int e1, int e2, int e3)`, creating the irrep;
- for a U(4) irrep type `T`: `std::array<int, 4> labels(const T&)`, returning
  {f1, f2, f3, f4}.

The types do not need to be default-constructible or comparable. The class
template checks these requirements with the concepts `untou3::U3IrrepType` and
`untou3::U4IrrepType`.

For example, the following makes `UNtoU3` work with the irrep types `u3::irrep`
and `u4::irrep` of the QCX libraries u3lib and u4lib, whose labels are accessed
by `n13()`, …, `n33()` and `n14()`, …, `n44()`:

```cpp
#include <U3basis.h>
#include <U4basis.h>
#include <UNtoU3.h>

template <>
struct untou3::IrrepTraits<u3::irrep>
{
  static std::array<int, 3> labels(const u3::irrep& i)
  {
    return {i.n13(), i.n23(), i.n33()};
  }

  static u3::irrep make(int e1, int e2, int e3) { return {e1, e2, e3}; }
};

template <>
struct untou3::IrrepTraits<u4::irrep>
{
  static std::array<int, 4> labels(const u4::irrep& i)
  {
    return {i.n14(), i.n24(), i.n34(), i.n44()};
  }
};

// UNtoU3 with u3lib and u4lib irrep types:
template <int L>
using UNtoU3 = untou3::UNtoU3<L, u3::irrep, u4::irrep>;
```

With this, `irreps()` returns `u3::irrep` objects, and `reduceU4()` and
`multiplicity()` accept `u4::irrep` and `u3::irrep`. `src/test_custom_irreps.cpp`
shows the same with simple user-defined types.

## Building the test programs

The repository contains four test/benchmark programs in `src/`:

- `test_10401.cpp` reduces the U(6) irrep [4,2,2,2,2,0] and prints the resulting
  U(3) irreps with their level dimensionalities.
- `test_u4_5511.cpp` performs the conjugated reduction of the U(4) irrep
  [5,5,1,1] for N = 6.
- `test_custom_irreps.cpp` performs both reductions with user-defined irrep
  types plugged in by `IrrepTraits` and checks that the results agree with the
  default types.
- `test_input.cpp` reads η and R = [r4,r3,r2,r1,r0] from the standard input,
  performs the reduction, and prints the total dimension of the resulting U(3)
  irreps.

The build is configured by environment variables: `UNTOU3_DIR` must point to the
repository root, and the compiler flags, including `-std=c++20`, are taken from
`CXXFLAGS`. The programs are built into `bin/`:

```bash
export UNTOU3_DIR=$PWD CXXFLAGS="-std=c++20 -O2"
make
bin/test_10401
bin/test_u4_5511
bin/test_custom_irreps
echo "2 1 0 4 0 1" | bin/test_input
```

The last command prints `U(3) irreps total dim = 405`. Optional settings:

- `HAVE_BOOST=1` (and `BOOST_ROOT`, if Boost is not in a default location) makes
  `test_input` also compute dim[f] analytically and print `U(N) irrep dim = 405`.
- Adding `-DMEASURE_TIME` to `CXXFLAGS` makes `test_input` print the runtime of the
  reduction.

Run `make clean` to remove the built programs.

## Repository layout

| Path | Contents |
|---|---|
| `include/UNtoU3.h` | The library (single header). |
| `src/` | Test and benchmark programs, with their `Makefile`. |
| `Makefile`, `config.mk` | Top-level build of the test programs. |

## Algorithm and publication

The method is described in:

> D. Langr and T. Dytrych: **Efficient Algorithm for U(N) to U(3) Representation
> Reduction in Isospin-Adapted Nuclear Structure Calculations.** In: R. Wyrzykowski
> et al. (eds.), *Parallel Processing and Applied Mathematics (PPAM 2024)*, Lecture
> Notes in Computer Science, vol. 15581, pp. 233–243, Springer, Cham, 2025.
> [doi:10.1007/978-3-031-85703-4_16](https://doi.org/10.1007/978-3-031-85703-4_16)

In brief, the reduction is based on generating all Gelfand patterns of the input
U(N) irrep, each of which labels a basis state and determines a U(3) weight. The
patterns are generated recursively row by row. The UNtoU3-4 algorithm:

- represents each Gelfand pattern row by the fixed-size array R of the counts of
  its labels, which avoids dynamic memory allocations;
- uses the "binary form" B of R (which of the counts are nonzero) as a key into a
  small, precomputed database of all possible differences between successive
  rows; for Λ = 4 there are only 2⁵ = 32 keys and 164 differences in total, so the
  database is cache-friendly;
- stores the contribution of each difference to the U(3) weight in the database
  as well, and accumulates partial inner products along the recursion instead of
  building complete weight vectors.

The paper evaluates the algorithm on a node of the Frontera supercomputer (TACC).
For Λ = 2 inputs, it is about 3× slower than the earlier, specialized UNtoU3-2
algorithm, which, however, cannot handle Λ = 4. For Λ = 4 inputs, it reduces U(N)
irreps with dimensions of up to about 7·10⁹ in a few minutes.

### Citing

If you use UNtoU3 in your work, please cite the paper above:

```bibtex
@InProceedings{Langr2025UNtoU3,
  author    = {Langr, Daniel and Dytrych, Tom{\'a}{\v{s}}},
  title     = {Efficient Algorithm for {U(N)} to {U(3)} Representation Reduction
               in Isospin-Adapted Nuclear Structure Calculations},
  booktitle = {Parallel Processing and Applied Mathematics (PPAM 2024)},
  editor    = {Wyrzykowski, Roman and others},
  series    = {Lecture Notes in Computer Science},
  volume    = {15581},
  pages     = {233--243},
  publisher = {Springer},
  address   = {Cham},
  year      = {2025},
  doi       = {10.1007/978-3-031-85703-4_16}
}
```

## Related work

- **UNtoU3-2**, the predecessor algorithm for Λ = 2 (protons and neutrons treated
  as distinct particles, as in the LSU3shell code): D. Langr, T. Dytrych,
  J. P. Draayer, K. D. Launey, P. Tvrdík, *Efficient algorithm for representations
  of U(3) in U(N)*, Comput. Phys. Commun. 244, 442–447 (2019),
  [doi:10.1016/j.cpc.2019.05.018](https://doi.org/10.1016/j.cpc.2019.05.018).
- **UNTOU3**, the original FORTRAN implementation of the U(N) → U(3) reduction:
  J. P. Draayer, Y. Leschber, S. C. Park, R. Lopez, *Representations of U(3) in
  U(N)*, Comput. Phys. Commun. 56, 279–290 (1989),
  [doi:10.1016/0010-4655(89)90024-6](https://doi.org/10.1016/0010-4655(89)90024-6).

## Branches

- `main`: the library described here.
- `qcx` (**deprecated**, to be removed): an earlier variant integrated with the
  QCX libraries, which uses their irrep types `u3::irrep` and `u4::irrep`. The
  same is now possible with `main` and `IrrepTraits` (see
  [Using other irrep types](#using-other-irrep-types)), which does not require
  a separate branch. The `qcx` branch is kept for a transition period only.

## License

UNtoU3 is released under the [MIT License](LICENSE).
