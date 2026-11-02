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
#include <cstdint>
#include <iostream>

#include <UNtoU3.h>

int main()
{
  // maximum U(N) label: 4
  untou3::UNtoU3<4> gen;

  // HO level eta = 2, i.e., N = (eta+1)(eta+2)/2 = 6
  gen.generateXYZ(2);

  // input U(6) irrep [4,2,2,2,2,0] given by R = [r4,r3,r2,r1,r0] = [1,0,4,0,1]
  gen.generateU3Weights({1, 0, 4, 0, 1}, 6);

  uint64_t total_dim = 0;

  // iterate over the generated U(3) weights; the level dimensionality is
  // nonzero only for those that are U(3) irreps
  for (const auto& [weight, mult] : gen.multMap())
    if (auto D_l = gen.getLevelDimensionality(weight))
    {
      uint64_t e1 = weight[0], e2 = weight[1], e3 = weight[2];
      uint64_t dim = (e1 - e2 + 1) * (e1 - e3 + 2) * (e2 - e3 + 1) / 2;

      std::cout << "[" << e1 << "," << e2 << "," << e3 << "] : D_l = " << D_l
                << ", dim = " << dim << "\n";

      total_dim += D_l * dim;
    }

  std::cout << "total dim = " << total_dim << "\n";
}
```

Compiled with `g++ -std=c++20 -O2 -I include example.cpp`, it prints the
following (the order of lines is unspecified; sorted here):

```text
[8,8,8] : D_l = 1, dim = 1
[9,8,7] : D_l = 1, dim = 8
[10,8,6] : D_l = 3, dim = 27
[10,9,5] : D_l = 1, dim = 35
[10,10,4] : D_l = 1, dim = 28
[11,7,6] : D_l = 1, dim = 35
[11,8,5] : D_l = 1, dim = 64
[12,6,6] : D_l = 1, dim = 28
[12,8,4] : D_l = 1, dim = 125
total dim = 405
```

which agrees with dim[4,2,2,2,2,0] = 405.

### API

All functionality is provided by the class template `untou3::UNtoU3<L>`, where
`L` is the maximum U(N) label.

| Member | Description |
|---|---|
| `void generateXYZ(int n)` | Prepares the HO quanta of the `n`-th HO level. Must be called first. |
| `void generateU3Weights(R, uint64_t N)` | Generates all U(3) weights of the U(N) irrep given by `R = {r_L, ..., r_0}` with their multiplicities; `N` must be equal to (n+1)(n+2)/2. |
| `const U3MultMap& multMap() const` | The generated U(3) weights (`std::array<uint32_t, 3>`) and their multiplicities. |
| `uint64_t getLevelDimensionality(const U3Weight& w) const` | Level dimensionality D_l of `w`, nonzero only if `w` is a U(3) irrep contained in [f]; `w` must be taken from `multMap()`. |

Notes:

- Results of successive `generateU3Weights` calls accumulate. Use a new
  `UNtoU3` object for each input U(N) irrep.
- `multMap()` is an unordered map; its iteration order is unspecified.
- `L = 4` is the case described in the paper and covered by the test programs.
  The template parameter allows other values, but these are untested.

## Building the test programs

The repository contains two test/benchmark programs in `src/`:

- `test_10401.cpp` reduces the U(6) irrep [4,2,2,2,2,0] and prints the resulting
  U(3) irreps with their level dimensionalities.
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

- `main`: the standalone library described here.
- `qcx`: a variant integrated with the QCX libraries (u3lib, u4lib). It adds a
  higher-level interface based on their U(3) and U(4) irrep types, including the
  conjugated U(4) → U(3) reduction.

## License

UNtoU3 is released under the [MIT License](LICENSE).
