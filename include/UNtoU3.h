#ifndef UNTOU3_UNTOU3_H
#define UNTOU3_UNTOU3_H

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numeric>
#include <ostream>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

namespace untou3
{
  // U(3) irrep [e1,e2,e3], where e1 >= e2 >= e3 >= 0
  struct U3Irrep
  {
    int e1, e2, e3;

    // irrep's dimension:
    int dim() const
    {
      return (e1 - e2 + 1) * (e1 - e3 + 2) * (e2 - e3 + 1) / 2;
    }

    auto operator<=>(const U3Irrep&) const = default;
  }; // struct U3Irrep

  inline std::ostream& operator<<(std::ostream& s, const U3Irrep& i)
  {
    return s << "[" << i.e1 << "," << i.e2 << "," << i.e3 << "]";
  }

  // U(4) irrep [f1,f2,f3,f4], where f1 >= f2 >= f3 >= f4 >= 0
  struct U4Irrep
  {
    int f1, f2, f3, f4;

    auto operator<=>(const U4Irrep&) const = default;
  }; // struct U4Irrep

  inline std::ostream& operator<<(std::ostream& s, const U4Irrep& i)
  {
    return s << "[" << i.f1 << "," << i.f2 << "," << i.f3 << "," << i.f4
             << "]";
  }

  // Traits that tell UNtoU3 how to work with a type representing U(3) or U(4)
  // irreps. To use UNtoU3 with other irrep types than U3Irrep and U4Irrep,
  // specialize IrrepTraits for them:
  //
  // - for a U(3) irrep type T:
  //     static std::array<int, 3> labels(const T& i); // returns {e1,e2,e3}
  //     static T make(int e1, int e2, int e3);        // creates [e1,e2,e3]
  //
  // - for a U(4) irrep type T:
  //     static std::array<int, 4> labels(const T& i); // returns {f1,f2,f3,f4}
  template <typename T>
  struct IrrepTraits;

  template <>
  struct IrrepTraits<U3Irrep>
  {
    static std::array<int, 3> labels(const U3Irrep& i)
    {
      return {i.e1, i.e2, i.e3};
    }

    static U3Irrep make(int e1, int e2, int e3) { return {e1, e2, e3}; }
  }; // struct IrrepTraits<U3Irrep>

  template <>
  struct IrrepTraits<U4Irrep>
  {
    static std::array<int, 4> labels(const U4Irrep& i)
    {
      return {i.f1, i.f2, i.f3, i.f4};
    }
  }; // struct IrrepTraits<U4Irrep>

  // type T can represent U(3) irreps in UNtoU3:
  template <typename T>
  concept U3IrrepType = requires(const T& i, int e) {
    { IrrepTraits<T>::labels(i) } -> std::convertible_to<std::array<int, 3>>;
    { IrrepTraits<T>::make(e, e, e) } -> std::convertible_to<T>;
  };

  // type T can represent U(4) irreps in UNtoU3:
  template <typename T>
  concept U4IrrepType = requires(const T& i) {
    { IrrepTraits<T>::labels(i) } -> std::convertible_to<std::array<int, 4>>;
  };

  // L = maximum value of U(N) label
  // U3, U4 = types representing U(3) and U(4) irreps (see IrrepTraits)
  template <int L, U3IrrepType U3 = U3Irrep, U4IrrepType U4 = U4Irrep>
  class UNtoU3
  {
  public:
    using u3_irrep = U3; // type representing U(3) irreps
    using u4_irrep = U4; // type representing U(4) irreps

    // Type for representation of a Gelfand pattern row.
    //
    // r[0] is the number of 'L's in the row.
    // r[1] is the nubmer of 'L-1's in the row.
    // ...
    // r[L] is the number of '0's in the row
    //
    // Example: L = 4, r = { 1, 0, 4, 0, 1 }
    // represesents Gelfand pattern row [ 4, 2, 2, 2, 2, 0 ]
    using GPR_rep = std::array<int16_t, L + 1>;

  private:
    using U3Weight = std::array<uint32_t, 3>; // representation of U(3) weight

    // helper haseher class; need to be here:
    struct U3Weight_hasher
    {
      // Boost-motivated:
      template <typename T>
      static void hash_combine(std::size_t& seed, const T& v)
      {
        seed ^= std::hash<T>{}(v) + 0x9E3779B9 + (seed << 6) + (seed >> 2);
      }

      size_t operator()(const U3Weight& key) const
      {
        size_t res = 0;

        for (auto w : key)
          hash_combine(res, w);

        return res;
      }
    }; // class U3Weight_hasher

    // type for table for storing U(3) weights and their multiplicities:
    using U3MultMap = std::unordered_map<U3Weight, uint32_t, U3Weight_hasher>;

    // type of elements in resulting table (difference + weight vector element)
    using Element_type = std::pair<GPR_rep, uint32_t>;

    // database of "differences" between Galfend pattern rows:
    class Diffs
    {
    public:
      Diffs() { generate(); }

      // set of differences for binary mask of Gelfand pattern row rep:
      std::span<Element_type> get_bin(uint64_t bin)
      {
        assert(bin < indexes_.size());
        assert(bin < counts_.size());

        auto first = indexes_[bin];
        auto count = counts_[bin];

        if (count == 0)
          return {};
        else
          return {differences_.data() + first, count};
      }

      // set of differences for Gelfand pattern row representation:
      std::span<Element_type> get_gpr(const GPR_rep& gpr)
      {
        auto bin = gpr_to_bin(gpr);
        return get_bin(bin);
      }

    private:
      int64_t get_weight(const GPR_rep& diff)
      {
        int64_t res = 0;

        for (int k = 0; k <= L; k++)
          res += diff[k] * (L - k) * (-1);

        assert(res >= 0);
        assert(res <= L);

        return res;
      }

      void generate_rules_recursive(const uint64_t bin, const GPR_rep& gpr,
                                    GPR_rep diff, int first)
      {
        // find second nonzero:
        int second = first + 1;
        while ((second < L + 1) && (gpr[second] == 0))
          second++;

        // no second nonzero:
        if (second == L + 1)
        {
          counts_[bin]++;
          differences_.emplace_back(diff, get_weight(diff));

          return;
        }

        for (int i = first; i <= second; i++)
        {
          diff[i]++;
          generate_rules_recursive(bin, gpr, diff, second);
          diff[i]--;
        }
      }

      void generate_rules(const uint64_t bin)
      {
        auto gpr = bin_to_gpr(bin);

        GPR_rep diff;
        for (int k = 0; k < L + 1; k++)
          diff[k] = gpr[k] ? -1 : 0;

        // find first nonzero:
        int first = 0;
        while ((first < L + 1) && (gpr[first] == 0))
          first++;

        // no nonzeros:
        if (first == L + 1)
        {
          indexes_[bin] = differences_.size();
          counts_[bin] = 0;

          return;
        }

        // find second nonzero:
        int second = first + 1;
        while ((second < L + 1) && (gpr[second] == 0))
          second++;

        // no second nonzero:
        if (second == L + 1)
        {
          indexes_[bin] = differences_.size();
          counts_[bin] = 1;
          differences_.emplace_back(diff, get_weight(diff));

          return;
        }

        indexes_[bin] = differences_.size();

        for (int i = first; i <= second; i++)
        {
          diff[i]++;
          generate_rules_recursive(bin, gpr, diff, second);
          diff[i]--;
        }
      }

      void generate()
      {
        uint64_t max = 1UL << (L + 1);

        indexes_.resize(max, 0);
        counts_.resize(max, 0);

        for (uint64_t bin = 0; bin < max; bin++)
          generate_rules(bin);
      }

      std::vector<Element_type> differences_; // differences data

      std::vector<size_t> indexes_; // index of first difference
      std::vector<size_t> counts_;  // count of differences
    }; // class Diffs

    enum { NZ, NX, NY }; // definition of the order of axis for weight vectors

    // HO level set by set_shell (-1 = not set yet):
    int n_ = -1;

    // HO quanta vectors generated by generateXYZ:
    std::array<std::vector<uint32_t>, 3> xyz_;

    // table of resulting U(3) irreps and their multiplicities:
    U3MultMap mult_;

    // differences for generation of lower Gelfand pattern rows:
    Diffs diffs_;

    // Generates HO quanta vectors for given nth HO level.
    // Need to be used before generateU3Weights() member function is called.
    void generateXYZ(int n)
    {
      for (auto& v : xyz_)
        v.clear();

      for (int k = 0; k <= n; k++)
      {
        uint32_t nz = n - k;
        for (int nx = k; nx >= 0; nx--)
        {
          xyz_[NZ].push_back(nz);
          xyz_[NX].push_back(nx);
          xyz_[NY].push_back(k - nx);
        }
      }
    }

    // Generates U(3) weights and their multiplicities for an input U(N) irrep
    // [f]. [f] is specified as a top row of Gelfand pattern [n_L,...,n_1,n_0]
    // N=n_L+...+n_0 must be equal to (n+1)*(n+2)/2, where n was used as an
    // argument of set_shell.
    void generateU3Weights(const GPR_rep& gpr, uint64_t N,
                           U3Weight weight = {0, 0, 0})
    {
      assert(static_cast<uint64_t>(
                 std::accumulate(std::begin(gpr), std::end(gpr), 0)) == N);

      auto span = diffs_.get_gpr(gpr);

      if (span.size() == 0)
      {
        mult_[weight]++;
        return;
      }

      for (const auto& e : span)
      {
        auto lower_gpr = add_gpr(gpr, e.first);
        generateU3Weights(lower_gpr, N - 1,
                          {weight[0] + e.second * xyz_[0][N - 1],
                           weight[1] + e.second * xyz_[1][N - 1],
                           weight[2] + e.second * xyz_[2][N - 1]});
      }
    }

    // Get level dimensionality for a given U(3) weight passed as an argument.
    // Returns 0 if the weight is not a U(3) irrep contained in the reduced
    // U(N) irrep.
    uint64_t getLevelDimensionality(const U3Weight& weight) const
    {
      uint32_t f1 = weight[0], f2 = weight[1], f3 = weight[2];
      if ((f1 < f2) || (f2 < f3))
        return 0;

      // if [f1,f2,f3] is not a weight, neither is any higher weight:
      auto u3_mult = mult_.find({f1, f2, f3});
      if (u3_mult == mult_.end())
        return 0;
      auto mult = u3_mult->second;

      u3_mult = mult_.find({f1 + 1, f2 + 1, f3 - 2});
      mult += (u3_mult == mult_.end()) ? 0 : u3_mult->second;

      u3_mult = mult_.find({f1 + 2, f2 - 1, f3 - 1});
      mult += (u3_mult == mult_.end()) ? 0 : u3_mult->second;

      u3_mult = mult_.find({f1 + 2, f2, f3 - 2});
      mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

      u3_mult = mult_.find({f1 + 1, f2 - 1, f3});
      mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

      u3_mult = mult_.find({f1, f2 + 1, f3 - 1});
      mult -= (u3_mult == mult_.end()) ? 0 : u3_mult->second;

      return mult;
    }

    // addition of 2 Gelfand pattern row representations:
    static inline GPR_rep add_gpr(const GPR_rep& a, const GPR_rep& b)
    {
      GPR_rep res;

      for (int k = 0; k <= L; k++)
        res[k] = a[k] + b[k];

      return res;
    }

    // transformation of binary mask into Gelfand pattern row representation:
    static inline GPR_rep bin_to_gpr(uint64_t n)
    {
      GPR_rep res;

      for (int k = 0; k <= L; k++)
        res[L - k] = (n >> k) & 0x01;

      return res;
    }

    // transformation of Gelfand pattern row representation to binary mask:
    static inline uint64_t gpr_to_bin(const GPR_rep& gpr)
    {
      uint64_t res = 0;

      for (int k = 0; k <= L; k++)
        if (gpr[L - k] != 0)
          res |= 1UL << k;

      return res;
    }

  public:
    // Sets the nth HO level, N = (n+1)*(n+2)/2. Needs to be called before
    // reduce() or reduceU4().
    void set_shell(int n)
    {
      generateXYZ(n);
      n_ = n;
    }

    // HO level set by set_shell() (-1 if not set yet).
    int shell() const { return n_; }

    // U(N) to U(3) reduction of an input U(N) irrep [f] given by its
    // representation gpr = {r_L,...,r_0}, where r_k is the number of labels
    // equal to k and r_L+...+r_0 = N. Replaces results of previous reductions.
    void reduce(const GPR_rep& gpr)
    {
      assert(n_ >= 0);

      mult_.clear();

      auto N = (n_ + 1) * (n_ + 2) / 2;

      assert(std::accumulate(std::begin(gpr), std::end(gpr), 0) == N);

      generateU3Weights(gpr, N);
    }

    // Conjugated U(4) to U(3) reduction: reduces the U(N) irrep conjugate to
    // the input U(4) irrep i (the transposed Young diagram), where N is given
    // by the HO level set by set_shell().
    void reduceU4(const U4& i)
      requires(L == 4)
    {
      assert(n_ >= 0);

      int N = (n_ + 1) * (n_ + 2) / 2;

      auto [f1, f2, f3, f4] = IrrepTraits<U4>::labels(i);

      assert((f1 >= f2) && (f2 >= f3) && (f3 >= f4) && (f4 >= 0));
      assert(f1 <= N);

      GPR_rep rep = {static_cast<int16_t>(f4), static_cast<int16_t>(f3 - f4),
                     static_cast<int16_t>(f2 - f3),
                     static_cast<int16_t>(f1 - f2),
                     static_cast<int16_t>(N - f1)};

      reduce(rep);
    }

    // The same as reduceU4(i), with N passed explicitly; N must be equal to
    // (n+1)*(n+2)/2 for the HO level n set by set_shell().
    void reduceU4(const U4& i, [[maybe_unused]] int N)
      requires(L == 4)
    {
      assert(N == (n_ + 1) * (n_ + 2) / 2);

      reduceU4(i);
    }

    // Level dimensionality (multiplicity) of the U(3) irrep i in the reduced
    // U(N) irrep. Returns 0 if i is not contained in it.
    int multiplicity(const U3& i) const
    {
      auto [e1, e2, e3] = IrrepTraits<U3>::labels(i);

      if (e3 < 0)
        return 0;

      U3Weight w = {static_cast<uint32_t>(e1), static_cast<uint32_t>(e2),
                    static_cast<uint32_t>(e3)};

      return getLevelDimensionality(w);
    }

    // Resulting U(3) irreps of the last reduction together with their level
    // dimensionalities (multiplicities), sorted by the U(3) irrep labels.
    std::vector<std::pair<U3, int>> irreps() const
    {
      std::vector<std::pair<U3Weight, int>> found;

      for (const auto& [w, m] : mult_)
        if (auto Dl = getLevelDimensionality(w); Dl > 0)
          found.emplace_back(w, Dl);

      // sorting by labels does not require comparison of U3 objects:
      std::sort(found.begin(), found.end());

      std::vector<std::pair<U3, int>> res;
      res.reserve(found.size());

      for (const auto& [w, Dl] : found)
        res.emplace_back(IrrepTraits<U3>::make(static_cast<int>(w[0]),
                                               static_cast<int>(w[1]),
                                               static_cast<int>(w[2])),
                         Dl);

      return res;
    }
  }; // class UNtoU3
} // namespace untou3

#endif
