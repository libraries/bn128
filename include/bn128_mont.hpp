#ifndef BN128_H_
#define BN128_H_

#include <intx/intx.hpp>

namespace bn128 {

// Maybe there is a better way to implement this macro, but this is enough for now.
#ifdef __riscv
#undef assert
#define assert(x)                                                                                                      \
  if (!x)                                                                                                              \
  exit(2)
#endif

using uint256 = intx::uint256;
using uint512 = intx::uint512;

inline uint256 _addmod(const uint256 &x, const uint256 &y, const uint256 &n) { return intx::addmod(x, y, n); }

inline uint256 _submod(const uint256 &x, const uint256 &y, const uint256 &n) { return _addmod(x, n - y, n); }

inline uint256 _negmod(const uint256 &x, const uint256 &n) { return n - x; }

inline uint256 _mulmod(const uint256 &x, const uint256 &y, const uint256 &n) { return intx::mulmod(x, y, n); }

// Extended euclidean algorithm to find modular inverses for integers.
// The return value INV_X satisfies: (X * INV_X) % N = 1
inline uint256 _invmod(const uint256 &x, const uint256 &n) {
  uint256 lm = 1, hm = 0;
  uint256 lo = x % n, hi = n;
  while (lo > 1) {
    uint256 r = hi / lo;
    uint256 nm = _submod(hm, _mulmod(lm, r, n), n);
    uint256 nw = hi - lo * r;
    hm = lm;
    lm = nm;
    hi = lo;
    lo = nw;
  }
  return lm % n;
}

inline uint256 _divmod(const uint256 &x, const uint256 &y, const uint256 &n) { return _mulmod(x, _invmod(y, n), n); }

inline uint256 _powmod(const uint256 &x, const uint256 &y, const uint256 &n) {
  if (y == 0) {
    return 1;
  } else if (y == 1) {
    return x;
  } else if (y & 1) {
    return _mulmod(_powmod(_mulmod(x, x, n), y >> 1, n), x, n);
  } else {
    return _powmod(_mulmod(x, x, n), y >> 1, n);
  }
}

constexpr inline uint256 h256(const char *s) { return intx::from_string<uint256>(s); }

// =====================================================================================================================
// EIP 196 👇
// =====================================================================================================================

// The prime modulus of the field.
#define HEX_FIELD_MODULUS "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47"
// R = 2 ** 256
// R_SQUARD = R * R % FIELD_MODULUS
#define HEX_R_SQUARD "0x06d89f71cab8351f47ab1eff0a417ff6b5e71911d44501fbf32cfc5b538afa89"
// R_CUBED = R_SQUARD * R % FIELD_MODULUS
#define HEX_R_CUBED "0x20fd6e902d592544ef7f0b0c0ada0afb62f210e6a7283db6b1cd6dafda1530df"
// R_PRIME * R % FIELD_MODULUS == 1
#define HEX_R_PRIME "0x2e67157159e5c639cf63e9cfb74492d9eb2022850278edf8ed84884a014afa37"
// FIELD_MODULUS_PRIME * (-FIELD_MODULUS) % R == 1
#define HEX_FIELD_MODULUS_PRIME "0xf57a22b791888c6bd8afcbd01833da809ede7d651eca6ac987d20782e4866389"
// FQ_ONE = mont_encode(1);
#define HEX_FQ_ONE "0x0e0a77c19a07df2f666ea36f7879462c0a78eb28f5c70b3dd35d438dc58f0d9d"
// FQ_NON_RESIDUE = mont_encode(FIELD_MODULUS - 1);
#define HEX_FQ_NON_RESIDUE "0x2259d6b14729c0fa51e1a247090812318d087f6872aabf4f68c3488912edefaa"
// G1_COEFF = mont_encode(3)
#define HEX_G1_COEFF "0x2a1f6744ce179d8e334bea4e696bd2841f6ac17ae15521b97a17caa950ad28d7"
#define HEX_G1_ONE_1 "0x1c14ef83340fbe5eccdd46def0f28c5814f1d651eb8e167ba6ba871b8b1e1b3a"

constexpr uint256 FIELD_MODULUS = h256(HEX_FIELD_MODULUS);
constexpr uint256 R_SQUARD = h256(HEX_R_SQUARD);
constexpr uint256 R_CUBED = h256(HEX_R_CUBED);
constexpr uint256 R_PRIME = h256(HEX_R_PRIME);
constexpr uint256 FIELD_MODULUS_PRIME = h256(HEX_FIELD_MODULUS_PRIME);

// Montgomery reduction, also known as REDC.
// REDC(T)=T*R' mod N(N>1)，
uint256 REDC(uint256 T) {
  uint512 t = uint512(T);
  uint512 n = uint512(FIELD_MODULUS);
  uint512 n_prime = uint512(FIELD_MODULUS_PRIME);

  uint512 m = uint512(intx::lo_half(t * n_prime));
  uint256 r = intx::hi_half(t + m * n);

  if (r >= FIELD_MODULUS) {
    r -= FIELD_MODULUS;
  }
  return r;
}

inline uint256 mont_encode(const uint256 &x) { return REDC(_mulmod(x, R_SQUARD, FIELD_MODULUS)); }

inline uint256 mont_decode(const uint256 &x) { return REDC(_mulmod(x, 1, FIELD_MODULUS)); }

struct FQ {
  uint256 c0;

  constexpr FQ() {}

  constexpr FQ(uint256 x) { c0 = x; }

  inline FQ inv() const { return REDC(_mulmod(_invmod(c0, FIELD_MODULUS), R_CUBED, FIELD_MODULUS)); }

  inline FQ pow(const uint256 &y) const { return FQ{c0 : _powmod(c0, y, FIELD_MODULUS)}; }

  inline FQ square() const { return FQ{c0 : REDC(_mulmod(c0, c0, FIELD_MODULUS))}; }
};

inline FQ operator+(const FQ &x, const FQ &y) { return FQ{c0 : _addmod(x.c0, y.c0, FIELD_MODULUS)}; }

inline FQ operator-(const FQ &x, const FQ &y) { return FQ{c0 : _submod(x.c0, y.c0, FIELD_MODULUS)}; }

inline FQ operator-(const FQ &x) { return FQ{c0 : _negmod(x.c0, FIELD_MODULUS)}; }

inline FQ operator*(const FQ &x, const FQ &y) { return FQ{c0 : REDC(_mulmod(x.c0, y.c0, FIELD_MODULUS))}; }

inline FQ operator/(const FQ &x, const FQ &y) { return FQ{c0 : _divmod(x.c0, y.c0, FIELD_MODULUS)}; }

inline bool operator==(const FQ &x, const FQ &y) { return x.c0 == y.c0; }

inline bool operator!=(const FQ &x, const FQ &y) { return x.c0 != y.c0; }

constexpr FQ FQ_ZERO = FQ(0);
constexpr FQ FQ_ONE = FQ(h256(HEX_FQ_ONE));
constexpr FQ FQ_NON_RESIDUE = FQ(h256(HEX_FQ_NON_RESIDUE));
constexpr FQ G1_COEFF_B = FQ(h256(HEX_G1_COEFF));

struct G1Affine;
struct G1;

struct G1Affine {
  FQ x;
  FQ y;

  G1 into() const;
};

struct G1 {
  FQ x;
  FQ y;
  FQ z;

  G1Affine affine() const;
  G1 doubl2() const;
  G1 mul(const uint256 &c) const;
};

constexpr G1 G1_ZERO = G1{
  x : FQ_ZERO,
  y : FQ_ONE,
  z : FQ_ZERO,
};

constexpr G1 G1_ONE = G1{
  x : FQ_ONE,
  y : FQ(h256(HEX_G1_ONE_1)),
  z : FQ_ONE,
};

G1 G1Affine::into() const {
  FQ a = x;
  FQ b = y;
  FQ c = (x == FQ_ZERO && y == FQ_ZERO) ? FQ_ZERO : FQ_ONE;
  return G1{x : a, y : b, z : c};
}

G1Affine G1::affine() const {
  if (z == FQ_ZERO) {
    return G1Affine{x : FQ_ZERO, y : FQ_ZERO};
  } else if (z == FQ_ONE) {
    return G1Affine{x : x, y : y};
  } else {
    FQ zinv = z.inv();
    FQ zinv_squared = zinv.square();
    return G1Affine{x : x * zinv_squared, y : y * (zinv_squared * zinv)};
  }
}

G1 G1::doubl2() const {
  FQ a = x.square();
  FQ b = y.square();
  FQ c = b.square();
  FQ d = (x + b).square() - a - c;
  d = d + d;
  FQ e = a + a + a;
  FQ f = e.square();
  FQ x3 = f - (d + d);
  FQ c8 = c + c;
  c8 = c8 + c8;
  c8 = c8 + c8;
  FQ yz = y * z;
  return G1{
    x : x3,
    y : e * (d - x3) - c8,
    z : yz + yz,
  };
}

G1 operator+(const G1 &p, const G1 &q) {
  FQ x1 = p.x, y1 = p.y, z1 = p.z;
  FQ x2 = q.x, y2 = q.y, z2 = q.z;
  if (z1 == FQ_ZERO) {
    return q;
  }
  if (z2 == FQ_ZERO) {
    return p;
  }
  FQ z1_squared = z1.square();
  FQ z2_squared = z2.square();
  FQ u1 = x1 * z2_squared;
  FQ u2 = x2 * z1_squared;
  FQ z1_cubed = z1 * z1_squared;
  FQ z2_cubed = z2 * z2_squared;
  FQ s1 = y1 * z2_cubed;
  FQ s2 = y2 * z1_cubed;
  if (u1 == u2 && s1 == s2) {
    return p.doubl2();
  }
  FQ h = u2 - u1;
  FQ s2_minus_s1 = s2 - s1;
  FQ i = (h + h).square();
  FQ j = h * i;
  FQ r = s2_minus_s1 + s2_minus_s1;
  FQ v = u1 * i;
  FQ s1_j = s1 * j;
  FQ x3 = r.square() - j - (v + v);
  return G1{
    x : x3,
    y : r * (v - x3) - (s1_j + s1_j),
    z : ((z1 + z2).square() - z1_squared - z2_squared) * h,
  };
}

G1 G1::mul(const uint256 &c) const {
  G1 r = G1_ZERO;
  bool found_one = 0;
  for (int i = 255; i > -1; i--) {
    if (found_one) {
      r = r.doubl2();
    }
    if (c & (uint256{1} << i)) {
      found_one = 1;
      r = r + *this;
    }
  }
  return r;
}

// =====================================================================================================================
// EIP 196 👆
// =====================================================================================================================
// EIP 197 👇
// =====================================================================================================================

// G2_COEFF_B0 = mont_encode(0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5)
// G2_COEFF_B1 = mont_encode(0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2)
#define HEX_G2_COEFF_B0 "0x2514c6324384a86d26b7edf049755260020b1b273633535d3bf938e377b802a8"
#define HEX_G2_COEFF_B1 "0x0141b9ce4a688d4dd749d0dd22ac00aa65f0b37d93ce0d3e38e7ecccd1dcff67"
#define HEX_G2_ONE_00 "0x19573841af96503bfbb8264797811adfdceb1935497b01728e83b5d102bc2026"
#define HEX_G2_ONE_01 "0x14fef0833aea7b6b09e950fc52a02f866043dd5a5802d8c4afb4737da84c6140"
#define HEX_G2_ONE_10 "0x28fd7eebae9e4206ff9e1a62231b7dfefe7fd297f59e9b78619dfa9d886be9f6"
#define HEX_G2_ONE_11 "0x0da4a0e693fd648255f935be33351076dc57f922327d3cbb64095b56c71856ee"

struct FQ2 {
  FQ c0;
  FQ c1;

  constexpr FQ2() {}

  constexpr FQ2(uint256 x, uint256 y) {
    c0 = FQ(x);
    c1 = FQ(y);
  }

  constexpr FQ2(FQ x, FQ y) {
    c0 = x;
    c1 = y;
  }

  FQ2 inv() const {
    FQ t = (c0.square() - c1.square() * FQ_NON_RESIDUE).inv();
    return FQ2{
      c0 : c0 * t,
      c1 : -(c1 * t),
    };
  }

  FQ2 square() const {
    FQ a = c0 * c1;
    return FQ2{
      c0 : (c1 * FQ_NON_RESIDUE + c0) * (c0 + c1) - a - a * FQ_NON_RESIDUE,
      c1 : a + a,
    };
  }
};

FQ2 operator+(const FQ2 &x, const FQ2 &y) {
  return FQ2{
    c0 : x.c0 + y.c0,
    c1 : x.c1 + y.c1,
  };
}

FQ2 operator-(const FQ2 &x, const FQ2 &y) {
  return FQ2{
    c0 : x.c0 - y.c0,
    c1 : x.c1 - y.c1,
  };
}

FQ2 operator-(const FQ2 &x) {
  return FQ2{
    c0 : -x.c0,
    c1 : -x.c1,
  };
}

FQ2 operator*(const FQ2 &x, const FQ2 &y) {
  FQ aa = x.c0 * y.c0;
  FQ bb = x.c1 * y.c1;
  return FQ2{
    c0 : bb * FQ_NON_RESIDUE + aa,
    c1 : (x.c0 + x.c1) * (y.c0 + y.c1) - aa - bb,
  };
}

bool operator==(const FQ2 &x, const FQ2 &y) { return x.c0 == y.c0 && x.c1 == y.c1; }

bool operator!=(const FQ2 &x, const FQ2 &y) { return x.c0 != y.c0 || x.c1 != y.c1; }

constexpr FQ2 FQ2_ZERO = FQ2{c0 : FQ_ZERO, c1 : FQ_ZERO};
constexpr FQ2 FQ2_ONE = FQ2{c0 : FQ_ONE, c1 : FQ_ZERO};
constexpr FQ2 G2_COEFF_B = FQ2{c0 : FQ(h256(HEX_G2_COEFF_B0)), c1 : FQ(h256(HEX_G2_COEFF_B1))};

struct G2Affine;
struct G2;

struct G2Affine {
  FQ2 x;
  FQ2 y;

  G2 into() const;
};

struct G2 {
  FQ2 x;
  FQ2 y;
  FQ2 z;

  G2Affine affine() const;
  G2 doubl2() const;
  G2 mul(const uint256 &c) const;
};

constexpr G2 G2_ZERO = G2{
  x : FQ2_ZERO,
  y : FQ2_ONE,
  z : FQ2_ZERO,
};

constexpr G2 G2_ONE = G2{
  x : FQ2(h256(HEX_G2_ONE_00), h256(HEX_G2_ONE_01)),
  y : FQ2(h256(HEX_G2_ONE_10), h256(HEX_G2_ONE_11)),
  z : FQ2_ONE,
};

G2 G2Affine::into() const {
  FQ2 a = x;
  FQ2 b = y;
  FQ2 c = (x == FQ2_ZERO && y == FQ2_ZERO) ? FQ2_ZERO : FQ2_ONE;
  return G2{x : a, y : b, z : c};
}

bool operator==(const G2Affine &x, const G2Affine &y) { return x.x == y.x && x.y == y.y; }
bool operator!=(const G2Affine &x, const G2Affine &y) { return x.x != y.x || x.y != y.y; }

G2Affine G2::affine() const {
  if (z == FQ2_ZERO) {
    return G2Affine{x : FQ2_ZERO, y : FQ2_ZERO};
  } else if (z == FQ2_ONE) {
    return G2Affine{x : x, y : y};
  } else {
    FQ2 zinv = z.inv();
    FQ2 zinv_squared = zinv.square();
    return G2Affine{x : x * zinv_squared, y : y * (zinv_squared * zinv)};
  }
}

G2 G2::doubl2() const {
  FQ2 a = x.square();
  FQ2 b = y.square();
  FQ2 c = b.square();
  FQ2 d = (x + b).square() - a - c;
  d = d + d;
  FQ2 e = a + a + a;
  FQ2 f = e.square();
  FQ2 x3 = f - (d + d);
  FQ2 c8 = c + c;
  c8 = c8 + c8;
  c8 = c8 + c8;
  FQ2 yz = y * z;
  return G2{
    x : x3,
    y : e * (d - x3) - c8,
    z : yz + yz,
  };
}

bool operator==(const G2 &x, const G2 &y) { return x.x == y.x && x.y == y.y && x.z == y.z; }
bool operator!=(const G2 &x, const G2 &y) { return x.x != y.x || x.y != y.y || x.z != y.z; }

G2 operator+(const G2 &p, const G2 &q) {
  FQ2 x1 = p.x, y1 = p.y, z1 = p.z;
  FQ2 x2 = q.x, y2 = q.y, z2 = q.z;
  if (z1 == FQ2_ZERO) {
    return q;
  }
  if (z2 == FQ2_ZERO) {
    return p;
  }
  FQ2 z1_squared = z1.square();
  FQ2 z2_squared = z2.square();
  FQ2 u1 = x1 * z2_squared;
  FQ2 u2 = x2 * z1_squared;
  FQ2 z1_cubed = z1 * z1_squared;
  FQ2 z2_cubed = z2 * z2_squared;
  FQ2 s1 = y1 * z2_cubed;
  FQ2 s2 = y2 * z1_cubed;
  if (u1 == u2 && s1 == s2) {
    return p.doubl2();
  }
  FQ2 h = u2 - u1;
  FQ2 s2_minus_s1 = s2 - s1;
  FQ2 i = (h + h).square();
  FQ2 j = h * i;
  FQ2 r = s2_minus_s1 + s2_minus_s1;
  FQ2 v = u1 * i;
  FQ2 s1_j = s1 * j;
  FQ2 x3 = r.square() - j - (v + v);
  return G2{
    x : x3,
    y : r * (v - x3) - (s1_j + s1_j),
    z : ((z1 + z2).square() - z1_squared - z2_squared) * h,
  };
}

G2 G2::mul(const uint256 &c) const {
  G2 r = G2_ZERO;
  bool found_one = 0;
  for (int i = 255; i > -1; i--) {
    if (found_one) {
      r = r.doubl2();
    }
    if (c & (uint256{1} << i)) {
      found_one = 1;
      r = r + *this;
    }
  }
  return r;
}

// =====================================================================================================================
// EIP 197 👆
// =====================================================================================================================

void alt_bn128_add(const uint256 p[2], const uint256 q[2], uint256 r[2]) {
  auto x_affine = G1Affine{
    x : FQ(mont_encode(p[0])),
    y : FQ(mont_encode(p[1])),
  };
  auto x = x_affine.into();
  auto y_affine = G1Affine{
    x : FQ(mont_encode(q[0])),
    y : FQ(mont_encode(q[1])),
  };
  auto y = y_affine.into();
  auto z = (x + y).affine();
  r[0] = mont_decode(z.x.c0);
  r[1] = mont_decode(z.y.c0);
}

void alt_bn128_mul(const uint256 p[2], const uint256 &n, uint256 r[2]) {
  auto x_affine = G1Affine{
    x : FQ(mont_encode(p[0])),
    y : FQ(mont_encode(p[1])),
  };
  auto x = x_affine.into();
  auto z = x.mul(n).affine();
  r[0] = mont_decode(z.x.c0);
  r[1] = mont_decode(z.y.c0);
}

} // namespace bn128

#endif /* BN128_H_ */
