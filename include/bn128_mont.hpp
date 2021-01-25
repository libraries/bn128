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

inline bool arreq(const uint256 *x, const uint256 *y, const int size) {
  switch (size) {
  case 12:
    if (y[11] != x[11])
      return 0;
  case 11:
    if (y[10] != x[10])
      return 0;
  case 10:
    if (y[9] != x[9])
      return 0;
  case 9:
    if (y[8] != x[8])
      return 0;
  case 8:
    if (y[7] != x[7])
      return 0;
  case 7:
    if (y[6] != x[6])
      return 0;
  case 6:
    if (y[5] != x[5])
      return 0;
  case 5:
    if (y[4] != x[4])
      return 0;
  case 4:
    if (y[3] != x[3])
      return 0;
  case 3:
    if (y[2] != x[2])
      return 0;
  case 2:
    if (y[1] != x[1])
      return 0;
  case 1:
    if (y[0] != x[0])
      return 0;
  }
  return 1;
}

// The prime modulus of the field.
#define FIELD_MODULUS_HEX "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47"
constexpr uint256 FIELD_MODULUS = h256(FIELD_MODULUS_HEX);
// R = 2 ** 256
// R_SQUARD = R * R % FIELD_MODULUS
// R_CUBED = R_SQUARD * R % FIELD_MODULUS
#define R_SQUARD_HEX "0x06d89f71cab8351f47ab1eff0a417ff6b5e71911d44501fbf32cfc5b538afa89"
#define R_CUBED_HEX "0x20fd6e902d592544ef7f0b0c0ada0afb62f210e6a7283db6b1cd6dafda1530df"
constexpr uint256 R_SQUARD = h256(R_SQUARD_HEX);
constexpr uint256 R_CUBED = h256(R_CUBED_HEX);
// R_PRIME * R % FIELD_MODULUS == 1
#define R_PRIME_HEX "0x2e67157159e5c639cf63e9cfb74492d9eb2022850278edf8ed84884a014afa37"
constexpr uint256 R_PRIME = h256(R_PRIME_HEX);
#define FIELD_MODULUS_PRIME_HEX "0xf57a22b791888c6bd8afcbd01833da809ede7d651eca6ac987d20782e4866389"
// FIELD_MODULUS_PRIME * (-FIELD_MODULUS) % R == 1
constexpr uint256 FIELD_MODULUS_PRIME = h256(FIELD_MODULUS_PRIME_HEX);

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

#define FQ_ONE_HEX "0x0e0a77c19a07df2f666ea36f7879462c0a78eb28f5c70b3dd35d438dc58f0d9d"
constexpr uint256 FQ_ONE = h256(FQ_ONE_HEX);

inline uint256 fq_add(const uint256 &x, const uint256 &y) { return _addmod(x, y, FIELD_MODULUS); }

inline uint256 fq_sub(const uint256 &x, const uint256 &y) { return _submod(x, y, FIELD_MODULUS); }

inline uint256 fq_mul(const uint256 &x, const uint256 &y) { return REDC(_mulmod(x, y, FIELD_MODULUS)); }

inline uint256 fq_inv(const uint256 &x) { return fq_mul(_invmod(x, FIELD_MODULUS), R_CUBED); }

inline uint256 fq_neg(const uint256 &x) { return _negmod(x, FIELD_MODULUS); }

// FQ_NON_RESIDUE = mont_encode(FIELD_MODULUS - 1);
#define FQ_NON_RESIDUE_HEX "0x2259d6b14729c0fa51e1a247090812318d087f6872aabf4f68c3488912edefaa"
constexpr uint256 FQ_NON_RESIDUE = h256(FQ_NON_RESIDUE_HEX);
constexpr uint256 FQ2_ONE[2] = {FQ_ONE, 0};

inline void fq2_add(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 a = fq_add(x[0], y[0]);
  uint256 b = fq_add(x[1], y[1]);
  r[0] = a;
  r[1] = b;
}

inline void fq2_sub(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 a = fq_sub(x[0], y[0]);
  uint256 b = fq_sub(x[1], y[1]);
  r[0] = a;
  r[1] = b;
}

inline void fq2_neg(const uint256 x[2], uint256 r[2]) {
  uint256 a = fq_neg(x[0]);
  uint256 b = fq_neg(x[1]);
  r[0] = a;
  r[1] = b;
}

inline void fq2_mul(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 aa = fq_mul(x[0], y[0]);
  uint256 bb = fq_mul(x[1], y[1]);
  uint256 c0 = fq_add(fq_mul(bb, FQ_NON_RESIDUE), aa);
  uint256 c1 = fq_sub(fq_sub(fq_mul(fq_add(x[0], x[1]), fq_add(y[0], y[1])), aa), bb);
  r[0] = c0;
  r[1] = c1;
}

inline void fq2_muc(const uint256 x[2], const uint256 &c, uint256 r[2]) {
  uint256 a = fq_mul(x[0], c);
  uint256 b = fq_mul(x[1], c);
  r[0] = a;
  r[1] = b;
}

inline void fq2_inv(const uint256 x[2], uint256 r[2]) {
  uint256 t = fq_inv(fq_sub(fq_mul(x[0], x[0]), fq_mul(fq_mul(x[1], x[1]), FQ_NON_RESIDUE)));
  r[0] = fq_mul(x[0], t);
  r[1] = fq_neg(fq_mul(x[1], t));
}

inline void fq2_squr(const uint256 x[2], uint256 r[2]) {
  uint256 a = fq_mul(x[0], x[1]);
  uint256 b = fq_mul(fq_add(fq_mul(x[1], FQ_NON_RESIDUE), x[0]), fq_add(x[0], x[1]));
  uint256 c = fq_sub(fq_sub(b, a), fq_mul(a, FQ_NON_RESIDUE));
  uint256 d = fq_add(a, a);
  r[0] = c;
  r[1] = d;
}

// G2_COEFF_B0 = mont_encode(0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5)
// G2_COEFF_B1 = mont_encode(0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2)
#define G2_COEFF_B0 "0x2514c6324384a86d26b7edf049755260020b1b273633535d3bf938e377b802a8"
#define G2_COEFF_B1 "0x0141b9ce4a688d4dd749d0dd22ac00aa65f0b37d93ce0d3e38e7ecccd1dcff67"
constexpr uint256 G2_COEFF_B[2] = {h256(G2_COEFF_B0), h256(G2_COEFF_B1)};

#define G2_ONE_00 "0x19573841af96503bfbb8264797811adfdceb1935497b01728e83b5d102bc2026"
#define G2_ONE_01 "0x14fef0833aea7b6b09e950fc52a02f866043dd5a5802d8c4afb4737da84c6140"
#define G2_ONE_10 "0x28fd7eebae9e4206ff9e1a62231b7dfefe7fd297f59e9b78619dfa9d886be9f6"
#define G2_ONE_11 "0x0da4a0e693fd648255f935be33351076dc57f922327d3cbb64095b56c71856ee"
constexpr uint256 G2_ONE[3][2] = {{h256(G2_ONE_00), h256(G2_ONE_01)}, {h256(G2_ONE_10), h256(G2_ONE_11)}, {FQ_ONE, 0}};
constexpr uint256 G2_ZERO[3][2] = {{0, 0}, {FQ2_ONE[0], FQ2_ONE[1]}, {0, 0}};

void g2_from_affine(const uint256 x[2][2], uint256 r[3][2]) {
  r[0][0] = x[0][0];
  r[0][1] = x[0][1];
  r[1][0] = x[1][0];
  r[1][1] = x[1][1];
  r[2][0] = FQ2_ONE[0];
  r[2][1] = FQ2_ONE[1];
}

void g2_from_jacobian(const uint256 x[3][2], uint256 r[2][2]) {
  if (x[2][0] == 0 && x[2][1] == 0) {
    assert(0);
  } else if (x[2][0] == 1 && x[2][1] == 0) {
    r[0][0] = x[0][0];
    r[0][1] = x[0][1];
    r[1][0] = x[1][0];
    r[1][1] = x[1][1];
  } else {
    uint256 zinv[2] = {};
    fq2_inv(x[2], zinv);

    uint256 zinv_squared[2] = {};
    fq2_squr(zinv, zinv_squared);

    uint256 zinv_squared_mul_zinv[2] = {};
    fq2_mul(zinv_squared, zinv, zinv_squared_mul_zinv);

    fq2_mul(x[0], zinv_squared, r[0]);
    fq2_mul(x[1], zinv_squared_mul_zinv, r[1]);
  }
}

void g2_double(const uint256 x[3][2], uint256 r[3][2]) {
  // A = X * X
  uint256 a[2];
  fq2_squr(x[0], a);
  // B = Y * Y
  uint256 b[2];
  fq2_squr(x[1], b);
  // C = B * B
  uint256 c[2];
  fq2_squr(b, c);
  // D = (X + B) ** 2 - A - C
  uint256 d[2];
  uint256 z[2];
  fq2_add(x[0], b, z);
  fq2_squr(z, d);
  fq2_sub(d, a, z);
  fq2_sub(z, c, d);
  // D = D + D
  fq2_add(d, d, d);
  // E = A + A + A
  uint256 e[2];
  fq2_add(a, a, z);
  fq2_add(z, a, e);
  // F = E * E
  uint256 f[2];
  fq2_squr(e, f);
  // X3 = F - (D + D)
  uint256 x3[2];
  fq2_add(d, d, z);
  fq2_sub(f, z, x3);
  // C8 = C + C + C + C + C + C + C + C
  uint256 c8[2];
  fq2_add(c, c, c8);
  fq2_add(c8, c8, c8);
  fq2_add(c8, c8, c8);
  // YX = Y * Z
  uint256 yz[2];
  fq2_mul(x[1], x[2], yz);
  // R[0] = X3
  r[0][0] = x3[0];
  r[0][1] = x3[1];
  // R[1] = E * (D - X3) - C8
  fq2_sub(d, x3, z);
  fq2_mul(e, z, z);
  fq2_sub(z, c8, r[1]);
  // R[2] = YZ + YZ
  fq2_add(yz, yz, r[2]);
}

void g2_add(const uint256 x[3][2], const uint256 y[3][2], uint256 r[3][2]) {
  if (x[2][0] == 0 && x[2][1] == 0) {
    r[0][0] = y[0][0];
    r[0][1] = y[0][1];
    r[1][0] = y[1][0];
    r[1][1] = y[1][1];
    r[2][0] = y[2][0];
    r[2][1] = y[2][1];
    return;
  }
  if (y[2][0] == 0 && y[2][1] == 0) {
    r[0][0] = x[0][0];
    r[0][1] = x[0][1];
    r[1][0] = x[1][0];
    r[1][1] = x[1][1];
    r[2][0] = x[2][0];
    r[2][1] = x[2][1];
    return;
  }

  uint256 z1_squared[2];
  fq2_squr(x[2], z1_squared);

  uint256 z2_squared[2];
  fq2_squr(y[2], z2_squared);

  uint256 u1[2];
  fq2_mul(x[0], z2_squared, u1);

  uint256 u2[2];
  fq2_mul(y[0], z1_squared, u2);

  uint256 z1_cubed[2];
  fq2_mul(x[2], z1_squared, z1_cubed);

  uint256 z2_cubed[2];
  fq2_mul(y[2], z2_squared, z2_cubed);

  uint256 s1[2];
  fq2_mul(x[1], z2_cubed, s1);

  uint256 s2[2];
  fq2_mul(y[1], z1_cubed, s2);

  if (arreq(u1, u2, 2) && arreq(s1, s2, 2)) {
    g2_double(x, r);
    return;
  }

  uint256 h[2];
  fq2_sub(u2, u1, h);

  uint256 s2_minus_s1[2];
  fq2_sub(s2, s1, s2_minus_s1);

  uint256 i[2];
  uint256 z[2];
  fq2_add(h, h, z);
  fq2_squr(z, i);

  uint256 j[2];
  fq2_mul(h, i, j);

  uint256 R[2];
  fq2_add(s2_minus_s1, s2_minus_s1, R);

  uint256 v[2];
  fq2_mul(u1, i, v);

  uint256 s1_j[2];
  fq2_mul(s1, j, s1_j);

  uint256 x3[2];
  fq2_squr(R, x3);
  fq2_sub(x3, j, x3);
  fq2_add(v, v, z);
  fq2_sub(x3, z, x3);

  r[0][0] = x3[0];
  r[0][1] = x3[1];
  uint256 v_minus_x3[2];
  fq2_sub(v, x3, v_minus_x3);
  fq2_mul(R, v_minus_x3, r[2]);
  fq2_add(s1_j, s1_j, z);
  fq2_sub(r[2], z, r[1]);

  fq2_add(x[2], y[2], r[2]);
  fq2_squr(r[2], z);
  fq2_sub(z, z1_squared, r[2]);
  fq2_sub(r[2], z2_squared, z);
  fq2_mul(z, h, r[2]);
}

void g2_mul(const uint256 x[3][2], const uint256 &y, uint256 r[3][2]) {
  uint256 a[3][2] = {{0, 0}, {FQ2_ONE[0], FQ2_ONE[1]}, {0, 0}}; // G2_ZERO
  uint256 b[3][2];

  bool found_one = 0;
  for (int i = 255; i > -1; i--) {
    if (found_one) {
      g2_double(a, b);
      a[0][0] = b[0][0];
      a[0][1] = b[0][1];
      a[1][0] = b[1][0];
      a[1][1] = b[1][1];
      a[2][0] = b[2][0];
      a[2][1] = b[2][1];
    }
    if (y & (uint256{1} << i)) {
      found_one = 1;
      g2_add(a, x, b);
      a[0][0] = b[0][0];
      a[0][1] = b[0][1];
      a[1][0] = b[1][0];
      a[1][1] = b[1][1];
      a[2][0] = b[2][0];
      a[2][1] = b[2][1];
    }
  }
  r[0][0] = a[0][0];
  r[0][1] = a[0][1];
  r[1][0] = a[1][0];
  r[1][1] = a[1][1];
  r[2][0] = a[2][0];
  r[2][1] = a[2][1];
}

} // namespace bn128

#endif /* BN128_H_ */
