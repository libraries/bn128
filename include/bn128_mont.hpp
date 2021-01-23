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

// The prime modulus of the field.
#define FIELD_MODULUS_HEX "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47"
constexpr uint256 FIELD_MODULUS = intx::from_string<uint256>(FIELD_MODULUS_HEX);
// R = 2 ** 256
// R_SQUARD = R * R % FIELD_MODULUS
// R_CUBED = R_SQUARD * R % FIELD_MODULUS
#define R_SQUARD_HEX "0x06d89f71cab8351f47ab1eff0a417ff6b5e71911d44501fbf32cfc5b538afa89"
#define R_CUBED_HEX "0x20fd6e902d592544ef7f0b0c0ada0afb62f210e6a7283db6b1cd6dafda1530df"
constexpr uint256 R_SQUARD = intx::from_string<uint256>(R_SQUARD_HEX);
constexpr uint256 R_CUBED = intx::from_string<uint256>(R_CUBED_HEX);
// R_PRIME * R % FIELD_MODULUS == 1
#define R_PRIME_HEX "0x2e67157159e5c639cf63e9cfb74492d9eb2022850278edf8ed84884a014afa37"
constexpr uint256 R_PRIME = intx::from_string<uint256>(R_PRIME_HEX);
#define FIELD_MODULUS_PRIME_HEX "0xf57a22b791888c6bd8afcbd01833da809ede7d651eca6ac987d20782e4866389"
// FIELD_MODULUS_PRIME * (-FIELD_MODULUS) % R == 1
constexpr uint256 FIELD_MODULUS_PRIME = intx::from_string<uint256>(FIELD_MODULUS_PRIME_HEX);

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

inline uint256 fq_add(const uint256 &x, const uint256 &y) { return _addmod(x, y, FIELD_MODULUS); }

inline uint256 fq_sub(const uint256 &x, const uint256 &y) { return _submod(x, y, FIELD_MODULUS); }

inline uint256 fq_mul(const uint256 &x, const uint256 &y) { return REDC(_mulmod(x, y, FIELD_MODULUS)); }

inline uint256 fq_inv(const uint256 &x) { return fq_mul(_invmod(x, FIELD_MODULUS), R_CUBED); }

inline uint256 fq_neg(const uint256 &x) { return _negmod(x, FIELD_MODULUS); }

// FQ_NON_RESIDUE = mont_encode(FIELD_MODULUS - 1);
#define FQ_NON_RESIDUE_HEX "0x2259d6b14729c0fa51e1a247090812318d087f6872aabf4f68c3488912edefaa"
constexpr uint256 FQ_NON_RESIDUE = intx::from_string<uint256>(FQ_NON_RESIDUE_HEX);

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

void fq2_neg(const uint256 x[2], uint256 r[2]) {
  uint256 a = fq_neg(x[0]);
  uint256 b = fq_neg(x[1]);
  r[0] = a;
  r[1] = b;
}

void fq2_mul(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 aa = fq_mul(x[0], y[0]);
  uint256 bb = fq_mul(x[1], y[1]);
  uint256 c0 = fq_add(fq_mul(bb, FQ_NON_RESIDUE), aa);
  uint256 c1 = fq_sub(fq_sub(fq_mul(fq_add(x[0], x[1]), fq_add(y[0], y[1])), aa), bb);
  r[0] = c0;
  r[1] = c1;
}

void fq2_muc(const uint256 x[2], const uint256 &c, uint256 r[2]) {
  uint256 a = fq_mul(x[0], c);
  uint256 b = fq_mul(x[1], c);
  r[0] = a;
  r[1] = b;
}

void fq2_inv(const uint256 x[2], uint256 r[2]) {
  uint256 t = fq_inv(fq_sub(fq_mul(x[0], x[0]), fq_mul(fq_mul(x[1], x[1]), FQ_NON_RESIDUE)));
  r[0] = fq_mul(x[0], t);
  r[1] = fq_neg(fq_mul(x[1], t));
}

} // namespace bn128

#endif /* BN128_H_ */
