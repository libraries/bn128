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
#define FIELD_MODULUS_HEXSTRING "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47"
constexpr uint256 FIELD_MODULUS = intx::from_string<uint256>(FIELD_MODULUS_HEXSTRING);

} // namespace bn128

#endif /* BN128_H_ */
