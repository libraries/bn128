#include <bn128_mont.hpp>
#include <intx/intx.hpp>

using namespace bn128;

int test_invmod() {
  constexpr uint256 x_case[8] = {
      intx::from_string<uint256>("0x0000000000000000000000000000000000000000000000000000000000000001"),
      intx::from_string<uint256>("0x00000000000000000000000012341234abcd0000000000000000000000000001"),
      intx::from_string<uint256>("0x0000000000000000000000000000000000000000000fffffffffffffffffffff"),
      intx::from_string<uint256>("0x0123456789a00000000000000000000000000000000000000000000000000001"),
      intx::from_string<uint256>("0x0000000000000000000000000000000000000000000000000000000000000007"),
      intx::from_string<uint256>("0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd07"),
      intx::from_string<uint256>("0x00644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47"),
      intx::from_string<uint256>("0x30644e72e131a029b82222228181585d97816a916871ca8d3c208c16d87cfd47"),
  };
  for (int i = 0; i < 8; i++) {
    uint256 a = x_case[i];
    uint256 r = _mulmod(a, _invmod(a, FIELD_MODULUS), FIELD_MODULUS);
    if (r != 1) {
      return 1;
    }
  }
  return 0;
}

int test_powmod() {
  uint256 r = _powmod(2, 10, FIELD_MODULUS);
  if (r != 1024) {
    return 1;
  }
  return 0;
}

int test_constexpr() {
  // R * R_PRIME % FIELD_MODULUS == 1
  if (uint512(1, 0) * uint512(R_PRIME) % uint512(FIELD_MODULUS) != 1) {
    return 1;
  }
  // -FIELD_MODULUS * FIELD_MODULUS_PRIME % R == 1
  if ((uint512(1, 0) - uint512(FIELD_MODULUS)) * uint512(FIELD_MODULUS_PRIME) % uint512(1, 0) != 1) {
    return 1;
  }
  return 0;
}

int test_mont_encode_decode() {
  // Taking from https://github.com/paritytech/bn
  uint256 a = mont_encode(1);
  if (a != intx::from_string<uint256>("0xe0a77c19a07df2f666ea36f7879462c0a78eb28f5c70b3dd35d438dc58f0d9d")) {
    return 1;
  }
  if (mont_decode(a) != 1) {
    return 1;
  }

  uint256 b = mont_encode(42);
  if (b != intx::from_string<uint256>("0x0903f860b6f71bd22a638bbbb1d55ed69dc595e76d525985dbc68430439c5c6e")) {
    return 1;
  }
  if (mont_decode(b) != 42) {
    return 1;
  }

  uint256 c = FIELD_MODULUS / 2 - 17;
  uint256 d = mont_encode(c);
  if (d != intx::from_string<uint256>("0x14707fbbcf072f27f529534d0bfd1a000a03b6d2f156954ed7d2e44ca56802cb")) {
    return 1;
  }
  if (mont_decode(d) != c) {
    return 1;
  }

  return 0;
}

int main() {
  if (test_invmod())
    return 1;
  if (test_powmod())
    return 1;
  if (test_constexpr())
    return 1;
  if (test_mont_encode_decode())
    return 1;
  return 0;
}
