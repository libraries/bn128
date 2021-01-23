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

int main() {
  if (test_invmod())
    return 1;
  return 0;
}
