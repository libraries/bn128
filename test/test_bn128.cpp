#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

bool eq12(const uint256 x[12], const uint256 y[12]) {
  for (int i = 0; i < 12; i++) {
    if (x[i] != y[i]) {
      return false;
    }
  }
  return true;
}

int main() {
  // Assert FQ(2) + FQ(2) == FQ(4)
  if (fq_add(2, 2) != 4) {
    return 1;
  }

  // Assert FQ(2) * FQ(2) == FQ(4)
  if (fq_mul(2, 2) != 4) {
    return 1;
  }

  // Assert FQ(2) / FQ(7) + FQ(9) / FQ(7) == FQ(11) / FQ(7)
  if (fq_add(fq_div(2, 7), fq_div(9, 7)) != fq_div(11, 7)) {
    return 1;
  }

  // Assert FQ(2) * FQ(7) + FQ(9) * FQ(7) == FQ(11) * FQ(7)
  if (fq_add(fq_mul(2, 7), fq_mul(9, 7)) != fq_mul(11, 7)) {
    return 1;
  }

  // Assert FQ(9) ** FIELD_MODULUS == FQ(9)
  if (fq_pow(9, FIELD_MODULUS) != 9) {
    return 1;
  }

  // Assert -1 == FIELD_MODULUS - 1
  if (fq_neg(1) != FIELD_MODULUS - 1) {
    return 1;
  }

  uint256 fq2_x[2] = {1, 0};
  uint256 fq2_f[2] = {1, 2};
  uint256 fq2_fpx[2] = {2, 2};
  uint256 fq2_one[2] = {1, 0};
  uint256 fq2_tmp[8][2] = {};

  // Assert x + f == fpx
  fq2_add(fq2_x, fq2_f, fq2_tmp[0]);
  if (!eq2(fq2_tmp[0], fq2_fpx)) {
    return 2;
  }

  // Assert f / f == one
  fq2_div(fq2_f, fq2_f, fq2_tmp[0]);
  if (!eq2(fq2_tmp[0], fq2_one)) {
    return 2;
  }

  // Assert one / f + x / f == (one + x) / f
  fq2_div(fq2_one, fq2_f, fq2_tmp[0]);
  fq2_div(fq2_x, fq2_f, fq2_tmp[1]);
  fq2_add(fq2_tmp[0], fq2_tmp[1], fq2_tmp[0]);
  fq2_add(fq2_one, fq2_x, fq2_tmp[1]);
  fq2_div(fq2_tmp[1], fq2_f, fq2_tmp[1]);
  if (!eq2(fq2_tmp[0], fq2_tmp[1])) {
    return 2;
  }

  // Assert one * f + x * f == (one + x) * f
  fq2_mul(fq2_one, fq2_f, fq2_tmp[0]);
  fq2_mul(fq2_x, fq2_f, fq2_tmp[1]);
  fq2_add(fq2_tmp[0], fq2_tmp[1], fq2_tmp[0]);
  fq2_add(fq2_one, fq2_x, fq2_tmp[1]);
  fq2_mul(fq2_tmp[1], fq2_f, fq2_tmp[1]);
  if (!eq2(fq2_tmp[0], fq2_tmp[1])) {
    return 2;
  }

  // Assert x ** (FIELD_MODULUS ** 2 - 1) == one
  fq2_pow(fq2_x, fq_sub(fq_pow(FIELD_MODULUS, 2), 1), fq2_tmp[0]);
  if (!eq2(fq2_tmp[0], fq2_one)) {
    return 2;
  }

  uint256 fq12_x[12] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  uint256 fq12_f[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  uint256 fq12_fpx[12] = {2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  uint256 fq12_one[12] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  uint256 fq12_tmp[8][12] = {};

  // Assert x + f == fpx
  fq12_add(fq12_x, fq12_f, fq12_tmp[0]);
  if (!eq12(fq12_tmp[0], fq12_fpx)) {
    return 3;
  }

  // Assert f / f == one
  fq12_div(fq12_f, fq12_f, fq12_tmp[0]);
  if (!eq12(fq12_tmp[0], fq12_one)) {
    return 3;
  }

  // Assert one / f + x / f == (one + x) / f
  fq12_div(fq12_one, fq12_f, fq12_tmp[0]);
  fq12_div(fq12_x, fq12_f, fq12_tmp[1]);
  fq12_add(fq12_tmp[0], fq12_tmp[1], fq12_tmp[2]);
  fq12_add(fq12_one, fq12_x, fq12_tmp[0]);
  fq12_div(fq12_tmp[0], fq12_f, fq12_tmp[1]);
  if (!eq12(fq12_tmp[1], fq12_tmp[2])) {
    return 3;
  }

  // Assert one * f + x * f == (one + x) * f
  fq12_mul(fq12_one, fq12_f, fq12_tmp[0]);
  fq12_mul(fq12_x, fq12_f, fq12_tmp[1]);
  fq12_add(fq12_tmp[0], fq12_tmp[1], fq12_tmp[2]);
  fq12_add(fq12_one, fq12_x, fq12_tmp[0]);
  fq12_mul(fq12_tmp[0], fq12_f, fq12_tmp[1]);
  if (!eq12(fq12_tmp[1], fq12_tmp[2])) {
    return 3;
  }

  // Assert x ** (FIELD_MODULUS ** 12 - 1) == one
  fq12_pow(fq12_x, fq_sub(fq_pow(FIELD_MODULUS, 12), 1), fq12_tmp[0]);
  if (!eq12(fq12_tmp[0], fq12_one)) {
    return 3;
  }

  // Assert FQ2([3, 0]) / FQ2([9, 1]) == B2
  fq2_tmp[0][0] = 3;
  fq2_tmp[0][1] = 0;
  fq2_tmp[1][0] = 9;
  fq2_tmp[1][1] = 1;
  fq2_div(fq2_tmp[0], fq2_tmp[1], fq2_tmp[0]);
  if (!eq2(fq2_tmp[0], B2)) {
    return 4;
  }

  // Assert G1 is on curve
  if (!g1::is_on_curve(G1)) {
    return 4;
  }

  // Assert G2 is on curve
  if (!g2::is_on_curve(G2)) {
    return 4;
  }

  // Assert Double(G1) == {0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd3, 0x15ed738c0e0a7c92e7845f96b2ae9c0a68a6a449e3538fc7ff3ebf7a5a18a2c4}
  g1::doubl2(G1, fq2_tmp[0]);
  fq2_tmp[1][0] = intx::from_string<uint256>("0x030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd3");
  fq2_tmp[1][1] = intx::from_string<uint256>("0x15ed738c0e0a7c92e7845f96b2ae9c0a68a6a449e3538fc7ff3ebf7a5a18a2c4");
  if (!eq2(fq2_tmp[0], fq2_tmp[1])) {
    return 4;
  }

  // Assert add(add(double(G1), G1), G1) == double(double(G1))
  g1::doubl2(G1, fq2_tmp[1]);
  g1::add(fq2_tmp[1], G1, fq2_tmp[0]);
  g1::add(fq2_tmp[0], G1, fq2_tmp[1]);
  g1::doubl2(G1, fq2_tmp[0]);
  g1::doubl2(fq2_tmp[0], fq2_tmp[2]);
  if (fq2_tmp[1][0] != fq2_tmp[2][0] || fq2_tmp[1][1] != fq2_tmp[2][1]) {
      return 4;
  }

  // Assert double(G1) != G1
  g1::doubl2(G1, fq2_tmp[0]);
  if (fq2_tmp[0][0] == G1[0] && fq2_tmp[0][1] == G1[1]) {
      return 4;
  }

  // Assert add(mul(G1, 9), mul(G1, 5)) == add(mul(G1, 12), mul(G1, 2))
  g1::mul(G1, 9, fq2_tmp[1]);
  g1::mul(G1, 5, fq2_tmp[2]);
  g1::add(fq2_tmp[1], fq2_tmp[2], fq2_tmp[0]);
  g1::mul(G1, 12, fq2_tmp[1]);
  g1::mul(G1, 2, fq2_tmp[2]);
  g1::add(fq2_tmp[1], fq2_tmp[2], fq2_tmp[3]);
  if (fq2_tmp[0][0] != fq2_tmp[3][0] || fq2_tmp[0][1] != fq2_tmp[3][1]) {
      return 4;
  }

  // Assert mul(G1, CURVE_ORDER) == inf
  g1::mul(G1, CURVE_ORDER, fq2_tmp[0]);
  if (!eq2(fq2_tmp[0], g1::INF)) {
      return 4;
  }

  return 0;
}
