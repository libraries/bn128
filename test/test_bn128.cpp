#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

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

  // Assert Double(G1) != ?
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

  uint256 pt2_tmp[8][2][2] = {};

  // Assert add(add(double(G2), G2), G2) == double(double(G2))
  g2::doubl2(G2, pt2_tmp[0]);
  g2::add(pt2_tmp[0], G2, pt2_tmp[1]);
  g2::add(pt2_tmp[1], G2, pt2_tmp[0]);
  g2::doubl2(G2, pt2_tmp[2]);
  g2::doubl2(pt2_tmp[2], pt2_tmp[1]);
  if (!(eq2(pt2_tmp[0][0], pt2_tmp[1][0]) && eq2(pt2_tmp[0][1], pt2_tmp[1][1]))) {
    return 5;
  }

  // Assert double(G2) != G2
  g2::doubl2(G2, pt2_tmp[0]);
  if (eq2(pt2_tmp[0][0], G2[0]) && eq2(pt2_tmp[0][1], G2[1])) {
    return 5;
  }

  // Assert add(mul(G2, 9), mul(G2, 5)) == add(mul(G2, 12), mul(G2, 2))
  g2::mul(G2, 9, pt2_tmp[2]);
  g2::mul(G2, 5, pt2_tmp[3]);
  g2::add(pt2_tmp[2], pt2_tmp[3], pt2_tmp[0]);
  g2::mul(G2, 12, pt2_tmp[2]);
  g2::mul(G2, 2, pt2_tmp[3]);
  g2::add(pt2_tmp[2], pt2_tmp[3], pt2_tmp[1]);
  if (!(eq2(pt2_tmp[0][0], pt2_tmp[1][0]) && eq2(pt2_tmp[0][1], pt2_tmp[1][1]))) {
    return 5;
  }

  // Assert mul(G2, CURVE_ORDER) == inf
  g2::mul(G2, CURVE_ORDER, pt2_tmp[0]);
  if (!g2::is_inf(pt2_tmp[0])) {
    return 5;
  }

  // Assert mul(G2, 2 * FIELD_MODULUS - CURVE_ORDER) != inf
  fq2_tmp[0][0] = fq_sub(fq_mul(2, FIELD_MODULUS), CURVE_ORDER);
  pt2_tmp[0][0][0] = fq_mul(G2[0][0], fq2_tmp[0][0]);
  pt2_tmp[0][0][1] = fq_mul(G2[0][1], fq2_tmp[0][0]);
  pt2_tmp[0][1][0] = fq_mul(G2[1][0], fq2_tmp[0][0]);
  pt2_tmp[0][1][1] = fq_mul(G2[1][1], fq2_tmp[0][0]);
  if (g2::is_inf(pt2_tmp[0])) {
    return 5;
  }

  // Assert is_on_curve(mul(G2, 9))
  g2::mul(G2, 9, pt2_tmp[0]);
  if (!g2::is_on_curve(pt2_tmp[0])) {
    return 5;
  }

  uint256 ptc_tmp[8][2][12] = {};

  // Assert twist(G2) == G12;
  g2::twist(G2, ptc_tmp[0]);
  if (!(eq12(ptc_tmp[0][0], G12[0]) && eq12(ptc_tmp[0][1], G12[1]))) {
    return 6;
  }

  // Assert add(add(double(G12), G12), G12) == double(double(G12))
  g12::doubl2(G12, ptc_tmp[0]);
  g12::add(ptc_tmp[0], G12, ptc_tmp[1]);
  g12::add(ptc_tmp[1], G12, ptc_tmp[0]);
  g12::doubl2(G12, ptc_tmp[2]);
  g12::doubl2(ptc_tmp[2], ptc_tmp[1]);
  if (!(eq12(ptc_tmp[0][0], ptc_tmp[1][0]) && eq12(ptc_tmp[0][1], ptc_tmp[1][1]))) {
    return 6;
  }

  // Assert double(G12) != G12
  g12::doubl2(G12, ptc_tmp[0]);
  if (eq12(ptc_tmp[0][0], G12[0]) && eq12(ptc_tmp[0][1], G12[1])) {
    return 6;
  }

  // Assert add(mul(G12, 9), mul(G12, 5)) == add(mul(G12, 12), mul(G12, 2))
  g12::mul(G12, 9, ptc_tmp[2]);
  g12::mul(G12, 5, ptc_tmp[3]);
  g12::add(ptc_tmp[2], ptc_tmp[3], ptc_tmp[0]);
  g12::mul(G12, 9, ptc_tmp[2]);
  g12::mul(G12, 5, ptc_tmp[3]);
  g12::add(ptc_tmp[2], ptc_tmp[3], ptc_tmp[1]);
  if (!(eq12(ptc_tmp[0][0], ptc_tmp[1][0]) && eq12(ptc_tmp[0][1], ptc_tmp[1][1]))) {
    return 6;
  }

  // Assert is_on_curve(mul(G12, 9))
  g12::mul(G12, 9, ptc_tmp[0]);
  if (!g12::is_on_curve(ptc_tmp[0])) {
    return 6;
  }

  // Assert mul(G12, CURVE_ORDER) == inf
  g12::mul(G12, CURVE_ORDER, ptc_tmp[0]);
  if (!g12::is_inf(ptc_tmp[0])) {
    return 6;
  }

  // Taking from https://github.com/xxuejie/benchmarking-wasm-ewasm-evm/blob/checkpoint/evmrace/ckbvm/bn256g2_test.cpp
  pt2_tmp[0][0][0] = intx::from_string<uint256>("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed");
  pt2_tmp[0][0][1] = intx::from_string<uint256>("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2");
  pt2_tmp[0][1][0] = intx::from_string<uint256>("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa");
  pt2_tmp[0][1][1] = intx::from_string<uint256>("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b");
  g2::mul(pt2_tmp[0], 0x2dddefa19, pt2_tmp[1]);
  pt2_tmp[2][0][0] = intx::from_string<uint256>("0x23997083c2c4409869ee3546806a544c8c16bc46cc88598c4e1c853eb81d45b0");
  pt2_tmp[2][0][1] = intx::from_string<uint256>("0x1142585a23028cbe57783f890d1a2f6837049fce43c9b3b5e8e14c40a43c617a");
  pt2_tmp[2][1][0] = intx::from_string<uint256>("0x215a23c8a96e1ca11d52cf6e2d6ada4ed01ee7e09b06dbc7f3315e7e6e73b919");
  pt2_tmp[2][1][1] = intx::from_string<uint256>("0x0edac9f3a977530e28d4a385e614bcb7a8f9c3c3cb65707c1b90b5ea86174512");
  if (!(eq2(pt2_tmp[1][0], pt2_tmp[2][0]) && eq2(pt2_tmp[1][1], pt2_tmp[2][1]))) {
    return 7;
  }

  return 0;
}
