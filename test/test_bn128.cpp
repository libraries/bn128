#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

int test_fq() {
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

  return 0;
}

int test_fq2() {
  uint256 x[2] = {1, 0};
  uint256 f[2] = {1, 2};
  uint256 fpx[2] = {2, 2};
  uint256 one[2] = {1, 0};
  uint256 tmp[8][2] = {};

  // Assert x + f == fpx
  fq2_add(x, f, tmp[0]);
  if (!eq2(tmp[0], fpx)) {
    return 1;
  }

  // Assert f / f == one
  fq2_div(f, f, tmp[0]);
  if (!eq2(tmp[0], one)) {
    return 1;
  }

  // Assert one / f + x / f == (one + x) / f
  fq2_div(one, f, tmp[0]);
  fq2_div(x, f, tmp[1]);
  fq2_add(tmp[0], tmp[1], tmp[0]);
  fq2_add(one, x, tmp[1]);
  fq2_div(tmp[1], f, tmp[1]);
  if (!eq2(tmp[0], tmp[1])) {
    return 1;
  }

  // Assert one * f + x * f == (one + x) * f
  fq2_mul(one, f, tmp[0]);
  fq2_mul(x, f, tmp[1]);
  fq2_add(tmp[0], tmp[1], tmp[0]);
  fq2_add(one, x, tmp[1]);
  fq2_mul(tmp[1], f, tmp[1]);
  if (!eq2(tmp[0], tmp[1])) {
    return 1;
  }

  // Assert x ** (FIELD_MODULUS ** 2 - 1) == one
  fq2_pow(x, fq_sub(fq_pow(FIELD_MODULUS, 2), 1), tmp[0]);
  if (!eq2(tmp[0], one)) {
    return 1;
  }

  // Assert FQ2([3, 0]) / FQ2([9, 1]) == B2
  tmp[0][0] = 3;
  tmp[0][1] = 0;
  tmp[1][0] = 9;
  tmp[1][1] = 1;
  fq2_div(tmp[0], tmp[1], tmp[0]);
  if (!eq2(tmp[0], B2)) {
    return 1;
  }

  return 0;
}

int test_fq12() {
  uint256 x[12] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  uint256 f[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  uint256 fpx[12] = {2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  uint256 one[12] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  uint256 tmp[8][12] = {};

  // Assert x + f == fpx
  fq12_add(x, f, tmp[0]);
  if (!eq12(tmp[0], fpx)) {
    return 1;
  }

  // Assert f / f == one
  fq12_div(f, f, tmp[0]);
  if (!eq12(tmp[0], one)) {
    return 1;
  }

  // Assert one / f + x / f == (one + x) / f
  fq12_div(one, f, tmp[0]);
  fq12_div(x, f, tmp[1]);
  fq12_add(tmp[0], tmp[1], tmp[2]);
  fq12_add(one, x, tmp[0]);
  fq12_div(tmp[0], f, tmp[1]);
  if (!eq12(tmp[1], tmp[2])) {
    return 1;
  }

  // Assert one * f + x * f == (one + x) * f
  fq12_mul(one, f, tmp[0]);
  fq12_mul(x, f, tmp[1]);
  fq12_add(tmp[0], tmp[1], tmp[2]);
  fq12_add(one, x, tmp[0]);
  fq12_mul(tmp[0], f, tmp[1]);
  if (!eq12(tmp[1], tmp[2])) {
    return 1;
  }

  // Assert x ** (FIELD_MODULUS ** 12 - 1) == one
  fq12_pow(x, fq_sub(fq_pow(FIELD_MODULUS, 12), 1), tmp[0]);
  if (!eq12(tmp[0], one)) {
    return 1;
  }

  return 0;
}

int test_g1() {
  uint256 tmp[4][2] = {};

  // Assert G1 is on curve
  if (!g1::is_on_curve(G1)) {
    return 1;
  }

  // Assert Double(G1)
  g1::doubl2(G1, tmp[0]);
  tmp[1][0] = intx::from_string<uint256>(
      "0x030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd3");
  tmp[1][1] = intx::from_string<uint256>(
      "0x15ed738c0e0a7c92e7845f96b2ae9c0a68a6a449e3538fc7ff3ebf7a5a18a2c4");
  if (!eq2(tmp[0], tmp[1])) {
    return 1;
  }

  // Assert add(add(double(G1), G1), G1) == double(double(G1))
  g1::doubl2(G1, tmp[1]);
  g1::add(tmp[1], G1, tmp[0]);
  g1::add(tmp[0], G1, tmp[1]);
  g1::doubl2(G1, tmp[0]);
  g1::doubl2(tmp[0], tmp[2]);
  if (!eq2(tmp[1], tmp[2])) {
    return 1;
  }

  // Assert double(G1) != G1
  g1::doubl2(G1, tmp[0]);
  if (eq2(tmp[0], G1)) {
    return 1;
  }

  // Assert add(mul(G1, 9), mul(G1, 5)) == add(mul(G1, 12), mul(G1, 2))
  g1::mul(G1, 9, tmp[1]);
  g1::mul(G1, 5, tmp[2]);
  g1::add(tmp[1], tmp[2], tmp[0]);
  g1::mul(G1, 12, tmp[1]);
  g1::mul(G1, 2, tmp[2]);
  g1::add(tmp[1], tmp[2], tmp[3]);
  if (!eq2(tmp[0], tmp[3])) {
    return 1;
  }

  // Assert mul(G1, CURVE_ORDER) == inf
  g1::mul(G1, CURVE_ORDER, tmp[0]);
  if (!eq2(tmp[0], g1::INF)) {
    return 1;
  }

  return 0;
}

int test_linefunc() {
  uint256 one[2] = {};
  uint256 two[2] = {};
  uint256 trd[2] = {};
  cp2(G1, one);
  g1::doubl2(G1, two);
  g1::mul(G1, 3, trd);

  uint256 neg_one[2] = {};
  uint256 neg_two[2] = {};
  uint256 neg_trd[2] = {};
  g1::mul(G1, CURVE_ORDER - 1, neg_one);
  g1::mul(G1, CURVE_ORDER - 2, neg_two);
  g1::mul(G1, CURVE_ORDER - 3, neg_trd);

  // Assert linefunc(one, two, one) == FQ(0)
  if (linefunc(one, two, one) != 0) {
    return 1;
  }

  // Assert linefunc(one, two, two) == FQ(0)
  if (linefunc(one, two, two) != 0) {
    return 1;
  }

  // Assert linefunc(one, two, three) != FQ(0)
  if (linefunc(one, two, trd) == 0) {
    return 1;
  }

  // Assert linefunc(one, two, negthree) == FQ(0)
  if (linefunc(one, two, neg_trd) != 0) {
    return 1;
  }

  // Assert linefunc(one, negone, one) == FQ(0)
  if (linefunc(one, neg_one, one) != 0) {
    return 1;
  }

  // Assert linefunc(one, negone, negone) == FQ(0)
  if (linefunc(one, neg_one, neg_one) != 0) {
    return 1;
  }

  // Assert linefunc(one, negone, two) != FQ(0)
  if (linefunc(one, neg_one, two) == 0) {
    return 1;
  }

  // Assert linefunc(one, one, one) == FQ(0)
  if (linefunc(one, one, one) != 0) {
    return 1;
  }

  // Assert linefunc(one, one, two) != FQ(0)
  if (linefunc(one, one, two) == 0) {
    return 1;
  }

  // Assert linefunc(one, one, negtwo) == FQ(0)
  if (linefunc(one, one, neg_two) != 0) {
    return 1;
  }

  return 0;
}

int main() {
  if (test_fq())
    return 1;
  if (test_fq2())
    return 1;
  if (test_fq12())
    return 1;
  if (test_g1())
    return 1;

  uint256 fq2_tmp[8][2] = {};

  // Assert G2 is on curve
  if (!g2::is_on_curve(G2)) {
    return 1;
  }

  uint256 pt2_tmp[8][2][2] = {};

  // Assert add(add(double(G2), G2), G2) == double(double(G2))
  g2::doubl2(G2, pt2_tmp[0]);
  g2::add(pt2_tmp[0], G2, pt2_tmp[1]);
  g2::add(pt2_tmp[1], G2, pt2_tmp[0]);
  g2::doubl2(G2, pt2_tmp[2]);
  g2::doubl2(pt2_tmp[2], pt2_tmp[1]);
  if (!(eq2(pt2_tmp[0][0], pt2_tmp[1][0]) &&
        eq2(pt2_tmp[0][1], pt2_tmp[1][1]))) {
    return 1;
  }

  // Assert double(G2) != G2
  g2::doubl2(G2, pt2_tmp[0]);
  if (eq2(pt2_tmp[0][0], G2[0]) && eq2(pt2_tmp[0][1], G2[1])) {
    return 1;
  }

  // Assert add(mul(G2, 9), mul(G2, 5)) == add(mul(G2, 12), mul(G2, 2))
  g2::mul(G2, 9, pt2_tmp[2]);
  g2::mul(G2, 5, pt2_tmp[3]);
  g2::add(pt2_tmp[2], pt2_tmp[3], pt2_tmp[0]);
  g2::mul(G2, 12, pt2_tmp[2]);
  g2::mul(G2, 2, pt2_tmp[3]);
  g2::add(pt2_tmp[2], pt2_tmp[3], pt2_tmp[1]);
  if (!(eq2(pt2_tmp[0][0], pt2_tmp[1][0]) &&
        eq2(pt2_tmp[0][1], pt2_tmp[1][1]))) {
    return 1;
  }

  // Assert mul(G2, CURVE_ORDER) == inf
  g2::mul(G2, CURVE_ORDER, pt2_tmp[0]);
  if (!g2::is_inf(pt2_tmp[0])) {
    return 1;
  }

  // Assert mul(G2, 2 * FIELD_MODULUS - CURVE_ORDER) != inf
  fq2_tmp[0][0] = fq_sub(fq_mul(2, FIELD_MODULUS), CURVE_ORDER);
  pt2_tmp[0][0][0] = fq_mul(G2[0][0], fq2_tmp[0][0]);
  pt2_tmp[0][0][1] = fq_mul(G2[0][1], fq2_tmp[0][0]);
  pt2_tmp[0][1][0] = fq_mul(G2[1][0], fq2_tmp[0][0]);
  pt2_tmp[0][1][1] = fq_mul(G2[1][1], fq2_tmp[0][0]);
  if (g2::is_inf(pt2_tmp[0])) {
    return 1;
  }

  // Assert is_on_curve(mul(G2, 9))
  g2::mul(G2, 9, pt2_tmp[0]);
  if (!g2::is_on_curve(pt2_tmp[0])) {
    return 1;
  }

  uint256 ptc_tmp[8][2][12] = {};

  // Assert twist(G2) == G12;
  g2::twist(G2, ptc_tmp[0]);
  if (!(eq12(ptc_tmp[0][0], G12[0]) && eq12(ptc_tmp[0][1], G12[1]))) {
    return 1;
  }

  // Assert add(add(double(G12), G12), G12) == double(double(G12))
  g12::doubl2(G12, ptc_tmp[0]);
  g12::add(ptc_tmp[0], G12, ptc_tmp[1]);
  g12::add(ptc_tmp[1], G12, ptc_tmp[0]);
  g12::doubl2(G12, ptc_tmp[2]);
  g12::doubl2(ptc_tmp[2], ptc_tmp[1]);
  if (!(eq12(ptc_tmp[0][0], ptc_tmp[1][0]) &&
        eq12(ptc_tmp[0][1], ptc_tmp[1][1]))) {
    return 1;
  }

  // Assert double(G12) != G12
  g12::doubl2(G12, ptc_tmp[0]);
  if (eq12(ptc_tmp[0][0], G12[0]) && eq12(ptc_tmp[0][1], G12[1])) {
    return 1;
  }

  // Assert add(mul(G12, 9), mul(G12, 5)) == add(mul(G12, 12), mul(G12, 2))
  g12::mul(G12, 9, ptc_tmp[2]);
  g12::mul(G12, 5, ptc_tmp[3]);
  g12::add(ptc_tmp[2], ptc_tmp[3], ptc_tmp[0]);
  g12::mul(G12, 9, ptc_tmp[2]);
  g12::mul(G12, 5, ptc_tmp[3]);
  g12::add(ptc_tmp[2], ptc_tmp[3], ptc_tmp[1]);
  if (!(eq12(ptc_tmp[0][0], ptc_tmp[1][0]) &&
        eq12(ptc_tmp[0][1], ptc_tmp[1][1]))) {
    return 1;
  }

  // Assert is_on_curve(mul(G12, 9))
  g12::mul(G12, 9, ptc_tmp[0]);
  if (!g12::is_on_curve(ptc_tmp[0])) {
    return 1;
  }

  // Assert mul(G12, CURVE_ORDER) == inf
  g12::mul(G12, CURVE_ORDER, ptc_tmp[0]);
  if (!g12::is_inf(ptc_tmp[0])) {
    return 1;
  }

  // Taking from
  // https://github.com/xxuejie/benchmarking-wasm-ewasm-evm/blob/checkpoint/evmrace/ckbvm/bn256g2_test.cpp
  pt2_tmp[0][0][0] = intx::from_string<uint256>(
      "0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed");
  pt2_tmp[0][0][1] = intx::from_string<uint256>(
      "0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2");
  pt2_tmp[0][1][0] = intx::from_string<uint256>(
      "0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa");
  pt2_tmp[0][1][1] = intx::from_string<uint256>(
      "0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b");
  g2::mul(pt2_tmp[0], 0x2dddefa19, pt2_tmp[1]);
  pt2_tmp[2][0][0] = intx::from_string<uint256>(
      "0x23997083c2c4409869ee3546806a544c8c16bc46cc88598c4e1c853eb81d45b0");
  pt2_tmp[2][0][1] = intx::from_string<uint256>(
      "0x1142585a23028cbe57783f890d1a2f6837049fce43c9b3b5e8e14c40a43c617a");
  pt2_tmp[2][1][0] = intx::from_string<uint256>(
      "0x215a23c8a96e1ca11d52cf6e2d6ada4ed01ee7e09b06dbc7f3315e7e6e73b919");
  pt2_tmp[2][1][1] = intx::from_string<uint256>(
      "0x0edac9f3a977530e28d4a385e614bcb7a8f9c3c3cb65707c1b90b5ea86174512");
  if (!(eq2(pt2_tmp[1][0], pt2_tmp[2][0]) &&
        eq2(pt2_tmp[1][1], pt2_tmp[2][1]))) {
    return 1;
  }

  if (test_linefunc() != 0) {
    return 1;
  }

  return 0;
}
