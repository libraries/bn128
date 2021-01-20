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
  uint256 tmp[2][2] = {};

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
  uint256 tmp[3][12] = {};

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
  uint256 tmp[4][3] = {};

  // Assert G1 is on curve
  if (!g1::is_on_curve(G1)) {
    return 1;
  }

  // Assert add(add(double(G1), G1), G1) == double(double(G1))
  g1::doubl2(G1, tmp[1]);
  g1::add(tmp[1], G1, tmp[0]);
  g1::add(tmp[0], G1, tmp[1]);
  g1::doubl2(G1, tmp[0]);
  g1::doubl2(tmp[0], tmp[2]);
  if (!g1::eq(tmp[1], tmp[2])) {
    return 1;
  }

  // Assert double(G1) != G1
  g1::doubl2(G1, tmp[0]);
  if (g1::eq(tmp[0], G1)) {
    return 1;
  }

  // Assert add(mul(G1, 9), mul(G1, 5)) == add(mul(G1, 12), mul(G1, 2))
  g1::mul(G1, 9, tmp[1]);
  g1::mul(G1, 5, tmp[2]);
  g1::add(tmp[1], tmp[2], tmp[0]);
  g1::mul(G1, 12, tmp[1]);
  g1::mul(G1, 2, tmp[2]);
  g1::add(tmp[1], tmp[2], tmp[3]);
  if (!g1::eq(tmp[0], tmp[3])) {
    return 1;
  }

  // Assert mul(G1, CURVE_ORDER) == inf
  g1::mul(G1, CURVE_ORDER, tmp[0]);
  if (!g1::is_inf(tmp[0])) {
    return 1;
  }

  return 0;
}

int test_g2() {
  uint256 tmp[4][3][2] = {};

  // Assert G2 is on curve
  if (!g2::is_on_curve(G2)) {
    return 1;
  }

  // Assert add(add(double(G2), G2), G2) == double(double(G2))
  g2::doubl2(G2, tmp[0]);
  g2::add(tmp[0], G2, tmp[1]);
  g2::add(tmp[1], G2, tmp[0]);
  g2::doubl2(G2, tmp[2]);
  g2::doubl2(tmp[2], tmp[1]);
  if (!g2::eq(tmp[0], tmp[1])) {
    return 1;
  }

  // Assert double(G2) != G2
  g2::doubl2(G2, tmp[0]);
  if (g2::eq(tmp[0], G2)) {
    return 1;
  }

  // Assert add(mul(G2, 9), mul(G2, 5)) == add(mul(G2, 12), mul(G2, 2))
  g2::mul(G2, 9, tmp[2]);
  g2::mul(G2, 5, tmp[3]);
  g2::add(tmp[2], tmp[3], tmp[0]);
  g2::mul(G2, 12, tmp[2]);
  g2::mul(G2, 2, tmp[3]);
  g2::add(tmp[2], tmp[3], tmp[1]);
  if (!g2::eq(tmp[0], tmp[1])) {
    return 1;
  }

  // Assert mul(G2, CURVE_ORDER) == inf
  g2::mul(G2, CURVE_ORDER, tmp[0]);
  if (!g2::is_inf(tmp[0])) {
    return 1;
  }

  // Assert mul(G2, 2 * FIELD_MODULUS - CURVE_ORDER) != inf
  uint256 n = fq_sub(fq_mul(2, FIELD_MODULUS), CURVE_ORDER);
  g2::mul(G2, n, tmp[0]);
  if (g2::is_inf(tmp[0])) {
    return 1;
  }

  // Assert is_on_curve(mul(G2, 9))
  g2::mul(G2, 9, tmp[0]);
  if (!g2::is_on_curve(tmp[0])) {
    return 1;
  }

  return 0;
}

int test_g12() {
  uint256 tmp[4][3][12] = {};

  // Assert twist(G2) == G12;
  g2::twist(G2, tmp[0]);
  if (!(eq12(tmp[0][0], G12[0]) && eq12(tmp[0][1], G12[1])) && eq12(tmp[0][2], G12[2])) {
    return 1;
  }

  // Assert add(add(double(G12), G12), G12) == double(double(G12))
  g12::doubl2(G12, tmp[0]);
  g12::add(tmp[0], G12, tmp[1]);
  g12::add(tmp[1], G12, tmp[0]);
  g12::doubl2(G12, tmp[2]);
  g12::doubl2(tmp[2], tmp[1]);
  if (!g12::eq(tmp[0], tmp[1])) {
    return 1;
  }

  // Assert double(G12) != G12
  g12::doubl2(G12, tmp[0]);
  if (g12::eq(tmp[0], G12)) {
    return 1;
  }

  // Assert add(mul(G12, 9), mul(G12, 5)) == add(mul(G12, 12), mul(G12, 2))
  g12::mul(G12, 9, tmp[2]);
  g12::mul(G12, 5, tmp[3]);
  g12::add(tmp[2], tmp[3], tmp[0]);
  g12::mul(G12, 9, tmp[2]);
  g12::mul(G12, 5, tmp[3]);
  g12::add(tmp[2], tmp[3], tmp[1]);
  if (!g12::eq(tmp[0], tmp[1])) {
    return 1;
  }

  // Assert is_on_curve(mul(G12, 9))
  g12::mul(G12, 9, tmp[0]);
  if (!g12::is_on_curve(tmp[0])) {
    return 1;
  }

  // Assert mul(G12, CURVE_ORDER) == inf
  g12::mul(G12, CURVE_ORDER, tmp[0]);
  if (!g12::is_inf(tmp[0])) {
    return 1;
  }

  return 0;
}

int test_linefunc() {
  uint256 one[3] = {};
  uint256 two[3] = {};
  uint256 trd[3] = {};
  uint256 out[2] = {};
  cp(G1, one, 3);
  g1::doubl2(G1, two);
  g1::mul(G1, 3, trd);

  uint256 neg_one[3] = {};
  uint256 neg_two[3] = {};
  uint256 neg_trd[3] = {};
  g1::mul(G1, CURVE_ORDER - 1, neg_one);
  g1::mul(G1, CURVE_ORDER - 2, neg_two);
  g1::mul(G1, CURVE_ORDER - 3, neg_trd);

  // Assert linefunc(one, two, one) == FQ(0)
  linefunc(one, two, one, out);
  if (out[0]) {
    return 1;
  }

  // Assert linefunc(one, two, two) == FQ(0)
  linefunc(one, two, two, out);
  if (out[0]) {
    return 1;
  }

  // Assert linefunc(one, two, three) != FQ(0)
  linefunc(one, two, trd, out);
  if (!out[0]) {
    return 1;
  }

  // Assert linefunc(one, two, negthree) == FQ(0)
  linefunc(one, two, neg_trd, out);
  if (out[0]) {
    return 1;
  }

  // Assert linefunc(one, negone, one) == FQ(0)
  linefunc(one, neg_one, one, out);
  if (out[0]) {
    return 1;
  }

  // Assert linefunc(one, negone, negone) == FQ(0)
  linefunc(one, neg_one, neg_one, out);
  if (out[0]) {
    return 1;
  }

  // Assert linefunc(one, negone, two) != FQ(0)
  linefunc(one, neg_one, two, out);
  if (!out[0]) {
    return 1;
  }

  // Assert linefunc(one, one, one) == FQ(0)
  linefunc(one, one, one, out);
  if (out[0]) {
    return 1;
  }

  // Assert linefunc(one, one, two) != FQ(0)
  linefunc(one, one, two, out);
  if (!out[0]) {
    return 1;
  }

  // Assert linefunc(one, one, negtwo) == FQ(0)
  linefunc(one, one, neg_two, out);
  if (out[0]) {
    return 1;
  }

  return 0;
}

int test_pairing() {
  uint256 p1[12] = {};
  uint256 pn1[12] = {};
  uint256 tmp[12] = {};
  uint256 neg_g1[3] = {G1[0], fq_neg(G1[1]), G1[2]};

  // p1 = pairing(G2, G1)
  // pn1 = pairing(G2, neg(G1))
  // assert p1 * pn1 == FQ12.one()
  _pairing(G2, G1, p1);
  _pairing(G2, neg_g1, pn1);
  fq12_mul(p1, pn1, tmp);
  if (!eq12(tmp, FQ12_ONE)) {
    return 1;
  }

  // np1 = pairing(neg(G2), G1)
  // assert p1 * np1 == FQ12.one()
  // assert pn1 == np1
  uint256 neg_g2[3][2] = {};
  uint256 np1[12] = {};
  cp2(G2[0], neg_g2[0]);
  fq2_neg(G2[1], neg_g2[1]);
  cp2(G2[2], neg_g2[2]);
  _pairing(neg_g2, G1, np1);
  fq12_mul(p1, np1, tmp);
  if (!eq12(tmp, FQ12_ONE)) {
    return 1;
  }
  if (!eq12(pn1, np1)) {
    return 1;
  }

  // assert p1 ** curve_order == FQ12.one()
  fq12_pow(p1, CURVE_ORDER, tmp);
  if (!eq12(tmp, FQ12_ONE)) {
    return 1;
  }

  uint256 g1_mul_2[3] = {};
  uint256 p2[12] = {};
  g1::mul(G1, 2, g1_mul_2);
  _pairing(G2, g1_mul_2, p2);
  fq12_mul(p1, p1, tmp);
  if (!eq12(tmp, p2)) {
    return 1;
  }

  // # assert p1 != p2 and p1 != np1 and p2 != np1
  if ((eq12(p1, p2) || eq12(p1, np1) || eq12(p2, np1))) {
    return 1;
  }

  uint256 g2_mul_2[3][2] = {};
  uint256 po2[12] = {};
  g2::mul(G2, 2, g2_mul_2);
  _pairing(g2_mul_2, G1, po2);
  fq12_mul(p1, p1, tmp);
  if (!eq12(tmp, po2)) {
    return 1;
  }

  uint256 g2_mul_27[3][2] = {};
  uint256 g1_mul_31[3] = {};
  uint256 p3[12] = {};
  g2::mul(G2, 27, g2_mul_27);
  g1::mul(G1, 37, g1_mul_31);
  _pairing(g2_mul_27, g1_mul_31, p3);
  uint256 g1_mul_999[3] = {};
  uint256 po3[12] = {};
  g1::mul(G1, 999, g1_mul_999);
  _pairing(G2, g1_mul_999, po3);
  if (!eq12(p3, po3)) {
    return 1;
  }

  return 0;
}

int test_misc() {
  uint256 tmp[3][2][2] = {};
  // Taking from
  // https://github.com/xxuejie/benchmarking-wasm-ewasm-evm/blob/checkpoint/evmrace/ckbvm/bn256g2_test.cpp
  tmp[0][0][0] = intx::from_string<uint256>(
      "0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed");
  tmp[0][0][1] = intx::from_string<uint256>(
      "0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2");
  tmp[0][1][0] = intx::from_string<uint256>(
      "0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa");
  tmp[0][1][1] = intx::from_string<uint256>(
      "0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b");
  mul(tmp[0], 0x2dddefa19, tmp[2]);
  tmp[1][0][0] = intx::from_string<uint256>(
      "0x23997083c2c4409869ee3546806a544c8c16bc46cc88598c4e1c853eb81d45b0");
  tmp[1][0][1] = intx::from_string<uint256>(
      "0x1142585a23028cbe57783f890d1a2f6837049fce43c9b3b5e8e14c40a43c617a");
  tmp[1][1][0] = intx::from_string<uint256>(
      "0x215a23c8a96e1ca11d52cf6e2d6ada4ed01ee7e09b06dbc7f3315e7e6e73b919");
  tmp[1][1][1] = intx::from_string<uint256>(
      "0x0edac9f3a977530e28d4a385e614bcb7a8f9c3c3cb65707c1b90b5ea86174512");
  if (!(eq2(tmp[1][0], tmp[2][0]) && eq2(tmp[1][1], tmp[2][1]))) {
    return 1;
  }

  // Taking from
  // https://github.com/ewasm/ethereum-bn128.rs/blob/master/src/lib.rs#L318
  uint256 out[3][12];
  tmp[0][0][0] = intx::from_string<uint256>(
      "0x2eca0c7238bf16e83e7a1e6c5d49540685ff51380f309842a98561558019fc02");
  tmp[0][0][1] = intx::from_string<uint256>(
      "0x03d3260361bb8451de5ff5ecd17f010ff22f5c31cdf184e9020b06fa5997db84");
  tmp[1][0][1] = intx::from_string<uint256>(
      "0x1213d2149b006137fcfb23036606f848d638d576a120ca981b5b1a5f9300b3ee");
  tmp[1][0][0] = intx::from_string<uint256>(
      "0x2276cf730cf493cd95d64677bbb75fc42db72513a4c1e387b476d056f80aa75f");
  tmp[1][1][1] = intx::from_string<uint256>(
      "0x21ee6226d31426322afcda621464d0611d226783262e21bb3bc86b537e986237");
  tmp[1][1][0] = intx::from_string<uint256>(
      "0x096df1f82dff337dd5972e32a8ad43e28a78a96a823ef1cd4debe12b6552ea5f");
  pairing(tmp[1], tmp[0][0], out[0]);
  tmp[0][0][0] = intx::from_string<uint256>(
      "0x06967a1237ebfeca9aaae0d6d0bab8e28c198c5a339ef8a2407e31cdac516db9");
  tmp[0][0][1] = intx::from_string<uint256>(
      "0x22160fa257a5fd5b280642ff47b65eca77e626cb685c84fa6d3b6882a283ddd1");
  tmp[1][0][1] = intx::from_string<uint256>(
      "0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2");
  tmp[1][0][0] = intx::from_string<uint256>(
      "0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed");
  tmp[1][1][1] = intx::from_string<uint256>(
      "0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b");
  tmp[1][1][0] = intx::from_string<uint256>(
      "0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa");
  pairing(tmp[1], tmp[0][0], out[1]);
  fq12_mul(out[0], out[1], out[2]);
  if (!eq12(out[2], FQ12_ONE)) {
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
  if (test_g2())
    return 1;
  if (test_g12())
    return 1;
  if (test_linefunc())
    return 1;
  if (test_pairing())
    return 1;
  if (test_misc())
    return 1;

  return 0;
}
