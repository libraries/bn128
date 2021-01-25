#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

int test_invmod() {
  constexpr uint256 x_case[8] = {
      h256("0x0000000000000000000000000000000000000000000000000000000000000001"),
      h256("0x00000000000000000000000012341234abcd0000000000000000000000000001"),
      h256("0x0000000000000000000000000000000000000000000fffffffffffffffffffff"),
      h256("0x0123456789a00000000000000000000000000000000000000000000000000001"),
      h256("0x0000000000000000000000000000000000000000000000000000000000000007"),
      h256("0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd07"),
      h256("0x00644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47"),
      h256("0x30644e72e131a029b82222228181585d97816a916871ca8d3c208c16d87cfd47"),
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
  if (!arreq(tmp[0], fpx, 2)) {
    return 1;
  }

  // Assert f / f == one
  fq2_div(f, f, tmp[0]);
  if (!arreq(tmp[0], one, 2)) {
    return 1;
  }

  // Assert one / f + x / f == (one + x) / f
  fq2_div(one, f, tmp[0]);
  fq2_div(x, f, tmp[1]);
  fq2_add(tmp[0], tmp[1], tmp[0]);
  fq2_add(one, x, tmp[1]);
  fq2_div(tmp[1], f, tmp[1]);
  if (!arreq(tmp[0], tmp[1], 2)) {
    return 1;
  }

  // Assert one * f + x * f == (one + x) * f
  fq2_mul(one, f, tmp[0]);
  fq2_mul(x, f, tmp[1]);
  fq2_add(tmp[0], tmp[1], tmp[0]);
  fq2_add(one, x, tmp[1]);
  fq2_mul(tmp[1], f, tmp[1]);
  if (!arreq(tmp[0], tmp[1], 2)) {
    return 1;
  }

  // Assert x ** (FIELD_MODULUS ** 2 - 1) == one
  fq2_pow(x, fq_sub(fq_pow(FIELD_MODULUS, 2), 1), tmp[0]);
  if (!arreq(tmp[0], one, 2)) {
    return 1;
  }

  // Assert FQ2([3, 0]) / FQ2([9, 1]) == B2
  tmp[0][0] = 3;
  tmp[0][1] = 0;
  tmp[1][0] = 9;
  tmp[1][1] = 1;
  fq2_div(tmp[0], tmp[1], tmp[0]);
  if (!arreq(tmp[0], B2, 2)) {
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
  if (!arreq(tmp[0], fpx, 12)) {
    return 1;
  }

  // Assert f / f == one
  fq12_div(f, f, tmp[0]);
  if (!arreq(tmp[0], one, 12)) {
    return 1;
  }

  // Assert one / f + x / f == (one + x) / f
  fq12_div(one, f, tmp[0]);
  fq12_div(x, f, tmp[1]);
  fq12_add(tmp[0], tmp[1], tmp[2]);
  fq12_add(one, x, tmp[0]);
  fq12_div(tmp[0], f, tmp[1]);
  if (!arreq(tmp[1], tmp[2], 12)) {
    return 1;
  }

  // Assert one * f + x * f == (one + x) * f
  fq12_mul(one, f, tmp[0]);
  fq12_mul(x, f, tmp[1]);
  fq12_add(tmp[0], tmp[1], tmp[2]);
  fq12_add(one, x, tmp[0]);
  fq12_mul(tmp[0], f, tmp[1]);
  if (!arreq(tmp[1], tmp[2], 12)) {
    return 1;
  }

  // Assert x ** (FIELD_MODULUS ** 12 - 1) == one
  fq12_pow(x, fq_sub(fq_pow(FIELD_MODULUS, 12), 1), tmp[0]);
  if (!arreq(tmp[0], one, 12)) {
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
  if (!(arreq(tmp[0][0], G12[0], 12) && arreq(tmp[0][1], G12[1], 2)) && arreq(tmp[0][2], G12[2], 2)) {
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
  arrcp(G1, one, 3);
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
  if (!arreq(tmp, FQ12_ONE, 12)) {
    return 1;
  }

  // np1 = pairing(neg(G2), G1)
  // assert p1 * np1 == FQ12.one()
  // assert pn1 == np1
  uint256 neg_g2[3][2] = {};
  uint256 np1[12] = {};
  arrcp(G2[0], neg_g2[0], 2);
  fq2_neg(G2[1], neg_g2[1]);
  arrcp(G2[2], neg_g2[2], 2);
  _pairing(neg_g2, G1, np1);
  fq12_mul(p1, np1, tmp);
  if (!arreq(tmp, FQ12_ONE, 12)) {
    return 1;
  }
  if (!arreq(pn1, np1, 12)) {
    return 1;
  }

  // assert p1 ** curve_order == FQ12.one()
  fq12_pow(p1, CURVE_ORDER, tmp);
  if (!arreq(tmp, FQ12_ONE, 12)) {
    return 1;
  }

  uint256 g1_mul_2[3] = {};
  uint256 p2[12] = {};
  g1::mul(G1, 2, g1_mul_2);
  _pairing(G2, g1_mul_2, p2);
  fq12_mul(p1, p1, tmp);
  if (!arreq(tmp, p2, 12)) {
    return 1;
  }

  // # assert p1 != p2 and p1 != np1 and p2 != np1
  if ((arreq(p1, p2, 12) || arreq(p1, np1, 12) || arreq(p2, np1, 12))) {
    return 1;
  }

  uint256 g2_mul_2[3][2] = {};
  uint256 po2[12] = {};
  g2::mul(G2, 2, g2_mul_2);
  _pairing(g2_mul_2, G1, po2);
  fq12_mul(p1, p1, tmp);
  if (!arreq(tmp, po2, 12)) {
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
  if (!arreq(p3, po3, 12)) {
    return 1;
  }

  return 0;
}

int test_misc() {
  // Taking from
  // https://github.com/xxuejie/benchmarking-wasm-ewasm-evm/blob/checkpoint/evmrace/ckbvm/bn256g2_test.cpp
  uint256 a[2][2];
  a[0][0] = h256("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed");
  a[0][1] = h256("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2");
  a[1][0] = h256("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa");
  a[1][1] = h256("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b");
  uint256 b[3][2];
  arrcp(a[0], b[0], 2);
  arrcp(a[1], b[1], 2);
  arrcp(FQ2_ONE, b[2], 2);
  uint256 c[3][2];
  g2::mul(b, 0x2dddefa19, c);
  uint256 d[2][2];
  g2::from_jacobian(c, d);
  if (d[0][0] != h256("0x23997083c2c4409869ee3546806a544c8c16bc46cc88598c4e1c853eb81d45b0")) {
    return 1;
  }
  if (d[0][1] != h256("0x1142585a23028cbe57783f890d1a2f6837049fce43c9b3b5e8e14c40a43c617a")) {
    return 1;
  }
  if (d[1][0] != h256("0x215a23c8a96e1ca11d52cf6e2d6ada4ed01ee7e09b06dbc7f3315e7e6e73b919")) {
    return 1;
  }
  if (d[1][1] != h256("0x0edac9f3a977530e28d4a385e614bcb7a8f9c3c3cb65707c1b90b5ea86174512")) {
    return 1;
  }

  uint256 tmp[3][2][2] = {};
  // Taking from
  // https://github.com/ewasm/ethereum-bn128.rs/blob/master/src/lib.rs#L318
  uint256 out[3][12];
  tmp[0][0][0] = h256("0x2eca0c7238bf16e83e7a1e6c5d49540685ff51380f309842a98561558019fc02");
  tmp[0][0][1] = h256("0x03d3260361bb8451de5ff5ecd17f010ff22f5c31cdf184e9020b06fa5997db84");
  tmp[1][0][1] = h256("0x1213d2149b006137fcfb23036606f848d638d576a120ca981b5b1a5f9300b3ee");
  tmp[1][0][0] = h256("0x2276cf730cf493cd95d64677bbb75fc42db72513a4c1e387b476d056f80aa75f");
  tmp[1][1][1] = h256("0x21ee6226d31426322afcda621464d0611d226783262e21bb3bc86b537e986237");
  tmp[1][1][0] = h256("0x096df1f82dff337dd5972e32a8ad43e28a78a96a823ef1cd4debe12b6552ea5f");
  alt_bn128_pairing(tmp[1], tmp[0][0], out[0]);
  tmp[0][0][0] = h256("0x06967a1237ebfeca9aaae0d6d0bab8e28c198c5a339ef8a2407e31cdac516db9");
  tmp[0][0][1] = h256("0x22160fa257a5fd5b280642ff47b65eca77e626cb685c84fa6d3b6882a283ddd1");
  tmp[1][0][1] = h256("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2");
  tmp[1][0][0] = h256("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed");
  tmp[1][1][1] = h256("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b");
  tmp[1][1][0] = h256("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa");
  alt_bn128_pairing(tmp[1], tmp[0][0], out[1]);
  fq12_mul(out[0], out[1], out[2]);
  if (!arreq(out[2], FQ12_ONE, 12)) {
    return 1;
  }

  return 0;
}

int test_alt_bn128_add() {
  // Taking from
  // https://github.com/ethereum/go-ethereum/blob/master/core/vm/testdata/precompiles/bn256Add.json
  uint256 a[2];
  uint256 b[2];
  uint256 r[2];

  a[0] = h256("0x18b18acfb4c2c30276db5411368e7185b311dd124691610c5d3b74034e093dc9");
  a[1] = h256("0x063c909c4720840cb5134cb9f59fa749755796819658d32efc0d288198f37266");
  b[0] = h256("0x07c2b7f58a84bd6145f00c9c2bc0bb1a187f20ff2c92963a88019e7c6a014eed");
  b[1] = h256("0x06614e20c147e940f2d70da3f74c9a17df361706a4485c742bd6788478fa17d7");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x2243525c5efd4b9c3d3c45ac0ca3fe4dd85e830a4ce6b65fa1eeaee202839703")) {
    return 1;
  }
  if (r[1] != h256("0x301d1d33be6da8e509df21cc35964723180eed7532537db9ae5e7d48f195c915")) {
    return 1;
  }

  a[0] = h256("0x2243525c5efd4b9c3d3c45ac0ca3fe4dd85e830a4ce6b65fa1eeaee202839703");
  a[1] = h256("0x301d1d33be6da8e509df21cc35964723180eed7532537db9ae5e7d48f195c915");
  b[0] = h256("0x18b18acfb4c2c30276db5411368e7185b311dd124691610c5d3b74034e093dc9");
  b[1] = h256("0x063c909c4720840cb5134cb9f59fa749755796819658d32efc0d288198f37266");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x2bd3e6d0f3b142924f5ca7b49ce5b9d54c4703d7ae5648e61d02268b1a0a9fb7")) {
    return 1;
  }
  if (r[1] != h256("0x21611ce0a6af85915e2f1d70300909ce2e49dfad4a4619c8390cae66cefdb204")) {
    return 1;
  }

  a[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  a[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  b[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  b[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x0000000000000000000000000000000000000000000000000000000000000000")) {
    return 1;
  }
  if (r[1] != h256("0x0000000000000000000000000000000000000000000000000000000000000000")) {
    return 1;
  }

  a[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  a[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  b[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000001");
  b[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000002");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x0000000000000000000000000000000000000000000000000000000000000001")) {
    return 1;
  }
  if (r[1] != h256("0x0000000000000000000000000000000000000000000000000000000000000002")) {
    return 1;
  }

  a[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000001");
  a[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000002");
  b[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  b[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000000");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x0000000000000000000000000000000000000000000000000000000000000001")) {
    return 1;
  }
  if (r[1] != h256("0x0000000000000000000000000000000000000000000000000000000000000002")) {
    return 1;
  }

  a[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000001");
  a[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000002");
  b[0] = h256("0x0000000000000000000000000000000000000000000000000000000000000001");
  b[1] = h256("0x0000000000000000000000000000000000000000000000000000000000000002");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd3")) {
    return 1;
  }
  if (r[1] != h256("0x15ed738c0e0a7c92e7845f96b2ae9c0a68a6a449e3538fc7ff3ebf7a5a18a2c4")) {
    return 1;
  }

  a[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  a[1] = h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c");
  b[0] = h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869");
  b[1] = h256("0x073a5ffcc6fc7a28c30723d6e58ce577356982d65b833a5a5c15bf9024b43d98");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x15bf2bb17880144b5d1cd2b1f46eff9d617bffd1ca57c37fb5a49bd84e53cf66")) {
    return 1;
  }
  if (r[1] != h256("0x049c797f9ce0d17083deb32b5e36f2ea2a212ee036598dd7624c168993d1355f")) {
    return 1;
  }

  a[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  a[1] = h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c");
  b[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  b[1] = h256("0x2e83f8d734803fc370eba25ed1f6b8768bd6d83887b87165fc2434fe11a830cb");
  alt_bn128_add(a, b, r);
  if (r[0] != h256("0x0000000000000000000000000000000000000000000000000000000000000000")) {
    return 1;
  }
  if (r[1] != h256("0x0000000000000000000000000000000000000000000000000000000000000000")) {
    return 1;
  }

  return 0;
}

int test_alt_bn128_mul() {
  // Taking from
  // https://github.com/ethereum/go-ethereum/blob/master/core/vm/testdata/precompiles/bn256ScalarMul.json
  uint256 a[2];
  uint256 b;
  uint256 r[2];

  a[0] = h256("0x2bd3e6d0f3b142924f5ca7b49ce5b9d54c4703d7ae5648e61d02268b1a0a9fb7");
  a[1] = h256("0x21611ce0a6af85915e2f1d70300909ce2e49dfad4a4619c8390cae66cefdb204");
  b = h256("0x00000000000000000000000000000000000000000000000011138ce750fa15c2");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x070a8d6a982153cae4be29d434e8faef8a47b274a053f5a4ee2a6c9c13c31e5c")) {
    return 1;
  }
  if (r[1] != h256("0x031b8ce914eba3a9ffb989f9cdd5b0f01943074bf4f0f315690ec3cec6981afc")) {
    return 1;
  }

  a[0] = h256("0x070a8d6a982153cae4be29d434e8faef8a47b274a053f5a4ee2a6c9c13c31e5c");
  a[1] = h256("0x031b8ce914eba3a9ffb989f9cdd5b0f01943074bf4f0f315690ec3cec6981afc");
  b = h256("0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd46");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x025a6f4181d2b4ea8b724290ffb40156eb0adb514c688556eb79cdea0752c2bb")) {
    return 1;
  }
  if (r[1] != h256("0x2eff3f31dea215f1eb86023a133a996eb6300b44da664d64251d05381bb8a02e")) {
    return 1;
  }

  a[0] = h256("0x025a6f4181d2b4ea8b724290ffb40156eb0adb514c688556eb79cdea0752c2bb");
  a[1] = h256("0x2eff3f31dea215f1eb86023a133a996eb6300b44da664d64251d05381bb8a02e");
  b = h256("0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x14789d0d4a730b354403b5fac948113739e276c23e0258d8596ee72f9cd9d323")) {
    return 1;
  }
  if (r[1] != h256("0x0af18a63153e0ec25ff9f2951dd3fa90ed0197bfef6e2a1a62b5095b9d2b4a27")) {
    return 1;
  }

  a[0] = h256("0x1a87b0584ce92f4593d161480614f2989035225609f08058ccfa3d0f940febe3");
  a[1] = h256("0x1a2f3c951f6dadcc7ee9007dff81504b0fcd6d7cf59996efdc33d92bf7f9f8f6");
  b = h256("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x2cde5879ba6f13c0b5aa4ef627f159a3347df9722efce88a9afbb20b763b4c41")) {
    return 1;
  }
  if (r[1] != h256("0x1aa7e43076f6aee272755a7f9b84832e71559ba0d2e0b17d5f9f01755e5b0d11")) {
    return 1;
  }

  a[0] = h256("0x1a87b0584ce92f4593d161480614f2989035225609f08058ccfa3d0f940febe3");
  a[1] = h256("0x1a2f3c951f6dadcc7ee9007dff81504b0fcd6d7cf59996efdc33d92bf7f9f8f6");
  b = h256("0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000000");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x1a87b0584ce92f4593d161480614f2989035225609f08058ccfa3d0f940febe3")) {
    return 1;
  }
  if (r[1] != h256("0x163511ddc1c3f25d396745388200081287b3fd1472d8339d5fecb2eae0830451")) {
    return 1;
  }

  a[0] = h256("0x1a87b0584ce92f4593d161480614f2989035225609f08058ccfa3d0f940febe3");
  a[1] = h256("0x1a2f3c951f6dadcc7ee9007dff81504b0fcd6d7cf59996efdc33d92bf7f9f8f6");
  b = h256("0x0000000000000000000000000000000100000000000000000000000000000000");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x1051acb0700ec6d42a88215852d582efbaef31529b6fcbc3277b5c1b300f5cf0")) {
    return 1;
  }
  if (r[1] != h256("0x135b2394bb45ab04b8bd7611bd2dfe1de6a4e6e2ccea1ea1955f577cd66af85b")) {
    return 1;
  }

  a[0] = h256("0x1a87b0584ce92f4593d161480614f2989035225609f08058ccfa3d0f940febe3");
  a[1] = h256("0x1a2f3c951f6dadcc7ee9007dff81504b0fcd6d7cf59996efdc33d92bf7f9f8f6");
  b = h256("0x0000000000000000000000000000000000000000000000000000000000000009");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x1dbad7d39dbc56379f78fac1bca147dc8e66de1b9d183c7b167351bfe0aeab74")) {
    return 1;
  }
  if (r[1] != h256("0x2cd757d51289cd8dbd0acf9e673ad67d0f0a89f912af47ed1be53664f5692575")) {
    return 1;
  }

  a[0] = h256("0x1a87b0584ce92f4593d161480614f2989035225609f08058ccfa3d0f940febe3");
  a[1] = h256("0x1a2f3c951f6dadcc7ee9007dff81504b0fcd6d7cf59996efdc33d92bf7f9f8f6");
  b = h256("0x0000000000000000000000000000000000000000000000000000000000000001");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x1a87b0584ce92f4593d161480614f2989035225609f08058ccfa3d0f940febe3")) {
    return 1;
  }
  if (r[1] != h256("0x1a2f3c951f6dadcc7ee9007dff81504b0fcd6d7cf59996efdc33d92bf7f9f8f6")) {
    return 1;
  }

  a[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  a[1] = h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c");
  b = h256("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x29e587aadd7c06722aabba753017c093f70ba7eb1f1c0104ec0564e7e3e21f60")) {
    return 1;
  }
  if (r[1] != h256("0x22b1143f6a41008e7755c71c3d00b6b915d386de21783ef590486d8afa8453b1")) {
    return 1;
  }

  a[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  a[1] = h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c");
  b = h256("0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000000");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9")) {
    return 1;
  }
  if (r[1] != h256("0x2e83f8d734803fc370eba25ed1f6b8768bd6d83887b87165fc2434fe11a830cb")) {
    return 1;
  }

  a[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  a[1] = h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c");
  b = h256("0x0000000000000000000000000000000100000000000000000000000000000000");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x221a3577763877920d0d14a91cd59b9479f83b87a653bb41f82a3f6f120cea7c")) {
    return 1;
  }
  if (r[1] != h256("0x2752c7f64cdd7f0e494bff7b60419f242210f2026ed2ec70f89f78a4c56a1f15")) {
    return 1;
  }

  a[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  a[1] = h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c");
  b = h256("0x0000000000000000000000000000000000000000000000000000000000000009");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x228e687a379ba154554040f8821f4e41ee2be287c201aa9c3bc02c9dd12f1e69")) {
    return 1;
  }
  if (r[1] != h256("0x1e0fd6ee672d04cfd924ed8fdc7ba5f2d06c53c1edc30f65f2af5a5b97f0a76a")) {
    return 1;
  }

  a[0] = h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9");
  a[1] = h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c");
  b = h256("0x0000000000000000000000000000000000000000000000000000000000000001");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x17c139df0efee0f766bc0204762b774362e4ded88953a39ce849a8a7fa163fa9")) {
    return 1;
  }
  if (r[1] != h256("0x01e0559bacb160664764a357af8a9fe70baa9258e0b959273ffc5718c6d4cc7c")) {
    return 1;
  }

  a[0] = h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869");
  a[1] = h256("0x073a5ffcc6fc7a28c30723d6e58ce577356982d65b833a5a5c15bf9024b43d98");
  b = h256("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x00a1a234d08efaa2616607e31eca1980128b00b415c845ff25bba3afcb81dc00")) {
    return 1;
  }
  if (r[1] != h256("0x242077290ed33906aeb8e42fd98c41bcb9057ba03421af3f2d08cfc441186024")) {
    return 1;
  }

  a[0] = h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869");
  a[1] = h256("0x073a5ffcc6fc7a28c30723d6e58ce577356982d65b833a5a5c15bf9024b43d98");
  b = h256("0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000000");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869")) {
    return 1;
  }
  if (r[1] != h256("0x2929ee761a352600f54921df9bf472e66217e7bb0cee9032e00acc86b3c8bfaf")) {
    return 1;
  }

  a[0] = h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869");
  a[1] = h256("0x073a5ffcc6fc7a28c30723d6e58ce577356982d65b833a5a5c15bf9024b43d98");
  b = h256("0x0000000000000000000000000000000100000000000000000000000000000000");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x1071b63011e8c222c5a771dfa03c2e11aac9666dd097f2c620852c3951a4376a")) {
    return 1;
  }
  if (r[1] != h256("0x2f46fe2f73e1cf310a168d56baa5575a8319389d7bfa6b29ee2d908305791434")) {
    return 1;
  }

  a[0] = h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869");
  a[1] = h256("0x073a5ffcc6fc7a28c30723d6e58ce577356982d65b833a5a5c15bf9024b43d98");
  b = h256("0x0000000000000000000000000000000000000000000000000000000000000009");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x19f75b9dd68c080a688774a6213f131e3052bd353a304a189d7a2ee367e3c258")) {
    return 1;
  }
  if (r[1] != h256("0x2612f545fb9fc89fde80fd81c68fc7dcb27fea5fc124eeda69433cf5c46d2d7f")) {
    return 1;
  }

  a[0] = h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869");
  a[1] = h256("0x073a5ffcc6fc7a28c30723d6e58ce577356982d65b833a5a5c15bf9024b43d98");
  b = h256("0x0000000000000000000000000000000000000000000000000000000000000001");
  alt_bn128_mul(a, b, r);
  if (r[0] != h256("0x039730ea8dff1254c0fee9c0ea777d29a9c710b7e616683f194f18c43b43b869")) {
    return 1;
  }
  if (r[1] != h256("0x073a5ffcc6fc7a28c30723d6e58ce577356982d65b833a5a5c15bf9024b43d98")) {
    return 1;
  }

  return 0;
}

int main() {
  if (test_invmod())
    return 1;
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
  if (test_alt_bn128_add())
    return 1;
  if (test_alt_bn128_mul())
    return 1;

  return 0;
}
