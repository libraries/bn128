#include <bn128_mont.hpp>
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
  if (mont_encode(1) != FQ_ONE) {
    return 1;
  }
  // G2_COEFF_B0 = mont_encode(0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5)
  if (mont_encode(h256("0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5")) != G2_COEFF_B[0]) {
    return 1;
  }
  // G2_COEFF_B1 = mont_encode(0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2)
  if (mont_encode(h256("0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2")) != G2_COEFF_B[1]) {
    return 1;
  }

  return 0;
}

int test_mont_encode_decode() {
  // Taking from https://github.com/paritytech/bn
  uint256 a = mont_encode(1);
  if (a != h256("0xe0a77c19a07df2f666ea36f7879462c0a78eb28f5c70b3dd35d438dc58f0d9d")) {
    return 1;
  }
  if (mont_decode(a) != 1) {
    return 1;
  }

  uint256 b = mont_encode(42);
  if (b != h256("0x0903f860b6f71bd22a638bbbb1d55ed69dc595e76d525985dbc68430439c5c6e")) {
    return 1;
  }
  if (mont_decode(b) != 42) {
    return 1;
  }

  uint256 c = FIELD_MODULUS / 2 - 17;
  uint256 d = mont_encode(c);
  if (d != h256("0x14707fbbcf072f27f529534d0bfd1a000a03b6d2f156954ed7d2e44ca56802cb")) {
    return 1;
  }
  if (mont_decode(d) != c) {
    return 1;
  }

  return 0;
}

int test_fq_inv() {
  uint256 a = mont_encode(10);
  uint256 b = fq_inv(a);
  // a = 0x2ba010aa41eb77868fb1d6edb1ba0cfd39b65a76c8e2db4fc9638b5c069c8d94
  // b = 0x0febf04f6c8facc474bcbedb7fb2ee53ae7fb77c84b60ab08d93172ea1809a25
  if (b != h256("0x0febf04f6c8facc474bcbedb7fb2ee53ae7fb77c84b60ab08d93172ea1809a25")) {
    return 1;
  }
  return 0;
}

int test_fq_neg() {
  uint256 a = mont_encode(42);
  uint256 b = fq_neg(a);
  // a = 0x0903f860b6f71bd22a638bbbb1d55ed69dc595e76d525985dbc68430439c5c6e
  // b = 0x276056122a3a84578decb9facfabf986f9bbd4a9fb1f7107605a07e694e0a0d9
  if (b != h256("0x276056122a3a84578decb9facfabf986f9bbd4a9fb1f7107605a07e694e0a0d9")) {
    return 1;
  }
  return 0;
}

int test_fq2_mul() {
  uint256 a[2] = {
      h256("0x0010b52d9fe70d08c967a97deeb9eb186da14c608196f376d63ca9589ca5990e"),
      h256("0x2f682d1f7dda8678b0d017978b3067b74807a5d49d2a41739659c6600a8bf018"),
  };
  uint256 b[2] = {
      h256("0x19015293267c8307ffb557fd4ad6052cc22e04f121f21da65bbe2a733c22c53d"),
      h256("0x0e2bf3144f8ca0808b1dfce33c2240c641ff2f7b8f2ca8c185201c5edcc67040"),
  };
  uint256 c[2];
  fq2_mul(a, b, c);
  if (c[0] != h256("0x1a458f1555acd5430609a64acd087c155541125d671ced7f65dcb5a4e48e7d6d")) {
    return 1;
  }
  if (c[1] != h256("0x8a9a40e08f6e1ee5c2997247f30e5d49235ef40a999734c3ae55e9726f7f1c8")) {
    return 1;
  }

  return 0;
}

int test_fq2_inv() {
  uint256 a[2] = {
      h256("0x0010b52d9fe70d08c967a97deeb9eb186da14c608196f376d63ca9589ca5990e"),
      h256("0x2f682d1f7dda8678b0d017978b3067b74807a5d49d2a41739659c6600a8bf018"),
  };
  uint256 b[2];
  fq2_inv(a, b);
  if (b[0] != h256("0x2e99d10b04272627f24f9f0c41e928004b3bdc31880830c2ff2bbe19546ec5d1")) {
    return 1;
  }
  if (b[1] != h256("0x2b5c2767a083b12aae7fe0d7b0682b5813f1a5566bc3e345b985aa42dd595820")) {
    return 1;
  }
  return 0;
}

int test_fq2_squr() {
  uint256 a[2] = {
      h256("0x0020b52d9fe70d08c967a97deeb9eb186da14c608196f376d63ca9589ca5970f"),
      h256("0x2b782d1f7dda8678b0d017978b3067b74807a5d49d2a41739659c6600a8bf015"),
  };
  uint256 b[2];
  fq2_squr(a, b);
  uint256 c[2];
  fq2_mul(a, a, c);
  if (!arrequ(b, c, 2)) {
    return 1;
  }

  uint256 d[2] = {10, 20};
  uint256 e[2];
  uint256 f[2];
  fq2_squr(d, e);
  fq2_mul(d, d, f);
  if (!arrequ(e, f, 2)) {
    return 1;
  }

  return 0;
}

int test_g2_jacobian_affine_conv() {
  uint256 a00 = mont_encode(h256("0x1ecfd2dff2aad18798b64bdb0c2b50c9d73e6c05619e04cbf5b448fd98726880"));
  uint256 a01 = mont_encode(h256("0x0e16c8d96362720af0916592be1b839a26f5e6b710f3ede0d8840d9a70eaf97f"));
  uint256 a10 = mont_encode(h256("0x2aa778acda9e7d4925c60ad84c12fb3b4f2b9539d5699934b0e6fdd10cc2c0e1"));
  uint256 a11 = mont_encode(h256("0x1e8f2c1f441fed039bb46d6bfb91236cf7ba240c75080cedbe40e049c46b26be"));
  uint256 a[2][2] = {{a00, a01}, {a10, a11}};
  uint256 b[3][2];
  uint256 c[2][2];
  g2_from_affine(a, b);
  if (b[0][0] != h256("0x180b9347a1a7d0a8fdaf0bced1c1762460053210c288b08a98bf18dbac6f2e8f")) {
    return 1;
  }
  if (b[0][1] != h256("0x20695a8c8acfdc037e47d786b0615dd17f3f17b3f63bb0ad9fd68d54d7bbf4ed")) {
    return 1;
  }
  if (b[1][0] != h256("0x29dad7762754564b72ac95ba67b28db571932de27a4d188d7704dd18f28cfcbe")) {
    return 1;
  }
  if (b[1][1] != h256("0x2c5d04ea5da8d50c233667e8e1086c089a67f531e25789bfbaeb22635e0a5f73")) {
    return 1;
  }
  if (b[2][0] != FQ2_ONE[0]) {
    return 1;
  }
  if (b[2][1] != FQ2_ONE[1]) {
    return 1;
  }
  g2_from_jacobian(b, c);
  if (!arrequ(a[0], c[0], 2) || !arrequ(a[1], c[1], 2)) {
    return 1;
  }

  return 0;
}

int test_g2_double() {
  uint256 a[2][2];
  a[0][0] = mont_encode(h256("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed"));
  a[0][1] = mont_encode(h256("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2"));
  a[1][0] = mont_encode(h256("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa"));
  a[1][1] = mont_encode(h256("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b"));
  uint256 b[3][2];
  g2_from_affine(a, b);
  uint256 c[3][2];
  g2_double(b, c);
  if (b[0][0] != h256("0x19573841af96503bfbb8264797811adfdceb1935497b01728e83b5d102bc2026")) {
    return 1;
  }
  if (b[0][1] != h256("0x14fef0833aea7b6b09e950fc52a02f866043dd5a5802d8c4afb4737da84c6140")) {
    return 1;
  }
  if (b[1][0] != h256("0x28fd7eebae9e4206ff9e1a62231b7dfefe7fd297f59e9b78619dfa9d886be9f6")) {
    return 1;
  }
  if (b[1][1] != h256("0xda4a0e693fd648255f935be33351076dc57f922327d3cbb64095b56c71856ee")) {
    return 1;
  }
  if (b[2][0] != h256("0xe0a77c19a07df2f666ea36f7879462c0a78eb28f5c70b3dd35d438dc58f0d9d")) {
    return 1;
  }
  if (b[2][1] != 0) {
    return 1;
  }
  if (c[0][0] != h256("0x190178401674fd1a6ce2b0feabb43601f71985df6e1d9ae392c552517113f04c")) {
    return 1;
  }
  if (c[0][1] != h256("0x22c5dc5d9fe664ed64074fdc1d7eba9d5be249895ca71ac4192c5eacbcc23fc5")) {
    return 1;
  }
  if (c[1][0] != h256("0x16b5ca02a81579a0b7c93f9ab34a255d29737922388c9abf69faf1822d82e9d2")) {
    return 1;
  }
  if (c[1][1] != h256("0x36e1261fdadb4eead4952efe416e4ba9a5d9892c8fce6864a14ebd2219f1b7")) {
    return 1;
  }
  if (c[2][0] != h256("0x2196af647c0ae3e446ebef0dc4b5a3a0657e3a9e82cb6c63871b6924385ad6a5")) {
    return 1;
  }
  if (c[2][1] != h256("0x1b4941cd27fac904abf26b7c666a20edb8aff24464fa7976c812b6ad8e30addc")) {
    return 1;
  }
  return 0;
}

int test_g2_add() {
  uint256 a[2][2];
  a[0][0] = mont_encode(h256("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed"));
  a[0][1] = mont_encode(h256("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2"));
  a[1][0] = mont_encode(h256("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa"));
  a[1][1] = mont_encode(h256("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b"));
  uint256 b[3][2];
  g2_from_affine(a, b);
  uint256 c[3][2];
  g2_add(b, G2_ONE, c);
  if (c[0][0] != h256("0x190178401674fd1a6ce2b0feabb43601f71985df6e1d9ae392c552517113f04c")) {
    return 1;
  }
  if (c[0][1] != h256("0x22c5dc5d9fe664ed64074fdc1d7eba9d5be249895ca71ac4192c5eacbcc23fc5")) {
    return 1;
  }
  if (c[1][0] != h256("0x16b5ca02a81579a0b7c93f9ab34a255d29737922388c9abf69faf1822d82e9d2")) {
    return 1;
  }
  if (c[1][1] != h256("0x36e1261fdadb4eead4952efe416e4ba9a5d9892c8fce6864a14ebd2219f1b7")) {
    return 1;
  }
  if (c[2][0] != h256("0x2196af647c0ae3e446ebef0dc4b5a3a0657e3a9e82cb6c63871b6924385ad6a5")) {
    return 1;
  }
  if (c[2][1] != h256("0x1b4941cd27fac904abf26b7c666a20edb8aff24464fa7976c812b6ad8e30addc")) {
    return 1;
  }

  uint256 d[3][2];
  g2_add(G2_ZERO, b, d);
  if (b[0][0] != d[0][0]) {
    return 1;
  }
  if (b[0][1] != d[0][1]) {
    return 1;
  }
  if (b[1][0] != d[1][0]) {
    return 1;
  }
  if (b[1][1] != d[1][1]) {
    return 1;
  }
  if (b[2][0] != d[2][0]) {
    return 1;
  }
  if (b[2][1] != d[2][1]) {
    return 1;
  }

  uint256 e[3][2];
  uint256 f[3][2];
  g2_add(b, b, e);
  g2_add(e, b, f);
  if (f[0][0] != h256("0x178ccafd45302ad9a9cc2c2619f354a546e17d9ef882176268874584435a174d")) {
    return 1;
  }
  if (f[0][1] != h256("0x22bd6dd71372cf61c0a45ba029d8286d2078e56e4dfb707c3be17d3082bd6313")) {
    return 1;
  }
  if (f[1][0] != h256("0x07f6555d95103d4a9ad6e833c966af819141b80fb9f0c736ba3155bb745fd473")) {
    return 1;
  }
  if (f[1][1] != h256("0x25472a0de9dc5cdbc176677542c925c01484160af28eb09b802aa2eb3a34c188")) {
    return 1;
  }
  if (f[2][0] != h256("0x0bc2acc98cf51700e0ece64c8e8888cd4fbee1ca0edf4a361eb498df3106e29f")) {
    return 1;
  }
  if (f[2][1] != h256("0x1ba16ff98d805609ceff8f1a56f25372514642a88c62b59b52558fb3a5b0397c")) {
    return 1;
  }

  return 0;
}

int test_g2_mul() {
  // Taking from
  // https://github.com/xxuejie/benchmarking-wasm-ewasm-evm/blob/checkpoint/evmrace/ckbvm/bn256g2_test.cpp
  uint256 a[2][2];
  a[0][0] = mont_encode(h256("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed"));
  a[0][1] = mont_encode(h256("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2"));
  a[1][0] = mont_encode(h256("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa"));
  a[1][1] = mont_encode(h256("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b"));
  uint256 b[3][2];
  g2_from_affine(a, b);
  uint256 c[3][2];
  g2_mul(b, 0x2dddefa19, c);
  if (c[0][0] != h256("0x1dafb20da22cc10ef4736831e8a0f1419a9753c43d8400f2949cb68a2d381719")) {
    return 1;
  }
  if (c[0][1] != h256("0x047d62b471a91cfa91fb08d4344792a0355cc0471784adb20e9a99bf11286133")) {
    return 1;
  }
  if (c[1][0] != h256("0x1b5694e95d08d04325f5ebb31dd1889a402e9e76de68700071cd0506326f8ec3")) {
    return 1;
  }
  if (c[1][1] != h256("0x19c7a2f5e6eedf85ec5adff3e0c2028a52b78e15a2e081513a5e5c6309cf7cc2")) {
    return 1;
  }
  if (c[2][0] != h256("0x1415a498e2d14c0157076be278709572d52e168d72ac62082e53aac308d50cbe")) {
    return 1;
  }
  if (c[2][1] != h256("0x1363bae5bce3ce99dbcc3a7419295dcb54434a925f3933a7342aa8d9077925a3")) {
    return 1;
  }
  uint256 d[2][2];
  g2_from_jacobian(c, d);
  if (mont_decode(d[0][0]) != h256("0x23997083c2c4409869ee3546806a544c8c16bc46cc88598c4e1c853eb81d45b0")) {
    return 1;
  }
  if (mont_decode(d[0][1]) != h256("0x1142585a23028cbe57783f890d1a2f6837049fce43c9b3b5e8e14c40a43c617a")) {
    return 1;
  }
  if (mont_decode(d[1][0]) != h256("0x215a23c8a96e1ca11d52cf6e2d6ada4ed01ee7e09b06dbc7f3315e7e6e73b919")) {
    return 1;
  }
  if (mont_decode(d[1][1]) != h256("0x0edac9f3a977530e28d4a385e614bcb7a8f9c3c3cb65707c1b90b5ea86174512")) {
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
  if (test_fq_inv())
    return 1;
  if (test_fq_neg())
    return 1;
  if (test_fq2_mul())
    return 1;
  if (test_fq2_inv())
    return 1;
  if (test_fq2_squr())
    return 1;
  if (test_g2_jacobian_affine_conv())
    return 1;
  if (test_g2_double())
    return 1;
  if (test_g2_add())
    return 1;
  if (test_g2_mul())
    return 1;
  return 0;
}
