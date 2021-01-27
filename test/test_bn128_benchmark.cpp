#include <bn128_mont.hpp>
#include <intx/intx.hpp>

using namespace bn128;

int main() {
  // Taking from
  // https://github.com/xxuejie/benchmarking-wasm-ewasm-evm/blob/checkpoint/evmrace/ckbvm/bn256g2_test.cpp
  G2Affine a = G2Affine{
    x : FQ2(FQ(mont_encode(h256("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed"))),
            FQ(mont_encode(h256("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2")))),
    y : FQ2(FQ(mont_encode(h256("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa"))),
            FQ(mont_encode(h256("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b")))),
  };
  G2 b = a.into();
  G2 c = b.mul(0x2dddefa19);
  if (c.x.c0 != FQ(h256("0x1dafb20da22cc10ef4736831e8a0f1419a9753c43d8400f2949cb68a2d381719"))) {
    return 1;
  }
  if (c.x.c1 != FQ(h256("0x047d62b471a91cfa91fb08d4344792a0355cc0471784adb20e9a99bf11286133"))) {
    return 1;
  }
  if (c.y.c0 != FQ(h256("0x1b5694e95d08d04325f5ebb31dd1889a402e9e76de68700071cd0506326f8ec3"))) {
    return 1;
  }
  if (c.y.c1 != FQ(h256("0x19c7a2f5e6eedf85ec5adff3e0c2028a52b78e15a2e081513a5e5c6309cf7cc2"))) {
    return 1;
  }
  if (c.z.c0 != FQ(h256("0x1415a498e2d14c0157076be278709572d52e168d72ac62082e53aac308d50cbe"))) {
    return 1;
  }
  if (c.z.c1 != FQ(h256("0x1363bae5bce3ce99dbcc3a7419295dcb54434a925f3933a7342aa8d9077925a3"))) {
    return 1;
  }
  return 0;
}
