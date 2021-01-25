#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

int main() {
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
  return 0;
}
