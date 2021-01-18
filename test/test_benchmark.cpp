#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

int test_fq12() {
  uint256 x[12] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  uint256 f[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  uint256 fpx[12] = {2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  uint256 one[12] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  uint256 tmp[8][12] = {};

  // Assert x + f == fpx
  // fq12_add(x, f, tmp[0]);
  // if (!eq12(tmp[0], fpx)) {
  //   return 1;
  // }

  // // Assert f / f == one
  // fq12_div(f, f, tmp[0]);
  // if (!eq12(tmp[0], one)) {
  //   return 1;
  // }

  // // Assert one / f + x / f == (one + x) / f
  // fq12_div(one, f, tmp[0]);
  // fq12_div(x, f, tmp[1]);
  // fq12_add(tmp[0], tmp[1], tmp[2]);
  // fq12_add(one, x, tmp[0]);
  // fq12_div(tmp[0], f, tmp[1]);
  // if (!eq12(tmp[1], tmp[2])) {
  //   return 1;
  // }

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
  // fq12_pow(x, fq_sub(fq_pow(FIELD_MODULUS, 12), 1), tmp[0]);
  // if (!eq12(tmp[0], one)) {
  //   return 1;
  // }

  return 0;
}

int main() {
  if (test_fq12())
    return 1;
  return 0;
}
