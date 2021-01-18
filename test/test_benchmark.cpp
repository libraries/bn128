#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

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

int main() {
  if (test_fq2())
    return 1;
  return 0;
}
