#include <bn128.hpp>
#include <intx/intx.hpp>

using namespace bn128;

int main() {
  int r = static_cast<int>(_powmod(100, 500, FIELD_MODULUS));
  assert(r == 1522359559);
  return 0;
}
