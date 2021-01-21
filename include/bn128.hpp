#ifndef BN128_H_
#define BN128_H_

#include <intx/intx.hpp>

namespace bn128 {

// Maybe there is a better way to implement this macro, but this is enough for
// now.
#ifdef __riscv
#undef assert
#define assert(x)                                                              \
  if (!x)                                                                      \
  exit(2)
#endif

using uint256 = intx::uint256;

inline bool eq2(const uint256 x[2], const uint256 y[2]) {
  return x[0] == y[0] && x[1] == y[1];
}

inline bool eq12(const uint256 x[12], const uint256 y[12]) {
  return x[0] == y[0] && x[1] == y[1] && x[2] == y[2] && x[3] == y[3] &&
         x[4] == y[4] && x[5] == y[5] && x[6] == y[6] && x[7] == y[7] &&
         x[8] == y[8] && x[9] == y[9] && x[10] == y[10] && x[11] == y[11];
}

inline void cp2(const uint256 x[2], uint256 r[2]) {
  r[0] = x[0];
  r[1] = x[1];
}

inline void cp3(const uint256 x[3], uint256 r[3]) {
  r[0] = x[0];
  r[1] = x[1];
  r[2] = x[2];
}

inline void cp12(const uint256 x[12], uint256 r[12]) {
  r[0] = x[0];
  r[1] = x[1];
  r[2] = x[2];
  r[3] = x[3];
  r[4] = x[4];
  r[5] = x[5];
  r[6] = x[6];
  r[7] = x[7];
  r[8] = x[8];
  r[9] = x[9];
  r[10] = x[10];
  r[11] = x[11];
}

inline void cp13(const uint256 x[13], uint256 r[13]) {
  cp12(x, r);
  r[12] = x[12];
}

// The prime modulus of the field.
constexpr uint256 FIELD_MODULUS = intx::from_string<uint256>(
    "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47");

inline uint256 _addmod(const uint256 &x, const uint256 &y, const uint256 &n) {
  return intx::addmod(x, y, n);
}

inline uint256 _submod(const uint256 &x, const uint256 &y, const uint256 &n) {
  return _addmod(x, n - y, n);
}

inline uint256 _negmod(const uint256 &x, const uint256 &n) {
  return n - x;
}

inline uint256 _mulmod(const uint256 &x, const uint256 &y, const uint256 &n) {
  return intx::mulmod(x, y, n);
}

// Extended euclidean algorithm to find modular inverses for integers.
inline uint256 _invmod(const uint256 &x, const uint256 &n) {
  uint256 t = 0;
  uint256 newt = 1;
  uint256 r = n;
  uint256 newr = x;
  while (newr) {
    uint256 q = r / newr;
    uint256 oldt = t;
    t = newt;
    newt = _submod(oldt, _mulmod(q, newt, n), n);
    uint256 oldr = r;
    r = newr;
    newr = oldr - q * newr;
  }
  return t;
}

inline uint256 _divmod(const uint256 &x, const uint256 &y, const uint256 &n) {
  return _mulmod(x, _invmod(y, FIELD_MODULUS), FIELD_MODULUS);
}

inline uint256 _powmod(const uint256 &x, const uint256 &y, const uint256 &n) {
  if (y == 0) {
    return 1;
  } else if (y == 1) {
    return x;
  } else if (y & 1) {
    return _mulmod(_powmod(_mulmod(x, x, n), y >> 1, n), x, n);
  } else {
    return _powmod(_mulmod(x, x, n), y >> 1, n);
  }
}

inline uint256 fq_add(const uint256 &x, const uint256 &y) {
  return _addmod(x, y, FIELD_MODULUS);
}

inline uint256 fq_sub(const uint256 &x, const uint256 &y) {
  return _submod(x, y, FIELD_MODULUS);
}

inline uint256 fq_neg(const uint256 &x) { return _negmod(x, FIELD_MODULUS); }

inline uint256 fq_mul(const uint256 &x, const uint256 &y) {
  return _mulmod(x, y, FIELD_MODULUS);
}

inline uint256 fq_inv(const uint256 &x) { return _invmod(x, FIELD_MODULUS); }

inline uint256 fq_div(const uint256 &x, const uint256 &y) {
  return _divmod(x, y, FIELD_MODULUS);
}

inline uint256 fq_pow(const uint256 &x, const uint256 &y) {
  return _powmod(x, y, FIELD_MODULUS);
}

// The quadratic extension field.
constexpr uint256 FQ2_ONE[2] = {1, 0};
constexpr uint256 FQ2_ZERO[2] = {0, 0};

void fq2_add(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 a = fq_add(x[0], y[0]);
  uint256 b = fq_add(x[1], y[1]);
  r[0] = a;
  r[1] = b;
}

void fq2_sub(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 a = fq_sub(x[0], y[0]);
  uint256 b = fq_sub(x[1], y[1]);
  r[0] = a;
  r[1] = b;
}

void fq2_neg(const uint256 x[2], uint256 r[2]) {
  uint256 a = fq_neg(x[0]);
  uint256 b = fq_neg(x[1]);
  r[0] = a;
  r[1] = b;
}

void fq2_mul(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 a = fq_sub(fq_mul(x[0], y[0]), fq_mul(x[1], y[1]));
  uint256 b = fq_add(fq_mul(x[0], y[1]), fq_mul(x[1], y[0]));
  r[0] = a;
  r[1] = b;
}

void fq2_muc(const uint256 x[2], const uint256 &c, uint256 r[2]) {
  uint256 a = fq_mul(x[0], c);
  uint256 b = fq_mul(x[1], c);
  r[0] = a;
  r[1] = b;
}

void fq2_inv(const uint256 x[2], uint256 r[2]) {
  uint256 i = fq_inv(fq_add(fq_mul(x[0], x[0]), fq_mul(x[1], x[1])));
  uint256 a = fq_mul(x[0], i);
  uint256 b = fq_sub(FIELD_MODULUS, fq_mul(x[1], i));
  r[0] = a;
  r[1] = b;
}

void fq2_div(const uint256 x[2], const uint256 y[2], uint256 r[2]) {
  uint256 t[2];
  fq2_inv(y, t);
  fq2_mul(x, t, r);
}

void fq2_dic(const uint256 x[2], const uint256 &c, uint256 r[2]) {
  uint256 a = fq_div(x[0], c);
  uint256 b = fq_div(x[1], c);
  r[0] = a;
  r[1] = b;
}

void fq2_pow(const uint256 x[2], const uint256 &y, uint256 r[2]) {
  if (y == 0) {
    cp2(FQ2_ONE, r);
  } else if (y == 1) {
    cp2(x, r);
  } else if (y & 1) {
    uint256 t[2];
    fq2_mul(x, x, t);
    fq2_pow(t, y >> 1, t);
    fq2_mul(x, t, r);
  } else {
    fq2_mul(x, x, r);
    fq2_pow(r, y >> 1, r);
  }
}

// The 12th-degree extension field.
constexpr uint256 FQ12_ONE[12] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
constexpr uint256 FQ12_ZERO[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// The modulus of the polynomial in this representation of FQ12.
constexpr uint256 FQ12_MODULUS_COEFFS[12] = {
    82, 0, 0, 0, 0, 0, FIELD_MODULUS - 18, 0, 0, 0, 0, 0};

void fq12_add(const uint256 x[12], const uint256 y[12], uint256 r[12]) {
  for (int i = 0; i < 12; i++) {
    r[i] = fq_add(x[i], y[i]);
  }
}

void fq12_sub(const uint256 x[12], const uint256 y[12], uint256 r[12]) {
  for (int i = 0; i < 12; i++) {
    r[i] = fq_sub(x[i], y[i]);
  }
}

void fq12_neg(const uint256 x[12], uint256 r[12]) {
  for (int i = 0; i < 12; i++) {
    r[i] = fq_neg(x[i]);
  }
}

void fq12_mul(const uint256 x[12], const uint256 y[12], uint256 r[12]) {
  uint256 b[23] = {0};
  uint256 eli;
  uint256 elj;
  uint256 top;
  for (int i = 0; i < 12; i++) {
    eli = x[i];
    for (int j = 0; j < 12; j++) {
      elj = y[j];
      b[i + j] = fq_add(b[i + j], fq_mul(eli, elj));
    }
  }
  for (int exp = 10; exp > -1; exp--) {
    top = b[exp + 12];
    b[exp] = fq_sub(b[exp], fq_mul(top, FQ12_MODULUS_COEFFS[0]));
    b[exp + 6] = fq_sub(b[exp + 6], fq_mul(top, FQ12_MODULUS_COEFFS[6]));
  }
  cp12(b, r);
}

void fq12_muc(const uint256 x[12], const uint256 &c, uint256 r[12]) {
  for (int i = 0; i < 12; i++) {
    r[i] = fq_mul(x[i], c);
  }
}

// Utility methods for polynomial math.
int _deg(const uint256 p[13]) {
  int d = 12;
  while (p[d] == 0 && d != 0) {
    d--;
  }
  return d;
}

void _poly_rounded_div(const uint256 a[13], const uint256 b[13],
                       uint256 r[13]) {
  int dega = _deg(a);
  int degb = _deg(b);
  uint256 t[13] = {};
  cp13(a, t);
  uint256 o[13] = {0};
  for (int i = dega - degb; i > -1; i--) {
    o[i] = fq_add(o[i], fq_div(t[degb + i], b[degb]));
    for (int c = 0; c < degb + 1; c++) {
      t[c + i] = fq_sub(t[c + i], o[c]);
    }
  }
  int n = _deg(o) + 1;
  for (int i = 0; i < n; i++) {
    r[i] = o[i];
  }
}

void fq12_inv(const uint256 x[12], uint256 r[12]) {
  uint256 lm[13] = {1, 0};
  uint256 hm[13] = {0};
  uint256 low[13] = {};
  cp12(x, low);
  low[12] = 0;
  uint256 high[13] = {};
  cp12(FQ12_MODULUS_COEFFS, high);
  high[12] = 1;

  uint256 temp[13] = {};
  uint256 nm[13] = {};
  uint256 news[13] = {};
  while (_deg(low)) {
    for (int i = 0; i < 13; i++) {
      temp[i] = 0;
      nm[i] = hm[i];
      news[i] = high[i];
    }
    _poly_rounded_div(high, low, temp);
    for (int i = 0; i < 13; i++) {
      for (int j = 0; j < 13 - i; j++) {
        nm[i + j] = fq_sub(nm[i + j], fq_mul(lm[i], temp[j]));
        news[i + j] = fq_sub(news[i + j], fq_mul(low[i], temp[j]));
      }
    }
    for (int i = 0; i < 13; i++) {
      hm[i] = lm[i];
      lm[i] = nm[i];
      high[i] = low[i];
      low[i] = news[i];
    }
  }
  for (int i = 0; i < 12; i++) {
    r[i] = fq_div(lm[i], low[0]);
  }
}

void fq12_div(const uint256 x[12], const uint256 y[12], uint256 r[12]) {
  uint256 t[12];
  fq12_inv(y, t);
  fq12_mul(x, t, r);
}

void fq12_dic(const uint256 x[12], const uint256 &c, uint256 r[12]) {
  for (int i = 0; i < 12; i++) {
    r[i] = fq_div(x[i], c);
  }
}

void fq12_pow(const uint256 x[12], const uint256 &y, uint256 r[12]) {
  if (y == 0) {
    cp12(FQ12_ONE, r);
  } else if (y == 1) {
    cp12(x, r);
  } else if (y & 1) {
    uint256 t[12];
    fq12_mul(x, x, t);
    fq12_pow(t, y >> 1, t);
    fq12_mul(x, t, r);
  } else {
    fq12_mul(x, x, r);
    fq12_pow(r, y >> 1, r);
  }
}

constexpr uint256 CURVE_ORDER = intx::from_string<uint256>(
    "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001");
// Curve is y**2 = x**3 + 3
constexpr uint256 B = 3;
// Twisted curve over FQ**2
constexpr uint256 B2[2] = {
    intx::from_string<uint256>(
        "0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5"),
    intx::from_string<uint256>(
        "0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2"),
};
// Extension curve over FQ**12; same b value as over FQ
constexpr uint256 B12[12] = {3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// Generator for curve over FQ
constexpr uint256 G1[3] = {1, 2, 1};
// Generator for twisted curve over FQ2
constexpr uint256 G2[3][2] = {
    {
        intx::from_string<uint256>("0x1800deef121f1e76426a00665e5c4479674322d4f"
                                   "75edadd46debd5cd992f6ed"),
        intx::from_string<uint256>("0x198e9393920d483a7260bfb731fb5d25f1aa49333"
                                   "5a9e71297e485b7aef312c2"),
    },
    {
        intx::from_string<uint256>("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690"
                                   "c43d37b4ce6cc0166fa7daa"),
        intx::from_string<uint256>("0x090689d0585ff075ec9e99ad690c3395bc4b31337"
                                   "0b38ef355acdadcd122975b"),
    }, {1, 0}};

// "Twist" a point in E(FQ2) into a point in E(FQ12)
constexpr uint256 W[12] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

constexpr uint256 G12[3][12] = {
    {0, 0,
     intx::from_string<uint256>(
         "0x23f336fd559fb538d6949f86240cb7f7ddcda4df1e9eaff81c78c659ed78407e"),
     0, 0, 0, 0, 0,
     intx::from_string<uint256>(
         "0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2"),
     0, 0, 0},
    {0, 0, 0,
     intx::from_string<uint256>(
         "0x2256233882903a1969b895d4df602107743001bce6d76207c214326bbdbd2605"),
     0, 0, 0, 0, 0,
     intx::from_string<uint256>(
         "0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b"),
     0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

namespace g1 {

// Check if a point is the point at infinity
inline bool is_inf(const uint256 pt[3]) { return pt[2] == 0; }

// Check that a point is on the curve defined by y**2 == x**3 + b
bool is_on_curve(const uint256 pt[3]) {
  if (is_inf(pt)) {
    return true;
  }
  uint256 x = pt[0], y = pt[1], z = pt[2];
  uint256 l = fq_sub(fq_mul(fq_mul(y, y), z), fq_mul(fq_mul(x, x), x));
  uint256 r = fq_mul(B, fq_mul(fq_mul(z, z), z));
  return l == r;
}

bool eq(const uint256 x[3], const uint256 y[3]) {
  uint256 x1 = x[0], y1 = x[1], z1 = x[2];
  uint256 x2 = y[0], y2 = y[1], z2 = y[2];
  return (fq_mul(x1, z2) == fq_mul(x2, z1)) && (fq_mul(y1, z2) == fq_mul(y2, z1));
}

void doubl2(const uint256 pt[3], uint256 r[3]) {
  uint256 x = pt[0], y = pt[1], z = pt[2];

  // W = 3 * x * x
  uint256 W = fq_mul(fq_mul(3, x), x);
  // S = y * z
  uint256 S = fq_mul(y, z);
  // B = x * y * S
  uint256 B = fq_mul(fq_mul(x, y), S);
  // H = W * W - 8 * B
  uint256 H = fq_sub(fq_mul(W, W), fq_mul(8, B));
  // S_squared = S * S
  uint256 S_squared = fq_mul(S, S);
  // newx = 2 * H * S
  uint256 newx = fq_mul(2, fq_mul(H, S));
  // newy = W * (4 * B - H) - 8 * y * y * S_squared
  uint256 newy_l = fq_mul(W, fq_sub(fq_mul(4, B), H));
  uint256 newy_r = fq_mul(fq_mul(fq_mul(8, y), y), S_squared);
  uint256 newy = fq_sub(newy_l, newy_r);
  // newz = 8 * S * S_squared
  uint256 newz = fq_mul(fq_mul(8, S), S_squared);
  r[0] = newx;
  r[1] = newy;
  r[2] = newz;
}

// Elliptic curve addition
void add(const uint256 p1[3], const uint256 p2[3], uint256 r[3]) {
  if (p1[2] == 0) {
    cp3(p2, r);
    return;
  }
  if (p2[2] == 0) {
    cp3(p1, r);
    return;
  }
  uint256 x1 = p1[0], y1 = p1[1], z1 = p1[2];
  uint256 x2 = p2[0], y2 = p2[1], z2 = p2[2];
  // U1 = y2 * z1
  uint256 U1 = fq_mul(y2, z1);
  // U2 = y1 * z2
  uint256 U2 = fq_mul(y1, z2);
  // V1 = x2 * z1
  uint256 V1 = fq_mul(x2, z1);
  // V2 = x1 * z2
  uint256 V2 = fq_mul(x1, z2);

  if (V1 == V2) {
    if (U1 == U2) {
      doubl2(p1, r);
    } else {
      r[0] = 1;
      r[1] = 1;
      r[2] = 0;
    }
    return;
  }

  // U = U1 - U2
  uint256 U = fq_sub(U1, U2);
  // V = V1 - V2
  uint256 V = fq_sub(V1, V2);
  // V_squared = V * V
  uint256 V_squared = fq_mul(V, V);
  // V_squared_times_V2 = V_squared * V2
  uint256 V_squared_times_V2 = fq_mul(V_squared, V2);
  // V_cubed = V * V_squared
  uint256 V_cubed = fq_mul(V, V_squared);
  // W = z1 * z2
  uint256 W = fq_mul(z1, z2);
  // A = U * U * W - V_cubed - 2 * V_squared_times_V2
  uint256 A = fq_sub(fq_sub(fq_mul(fq_mul(U, U), W), V_cubed), fq_mul(2, V_squared_times_V2));
  // newx = V * A
  uint256 newx = fq_mul(V, A);
  // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
  uint256 newy = fq_sub(fq_mul(U, fq_sub(V_squared_times_V2, A)), fq_mul(V_cubed, U2));
  // newz = V_cubed * W
  uint256 newz = fq_mul(V_cubed, W);

  r[0] = newx;
  r[1] = newy;
  r[2] = newz;
}

// Elliptic curve point multiplication.
void mul(const uint256 pt[3], const uint256 &n, uint256 r[3]) {
  if (n == 0) {
    r[0] = 1;
    r[1] = 1;
    r[2] = 0;
  } else if (n == 1) {
    cp3(pt, r);
  } else if (n & 1) {
    uint256 t[2][3];
    doubl2(pt, t[0]);
    mul(t[0], n >> 1, t[1]);
    add(t[1], pt, r);
  } else {
    uint256 t[3];
    doubl2(pt, t);
    mul(t, n >> 1, r);
  }
}

} // namespace g1

namespace g2 {

// Check if a point is the point at infinity
inline bool is_inf(const uint256 pt[3][2]) {
  return eq2(pt[2], FQ2_ZERO);
}

// Check that a point is on the curve defined by y**2 == x**3 + b.
bool is_on_curve(const uint256 pt[3][2]) {
  if (is_inf(pt)) {
    return true;
  }
  auto x = pt[0], y = pt[1], z = pt[2];
  uint256 t[3][2];
  fq2_mul(y, y, t[0]);
  fq2_mul(t[0], z, t[1]);
  fq2_mul(x, x, t[0]);
  fq2_mul(t[0], x, t[2]);
  fq2_sub(t[1], t[2], t[0]);

  fq2_mul(z, z, t[1]);
  fq2_mul(t[1], z, t[2]);
  fq2_mul(t[2], B2, t[1]);

  return eq2(t[0], t[1]);
}

bool eq(const uint256 x[3][2], const uint256 y[3][2]) {
  auto x1 = x[0], y1 = x[1], z1 = x[2];
  auto x2 = y[0], y2 = y[1], z2 = y[2];
  uint256 t[2][2];
  fq2_mul(x1, z2, t[0]);
  fq2_mul(x2, z1, t[1]);
  if (!eq2(t[0], t[1])) {
    return 0;
  }
  fq2_mul(y1, z2, t[0]);
  fq2_mul(y2, z1, t[1]);
  if (!eq2(t[0], t[1])) {
    return 0;
  }
  return 1;
}

void from_jacobian(const uint256 pt1[3][2], uint256 pt2[2][2]) {
  uint256 invz[2];
  fq2_inv(pt1[2], invz);
  fq2_mul(pt1[0], invz, pt2[0]);
  fq2_mul(pt1[1], invz, pt2[1]);
}

void doubl2(const uint256 pt[3][2], uint256 r[3][2]) {
  uint256 temp1[3][2];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      temp1[i][j] = pt[i][j];
      r[i][j] = 0;
    }
  }
  fq2_muc(temp1[0], 3, r[0]);
  fq2_mul(r[0], temp1[0], r[0]);
  fq2_mul(temp1[1], temp1[2], temp1[2]);
  fq2_mul(temp1[0], temp1[1], r[1]);
  fq2_mul(r[1], temp1[2], r[1]);
  fq2_mul(r[0], r[0], temp1[0]);
  fq2_muc(r[1], 8, r[2]);
  fq2_sub(temp1[0], r[2], temp1[0]);
  fq2_mul(temp1[2], temp1[2], r[2]);
  fq2_muc(r[1], 4, r[1]);
  fq2_sub(r[1], temp1[0], r[1]);
  fq2_mul(r[1], r[0], r[1]);
  fq2_muc(temp1[1], 8, r[0]);
  fq2_mul(r[0], temp1[1], r[0]);
  fq2_mul(r[0], r[2], r[0]);
  fq2_sub(r[1], r[0], r[1]);
  fq2_muc(temp1[0], 2, r[0]);
  fq2_mul(r[0], temp1[2], r[0]);
  fq2_mul(temp1[2], r[2], r[2]);
  fq2_muc(r[2], 8, r[2]);
}

void add(const uint256 p1[3][2], const uint256 p2[3][2], uint256 r[3][2]) {

  if (p1[2][0] == 0 && p1[2][1] == 0) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        r[i][j] = p2[i][j];
      }
    }
    return;
  } else if (p2[2][0] == 0 && p2[2][1] == 0) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        r[i][j] = p1[i][j];
      }
    }
    return;
  }

  uint256 temp1[3][2], temp2[3][2];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      temp1[i][j] = p1[i][j];
      temp2[i][j] = p2[i][j];
    }
  }

  fq2_mul(temp2[1], temp1[2], temp2[1]);
  fq2_mul(temp1[1], temp2[2], r[1]);
  fq2_mul(temp2[0], temp1[2], temp2[0]);
  fq2_mul(temp1[0], temp2[2], r[2]);

  if (temp2[0][0] == r[2][0] && temp2[0][1] == r[2][1]) {
    if (temp2[1][0] == r[1][0] && temp2[1][1] == r[1][1]) {
      doubl2(temp1, r);
      return;
    }

    r[0][0] = 1;
    r[0][1] = 0;
    r[1][0] = 1;
    r[1][1] = 0;
    r[2][0] = 0;
    r[2][1] = 0;
    return;
  }

  // W = z1 * z2
  fq2_mul(temp1[2], temp2[2], temp2[2]);
  // U = U1 - U2
  fq2_sub(temp2[1], r[1], temp1[0]);
  // V = V1 - V2
  fq2_sub(temp2[0], r[2], temp1[1]);
  // V_squared = V * V
  fq2_mul(temp1[1], temp1[1], temp1[2]);
  // V_squared_times_V2 = V_squared * V2
  fq2_mul(temp1[2], r[2], temp2[1]);
  // V_cubed = V * V_squared
  fq2_mul(temp1[2], temp1[1], temp1[2]);
  // newz = V_cubed * W
  fq2_mul(temp1[2], temp2[2], r[2]);
  // U * U
  fq2_mul(temp1[0], temp1[0], temp2[0]);
  // U * U * W
  fq2_mul(temp2[0], temp2[2], temp2[0]);
  // U * U * U - V_cubed
  fq2_sub(temp2[0], temp1[2], temp2[0]);
  // 2 * V_squared_times_V2
  fq2_muc(temp2[1], 2, temp2[2]);
  // A = U * U * W - V_cubed - 2 * V_squared_times_V2
  fq2_sub(temp2[0], temp2[2], temp2[0]);
  // newx = V * A
  fq2_mul(temp1[1], temp2[0], r[0]);
  // V_squared_times_V2 - A
  fq2_sub(temp2[1], temp2[0], temp1[1]);
  // U * (V_squared_times_V2 - A)
  fq2_mul(temp1[0], temp1[1], temp1[1]);
  // V_cubed * U2
  fq2_mul(temp1[2], r[1], temp1[0]);
  // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
  fq2_sub(temp1[1], temp1[0], r[1]);
}

// Elliptic curve point multiplication.
void mul(const uint256 pt[3][2], const uint256 &n, uint256 r[3][2]) {
  if (n == 0) {
    r[0][0] = 1;
    r[0][1] = 0;
    r[1][0] = 1;
    r[1][1] = 0;
    r[1][0] = 0;
    r[1][1] = 0;
  } else if (n == 1) {
    r[0][0] = pt[0][0];
    r[0][1] = pt[0][1];
    r[1][0] = pt[1][0];
    r[1][1] = pt[1][1];
    r[2][0] = pt[2][0];
    r[2][1] = pt[2][1];
  } else if (n & 1) {
    uint256 t[2][3][2];
    doubl2(pt, t[0]);
    mul(t[0], n >> 1, t[1]);
    add(t[1], pt, r);
  } else {
    uint256 t[3][2];
    doubl2(pt, t);
    mul(t, n >> 1, r);
  }
}

// Elliptic curve point multiplication.
void twist(const uint256 pt[3][2], uint256 r[3][12]) {
  auto _x = pt[0], _y = pt[1], _z = pt[2];
  uint256 xcoeffs[2] = {fq_sub(_x[0], fq_mul(_x[1], 9)), _x[1]};
  uint256 ycoeffs[2] = {fq_sub(_y[0], fq_mul(_y[1], 9)), _y[1]};
  uint256 zcoeffs[2] = {fq_sub(_z[0], fq_mul(_z[1], 9)), _z[1]};

  uint256 x[2], y[2], z[2];
  fq2_muc(_y, 9, x);
  fq2_sub(_x, x, x);
  cp2(_y, y);
  cp2(_z, z);

  uint256 nx[12] = {xcoeffs[0], 0, 0, 0, 0, 0, xcoeffs[1], 0, 0, 0, 0, 0};
  uint256 ny[12] = {ycoeffs[0], 0, 0, 0, 0, 0, ycoeffs[1], 0, 0, 0, 0, 0};
  uint256 nz[12] = {zcoeffs[0], 0, 0, 0, 0, 0, zcoeffs[1], 0, 0, 0, 0, 0};

  fq12_mul(W, W, r[2]);
  fq12_mul(nx, r[2], r[0]);
  fq12_mul(W, W, r[1]);
  fq12_mul(W, r[1], r[2]);
  fq12_mul(ny, r[2], r[1]);
  cp12(nz, r[2]);
}

} // namespace g2

namespace g12 {

// Check if a point is the point at infinity
inline bool is_inf(const uint256 pt[3][12]) {
  return eq12(pt[2], FQ12_ZERO);
}

// Check that a point is on the curve defined by y**2 == x**3 + b.
bool is_on_curve(const uint256 pt[3][12]) {
  if (is_inf(pt)) {
    return true;
  }
  auto x = pt[0], y = pt[1], z = pt[2];
  uint256 t[3][12];
  fq12_mul(y, y, t[0]);
  fq12_mul(t[0], z, t[1]);
  fq12_mul(x, x, t[0]);
  fq12_mul(t[0], x, t[2]);
  fq12_sub(t[1], t[2], t[0]);

  fq12_mul(z, z, t[1]);
  fq12_mul(t[1], z, t[2]);
  fq12_mul(t[2], B12, t[1]);

  return eq12(t[0], t[1]);
}

bool eq(const uint256 x[3][12], const uint256 y[3][12]) {
  auto x1 = x[0], y1 = x[1], z1 = x[2];
  auto x2 = y[0], y2 = y[1], z2 = y[2];
  uint256 t[2][12];
  fq12_mul(x1, z2, t[0]);
  fq12_mul(x2, z1, t[1]);
  if (!eq12(t[0], t[1])) {
    return 0;
  }
  fq12_mul(y1, z2, t[0]);
  fq12_mul(y2, z1, t[1]);
  if (!eq12(t[0], t[1])) {
    return 0;
  }
  return 1;
}

void doubl2(const uint256 pt[3][12], uint256 r[3][12]) {
  uint256 temp1[3][12];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 12; j++) {
      temp1[i][j] = pt[i][j];
      r[i][j] = 0;
    }
  }
  fq12_muc(temp1[0], 3, r[0]);
  fq12_mul(r[0], temp1[0], r[0]);
  fq12_mul(temp1[1], temp1[2], temp1[2]);
  fq12_mul(temp1[0], temp1[1], r[1]);
  fq12_mul(r[1], temp1[2], r[1]);
  fq12_mul(r[0], r[0], temp1[0]);
  fq12_muc(r[1], 8, r[2]);
  fq12_sub(temp1[0], r[2], temp1[0]);
  fq12_mul(temp1[2], temp1[2], r[2]);
  fq12_muc(r[1], 4, r[1]);
  fq12_sub(r[1], temp1[0], r[1]);
  fq12_mul(r[1], r[0], r[1]);
  fq12_muc(temp1[1], 8, r[0]);
  fq12_mul(r[0], temp1[1], r[0]);
  fq12_mul(r[0], r[2], r[0]);
  fq12_sub(r[1], r[0], r[1]);
  fq12_muc(temp1[0], 2, r[0]);
  fq12_mul(r[0], temp1[2], r[0]);
  fq12_mul(temp1[2], r[2], r[2]);
  fq12_muc(r[2], 8, r[2]);
}

void add(const uint256 p1[3][12], const uint256 p2[3][12], uint256 r[3][12]) {
  if (is_inf(p1)) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 12; j++) {
        r[i][j] = p2[i][j];
      }
    }
    return;
  } else if (is_inf(p2)) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        r[i][j] = p1[i][j];
      }
    }
    return;
  }

  uint256 temp1[3][12], temp2[3][12];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 12; j++) {
      temp1[i][j] = p1[i][j];
      temp2[i][j] = p2[i][j];
    }
  }

  fq12_mul(temp2[1], temp1[2], temp2[1]);
  fq12_mul(temp1[1], temp2[2], r[1]);
  fq12_mul(temp2[0], temp1[2], temp2[0]);
  fq12_mul(temp1[0], temp2[2], r[2]);

  if (eq12(temp2[0], r[2])) {
    if (eq12(temp2[1],r[1])) {
      doubl2(temp1, r);
      return;
    }

    cp12(FQ12_ONE, r[0]);
    cp12(FQ12_ONE, r[1]);
    cp12(FQ12_ZERO, r[2]);
    return;
  }

  // W = z1 * z2
  fq12_mul(temp1[2], temp2[2], temp2[2]);
  // U = U1 - U2
  fq12_sub(temp2[1], r[1], temp1[0]);
  // V = V1 - V2
  fq12_sub(temp2[0], r[2], temp1[1]);
  // V_squared = V * V
  fq12_mul(temp1[1], temp1[1], temp1[2]);
  // V_squared_times_V2 = V_squared * V2
  fq12_mul(temp1[2], r[2], temp2[1]);
  // V_cubed = V * V_squared
  fq12_mul(temp1[2], temp1[1], temp1[2]);
  // newz = V_cubed * W
  fq12_mul(temp1[2], temp2[2], r[2]);
  // U * U
  fq12_mul(temp1[0], temp1[0], temp2[0]);
  // U * U * W
  fq12_mul(temp2[0], temp2[2], temp2[0]);
  // U * U * U - V_cubed
  fq12_sub(temp2[0], temp1[2], temp2[0]);
  // 2 * V_squared_times_V2
  fq12_muc(temp2[1], 2, temp2[2]);
  // A = U * U * W - V_cubed - 2 * V_squared_times_V2
  fq12_sub(temp2[0], temp2[2], temp2[0]);
  // newx = V * A
  fq12_mul(temp1[1], temp2[0], r[0]);
  // V_squared_times_V2 - A
  fq12_sub(temp2[1], temp2[0], temp1[1]);
  // U * (V_squared_times_V2 - A)
  fq12_mul(temp1[0], temp1[1], temp1[1]);
  // V_cubed * U2
  fq12_mul(temp1[2], r[1], temp1[0]);
  // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
  fq12_sub(temp1[1], temp1[0], r[1]);
}

// Elliptic curve point multiplication.
void mul(const uint256 pt[3][12], const uint256 &n, uint256 r[3][12]) {
  if (n == 0) {
    cp12(FQ12_ONE, r[0]);
    cp12(FQ12_ONE, r[1]);
    cp12(FQ12_ZERO, r[2]);
  } else if (n == 1) {
    cp12(pt[0], r[0]);
    cp12(pt[1], r[1]);
    cp12(pt[2], r[2]);
  } else if (n & 1) {
    uint256 t[2][3][12];
    doubl2(pt, t[0]);
    mul(t[0], n >> 1, t[1]);
    add(t[1], pt, r);
  } else {
    uint256 t[3][12];
    doubl2(pt, t);
    mul(t, n >> 1, r);
  }
}

} // namespace g12

constexpr uint256 ATE_LOOP_COUNT =
    intx::from_string<uint256>("0x19d797039be763ba8");
constexpr int LOG_ATE_LOOP_COUNT = 63;

// Create a function representing the line between P1 and P2, and evaluate it at
// T
void linefunc(const uint256 p1[3], const uint256 p2[3],
                 const uint256 pt[3], uint256 r[2]) {
  uint256 x1 = p1[0], y1 = p1[1], z1 = p1[2];
  uint256 x2 = p2[0], y2 = p2[1], z2 = p2[2];
  uint256 xt = pt[0], yt = pt[1], zt = pt[2];

  uint256 m_numerator = fq_sub(fq_mul(y2, z1), fq_mul(y1, z2));
  uint256 m_denominator = fq_sub(fq_mul(x2, z1), fq_mul(x1, z2));

  if (m_denominator != 0) {
    r[0] = fq_sub(fq_mul(m_numerator, fq_sub(fq_mul(xt, z1), fq_mul(x1, zt))), fq_mul(m_denominator, fq_sub(fq_mul(yt, z1), fq_mul(y1, zt))));
    r[1] = fq_mul(fq_mul(m_denominator, zt), z1);
  } else if (m_numerator == 0) {
    m_numerator = fq_mul(fq_mul(3, x1), x1);
    m_denominator = fq_mul(fq_mul(2, y1), z1);
    r[0] = fq_sub(fq_mul(m_numerator, fq_sub(fq_mul(xt, z1), fq_mul(x1, zt))), fq_mul(m_denominator, fq_sub(fq_mul(yt, z1), fq_mul(y1, zt))));
    r[1] = fq_mul(fq_mul(m_denominator, zt), z1);
  } else {
    r[0] = fq_sub(fq_mul(xt, z1), fq_mul(x1, zt));
    r[1] = fq_mul(z1, zt);
  }
}

// Create a function representing the line between P1 and P2, and evaluate it at
// T
void linefunc12(const uint256 p1[3][12], const uint256 p2[3][12],
                const uint256 pt[3][12], uint256 r[2][12]) {
    auto x1 = p1[0], y1 = p1[1], z1 = p1[2];
    auto x2 = p2[0], y2 = p2[1], z2 = p2[2];
    auto xt = pt[0], yt = pt[1], zt = pt[2];
    uint256 m_numerator[12] = {};
    uint256 m_denominator[12] = {};
    uint256 temp[4][12] = {};
    fq12_mul(y2, z1, temp[0]);
    fq12_mul(y1, z2, temp[1]);
    fq12_sub(temp[0], temp[1], m_numerator);
    fq12_mul(x2, z1, temp[0]);
    fq12_mul(x1, z2, temp[1]);
    fq12_sub(temp[0], temp[1], m_denominator);
    if (!eq12(m_denominator, FQ12_ZERO)) {
      fq12_mul(xt, z1, temp[0]);
      fq12_mul(x1, zt, temp[1]);
      fq12_sub(temp[0], temp[1], temp[2]);
      fq12_mul(m_numerator, temp[2], temp[3]);
      fq12_mul(yt, z1, temp[0]);
      fq12_mul(y1, zt, temp[1]);
      fq12_sub(temp[0], temp[1], temp[2]);
      fq12_mul(m_denominator, temp[2], temp[0]);
      fq12_sub(temp[3], temp[0], r[0]);
      fq12_mul(zt, z1, temp[0]);
      fq12_mul(temp[0], m_denominator, r[1]);
    } else if (eq12(m_numerator, FQ12_ZERO)) {
        fq12_mul(x1, x1, temp[0]);
        fq12_muc(temp[0], 3, m_numerator);
        fq12_mul(y1, z1, temp[0]);
        fq12_muc(temp[0], 2, m_denominator);
        fq12_mul(xt, z1, temp[0]);
        fq12_mul(x1, zt, temp[1]);
        fq12_sub(temp[0], temp[1], temp[2]);
        fq12_mul(m_numerator, temp[2], temp[3]);
        fq12_mul(yt, z1, temp[0]);
        fq12_mul(y1, zt, temp[1]);
        fq12_sub(temp[0], temp[1], temp[2]);
        fq12_mul(m_denominator, temp[2], temp[0]);
        fq12_sub(temp[3], temp[0], r[0]);
        fq12_mul(zt, z1, temp[0]);
        fq12_mul(temp[0], m_denominator, r[1]);
    } else {
      fq12_mul(xt, z1, temp[0]);
      fq12_mul(x1, zt, temp[1]);
      fq12_sub(temp[0], temp[1], r[0]);
      fq12_mul(z1, zt, r[1]);
    }
}

constexpr int PSEUDO_BINARY_ENCODING[65] = {0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0,
                          0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, 1,
                          1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1,
                          1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, 1, 1};

void final_exponentiate(const uint256 x[12], const intx::uint<4096> &y,
                        uint256 r[12]) {
  uint256 o[12] = {1, 0};
  uint256 t[12] = {};
  cp12(x, t);
  intx::uint<4096> other = y;
  while (other > 0) {
    if (other & 1) {
      fq12_mul(o, t, o);
    }
    other >>= 1;
    fq12_mul(t, t, t);
  }
  cp12(o, r);
}

// Main miller loop
void _miller_loop(const uint256 Q[3][12], const uint256 P[3][12], uint256 r[12]) {
  if (g12::is_inf(Q) || g12::is_inf(P)) {
    cp12(FQ12_ONE, r);
    return;
  }

  uint256 R[3][12] = {};
  cp12(Q[0], R[0]);
  cp12(Q[1], R[1]);
  cp12(Q[2], R[2]);

  uint256 f_num[12] = {};
  uint256 f_den[12] = {};
  cp12(FQ12_ONE, f_num);
  cp12(FQ12_ONE, f_den);

  uint256 nd1[2][12] = {};
  uint256 nd2[2][12] = {};

  uint256 temp[3][12] = {};

  for (int i = 63; i > -1; i--) {
    int b = PSEUDO_BINARY_ENCODING[i];

    linefunc12(R, R, P, nd1);
    fq12_mul(f_num, f_num, nd2[0]);
    fq12_mul(nd2[0], nd1[0], f_num);
    fq12_mul(f_den, f_den, nd2[0]);
    fq12_mul(nd2[0], nd1[1], f_den);

    g12::doubl2(R, temp);
    cp12(temp[0], R[0]);
    cp12(temp[1], R[1]);
    cp12(temp[2], R[2]);

    if (b == 1) {
      linefunc12(R, Q, P, nd1);

      fq12_mul(f_num, nd1[0], nd2[0]);
      cp12(nd2[0], f_num);

      fq12_mul(f_den, nd1[1], nd2[0]);
      cp12(nd2[0], f_den);

      g12::add(R, Q, temp);
      cp12(temp[0], R[0]);
      cp12(temp[1], R[1]);
      cp12(temp[2], R[2]);
    } else if (b == -1) {
      uint256 nQ[3][12] = {};
      cp12(Q[0], nQ[0]);
      fq12_neg(Q[1], nQ[1]);
      cp12(Q[2], nQ[2]);

      linefunc12(R, nQ, P, nd1);
      fq12_mul(f_num, nd1[0], nd2[0]);
      cp12(nd2[0], f_num);

      fq12_mul(f_den, nd1[1], nd2[0]);
      cp12(nd2[0], f_den);

      g12::add(R, nQ, temp);
      cp12(temp[0], R[0]);
      cp12(temp[1], R[1]);
      cp12(temp[2], R[2]);
    }
  }

  uint256 Q1[3][12] = {};
  fq12_pow(Q[0], FIELD_MODULUS, Q1[0]);
  fq12_pow(Q[1], FIELD_MODULUS, Q1[1]);
  fq12_pow(Q[2], FIELD_MODULUS, Q1[2]);

  uint256 nQ2[3][12] = {};
  fq12_pow(Q1[0], FIELD_MODULUS, nQ2[0]);
  fq12_pow(Q1[1], FIELD_MODULUS, nQ2[2]);
  fq12_neg(nQ2[2], nQ2[1]);
  fq12_pow(Q1[2], FIELD_MODULUS, nQ2[2]);

  linefunc12(R, Q1, P, nd1);
  g12::add(R, Q1, temp);
  cp12(temp[0], R[0]);
  cp12(temp[1], R[1]);
  cp12(temp[2], R[2]);

  linefunc12(R, nQ2, P, nd2);
  fq12_mul(f_den, nd1[1], temp[0]);
  fq12_mul(temp[0], nd2[1], temp[1]);
  fq12_div(nd2[0], temp[1], temp[0]);
  fq12_mul(temp[0], nd1[0], temp[1]);
  fq12_mul(temp[1], f_num, temp[0]);

  intx::uint<4096> n =
      (intx::uint<4096>{FIELD_MODULUS} * intx::uint<4096>{FIELD_MODULUS} *
           intx::uint<4096>{FIELD_MODULUS} * intx::uint<4096>{FIELD_MODULUS} *
           intx::uint<4096>{FIELD_MODULUS} * intx::uint<4096>{FIELD_MODULUS} *
           intx::uint<4096>{FIELD_MODULUS} * intx::uint<4096>{FIELD_MODULUS} *
           intx::uint<4096>{FIELD_MODULUS} * intx::uint<4096>{FIELD_MODULUS} *
           intx::uint<4096>{FIELD_MODULUS} * intx::uint<4096>{FIELD_MODULUS} -
       1) /
      intx::uint<4096>{CURVE_ORDER};
  final_exponentiate(temp[0], n, r);
}

// Pairing computation
void _pairing(const uint256 Q[3][2], const uint256 P[3], uint256 r[12]) {
  if (g2::is_inf(Q) || g1::is_inf(P)) {
    cp12(FQ12_ONE, r);
    return;
  }

  uint256 twist_q[3][12];
  g2::twist(Q, twist_q);
  uint256 fq12_p[3][12] = {{P[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {P[1], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {P[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           };
  _miller_loop(twist_q, fq12_p, r);
}

void add(const uint256 p1[2][2], const uint256 p2[2][2], uint256 r[2][2]) {
  uint256 x[3][2] = {{p1[0][0], p1[0][1]}, {p1[1][0], p1[1][1]}, {1, 0}};
  uint256 y[3][2] = {{p2[0][0], p2[0][1]}, {p2[1][0], p2[1][1]}, {1, 0}};
  uint256 o[3][2] = {};
  g2::add(x, y, o);
  g2::from_jacobian(o, r);
}

void mul(const uint256 pt[2][2], const uint256 &n, uint256 r[2][2]) {
  uint256 x[3][2] = {{pt[0][0], pt[0][1]}, {pt[1][0], pt[1][1]}, {1, 0}};
  uint256 o[3][2] = {};
  g2::mul(x, n, o);
  g2::from_jacobian(o, r);
}

void pairing(const uint256 Q[2][2], const uint256 P[2], uint256 r[12]) {
  uint256 x[3][2] = {{Q[0][0], Q[0][1]}, {Q[1][0], Q[1][1]}, {1, 0}};
  uint256 y[3] = {P[0], P[1], 1};
  _pairing(x, y, r);
}

} // namespace bn128

#endif /* BN128_H_ */
