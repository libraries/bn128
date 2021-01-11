#include <intx/intx.hpp>

namespace bn128 {

using uint256 = intx::uint256;

inline bool eq1(const uint256 &x, const uint256 &y) {
  return x == y;
}

inline bool eq2(const uint256 x[2], const uint256 y[2]) {
  return x[0] == y[0] && x[1] == y[1];
}

// The prime modulus of the field.
constexpr uint256 FIELD_MODULUS = intx::from_string<uint256>("0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47");

inline uint256 _addmod(const uint256 &x, const uint256 &y, const uint256 &n) {
  return n != 0 ? intx::addmod(x, y, n) : 0;
}

inline uint256 _submod(const uint256 &x, const uint256 &y, const uint256 &n) {
  return _addmod(x, n - y, n);
}

inline uint256 _negmod(const uint256 &x, const uint256 &n) {
  return _submod(n, x, n);
}

inline uint256 _mulmod(const uint256 &x, const uint256 &y, const uint256 &n) {
  return n != 0 ? intx::mulmod(x, y, n) : 0;
}

// Extended euclidean algorithm to find modular inverses for integers.
inline uint256 _invmod(const uint256 &x, const uint256 &n) {
  uint256 t = 0;
  uint256 newt = 1;
  uint256 r = n;
  uint256 newr = x;
  while (newr != 0) {
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
  } else if (y % 2 == 0) {
    return _powmod(_mulmod(x, x, n), y >> 1, n);
  } else {
    return _mulmod(_powmod(_mulmod(x, x, n), y >> 1, n), x, n);
  }
}

inline uint256 fq_add(const uint256 &x, const uint256 &y) {
  return _addmod(x, y, FIELD_MODULUS);
}

inline uint256 fq_sub(const uint256 &x, const uint256 &y) {
  return _submod(x, y, FIELD_MODULUS);
}

inline uint256 fq_neg(const uint256 &x) {
  return _negmod(x, FIELD_MODULUS);
}

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

void fq2_pow(const uint256 x[2], const uint256 &y, uint256 r[2]) {
  if (y == 0) {
    r[0] = 1;
    r[1] = 0;
  } else if (y == 1) {
    r[0] = x[0];
    r[1] = x[1];
  } else if (y % 2 == 0) {
    uint256 t[2];
    fq2_mul(x, x, t);
    fq2_pow(t, y >> 1, r);
  } else {
    uint256 t[2];
    fq2_mul(x, x, t);
    fq2_pow(t, y >> 1, t);
    fq2_mul(x, t, r);
  }
}

// The 12th-degree extension field.

// The modulus of the polynomial in this representation of FQ12.
constexpr uint256 FQ12_MODULUS_COEFFS[12] = {82, 0, 0, 0, 0, 0, FIELD_MODULUS - 18, 0, 0, 0, 0, 0};

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
  const int degree = 12;
  uint256 b[2 * degree - 1] = {0};
  for (int i = 0; i < degree; i++) {
    for (int j = 0; j < degree; j++) {
      b[i + j] = fq_add(b[i + j], fq_mul(x[i], y[j]));
    }
  }
  int lenb = 2 * degree - 1;
  while (lenb > degree) {
    int exp = lenb - degree - 1;
    uint256 top = b[lenb - 1];
    lenb-=1;
    for (int i =0; i < degree; i++) {
      b[exp + i] = fq_sub(b[exp + i], fq_mul(top, FQ12_MODULUS_COEFFS[i]));
    }
  }
  for (int i = 0; i < degree; i++) {
    r[i] = b[i];
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

void _poly_rounded_div(const uint256 a[13], const uint256 b[13], uint256 r[13]) {
  int dega = _deg(a);
  int degb = _deg(b);
  uint256 temp[13] = {};
  for (int i = 0; i < 13; i++) {
    temp[i] = a[i];
  }
  uint256 o[13] = { 0 };
  for (int i = dega - degb; i > -1; i--) {
     o[i] = fq_add(o[i], fq_div(temp[degb + i], b[degb]));
     for (int c = 0; c < degb + 1; c++) {
       temp[c + i] = fq_sub(temp[c + i], o[c]);
     }
  }
  int n = _deg(o) + 1;
  for (int i = 0; i < n; i++) {
    r[i] = o[i];
  }
}

void fq12_inv(const uint256 x[12], uint256 r[12]) {
  const int degree = 12;
  uint256 lm[degree + 1] = { lm[0]=1, 0 };
  uint256 hm[degree + 1] = { 0 };
  uint256 low[degree + 1] = {};
  for (int i = 0; i < degree ; i++) {
    low[i] = x[i];
  }
  low[degree] = 0;
  uint256 high[degree + 1] = {};
  for (int i = 0; i < degree ; i++) {
    high[i] = FQ12_MODULUS_COEFFS[i];
  }
  high[degree] = 1;

  uint256 temp[degree + 1] = {};
  uint256 nm[degree + 1] = {};
  uint256 news[degree + 1] = {};
  while (_deg(low) != 0) {
    for (int i = 0; i < degree + 1; i ++) {
      temp[i] = 0;
      nm[i] = hm[i];
      news[i] = high[i];
    }
    _poly_rounded_div(high, low, temp);
    for (int i = 0; i < degree + 1; i++) {
      for (int j = 0; j < degree + 1 - i; j++) {
        nm[i + j] = fq_sub(nm[i + j], fq_mul(lm[i], temp[j]));
        news[i + j] = fq_sub(news[i + j], fq_mul(low[i], temp[j]));
      }
    }
    for (int i = 0; i < degree + 1; i++) {
      hm[i] = lm[i];
      lm[i] = nm[i];
      high[i] = low[i];
      low[i] = news[i];
    }
  }
  for (int i = 0; i < degree; i++) {
    r[i] = fq_div(lm[i], low[0]);
  }
}

void fq12_div(const uint256 x[12], const uint256 y[12], uint256 r[12]) {
  uint256 temp[12];
  fq12_inv(y, temp);
  fq12_mul(x, temp, r);
}

void fq12_pow(const uint256 x[12], const uint256& y, uint256 r[12]) {
  if (y == 0) {
    r[0] = 1;
    for (int i = 1; i < 12; i++) {
      r[i] = 0;
    }
  } else if (y == 1) {
    for (int i = 0; i < 12; i++) {
      r[i] = x[i];
    }
  } else if (y % 2 == 0) {
    uint256 t[12];
    fq12_mul(x, x, t);
    fq12_pow(t, y >> 1, r);
  } else {
    uint256 t[12];
    fq12_mul(x, x, t);
    fq12_pow(t, y >> 1, r);
    fq12_mul(x, r, r);
  }
}

constexpr uint256 CURVE_ORDER = intx::from_string<uint256>("0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001");
// Curve is y**2 = x**3 + 3
constexpr uint256 B = 3;
// Twisted curve over FQ**2
constexpr uint256 B2[2] = {
  intx::from_string<uint256>("0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5"),
  intx::from_string<uint256>("0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2"),
};
// Extension curve over FQ**12; same b value as over FQ
constexpr uint256 B12[12] = {3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// Generator for curve over FQ
constexpr uint256 G1[2] = {1, 2};
// Generator for twisted curve over FQ2
constexpr uint256 G2[2][2] = {
  {
    intx::from_string<uint256>("0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed"),
    intx::from_string<uint256>("0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2"),
  },
  {
    intx::from_string<uint256>("0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa"),
    intx::from_string<uint256>("0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b") ,
  }
};

namespace g1 {

constexpr uint256 INF[2] = {-1, -1};

// Check if a point is the point at infinity
inline bool is_inf(const uint256 pt[2]) {
  return eq2(pt, INF);
}

// Check that a point is on the curve defined by y**2 == x**3 + b
bool is_on_curve(const uint256 pt[2]) {
    if (is_inf(pt)) {
      return true;
    }
    uint256 l = fq_mul(pt[1], pt[1]);
    uint256 r = fq_add(fq_mul(fq_mul(pt[0], pt[0]), pt[0]), B);
    return l == r;
}

void doubl2(const uint256 pt[2], uint256 r[2]) {
  uint256 x = pt[0];
  uint256 y = pt[1];
  uint256 l = fq_div(fq_mul(fq_mul(x, x), 3), fq_mul(y, 2));
  uint256 newx = fq_sub(fq_mul(l, l), fq_mul(x, 2));
  uint256 newy = fq_sub(fq_add(fq_mul(fq_sub(FIELD_MODULUS, l), newx), fq_mul(l , x)), y);
  r[0] = newx;
  r[1] = newy;
}

// Elliptic curve addition
void add(const uint256 p1[2], const uint256 p2[2], uint256 r[2]) {
  if (is_inf(p1)) {
    r[0] = p2[0];
    r[1] = p2[1];
  } else if (is_inf(p2)) {
    r[0] = p1[0];
    r[1] = p1[1];
  } else if (eq1(p2[0], p1[0]) && eq1(p2[1], p1[1])) {
    doubl2(p1, r);
  } else if (eq1(p2[0], p1[0])) {
    r[0] = INF[0];
    r[1] = INF[1];
  } else {
    uint256 l = fq_div(fq_sub(p2[1], p1[1]), fq_sub(p2[0], p1[0]));
    uint256 newx = fq_sub(fq_sub(fq_mul(l, l), p1[0]), p2[0]);
    uint256 newy = fq_sub(fq_add(fq_mul(fq_neg(l), newx),  fq_mul(l, p1[0])), p1[1]);
    r[0] = newx;
    r[1] = newy;
    return;
  }
}

// Elliptic curve point multiplication.
void mul(const uint256 pt[2], const uint256 &n, uint256 r[2]) {
  if (n == 0) {
    r[0] = INF[0];
    r[1] = INF[1];
  } else if (n == 1) {
    r[0] = pt[0];
    r[1] = pt[1];
  } else if (n % 2 == 0) {
    uint256 t[2];
    doubl2(pt, t);
    mul(t, n >> 1, r);
  } else {
    uint256 t1[2];
    uint256 t2[2];
    doubl2(pt, t1);
    mul(t1, n >> 1, t2);
    add(t2, pt, r);
  }
}

}

namespace g2 {

constexpr uint256 INF[2][2] = {
  {
    intx::from_string<uint256>("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"),
    intx::from_string<uint256>("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"),
  },
  {
    intx::from_string<uint256>("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"),
    intx::from_string<uint256>("0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff") ,
  }
};

bool is_inf(const uint256 pt[2][2]) {
  return pt[0][0] == INF[0][0] && pt[0][1] == INF[0][1] && pt[1][0] == INF[1][0] && pt[1][1] == INF[1][1];
}

// Check that a point is on the curve defined by y**2 == x**3 + b.
bool is_on_curve(const uint256 pt[2][2]) {
    if (is_inf(pt)) {
      return true;
    }
    uint256 t1[2];
    uint256 t2[2];
    fq2_mul(pt[0], pt[0], t2);
    fq2_mul(pt[0], t2, t1);
    fq2_add(t1, B2, t2);
    fq2_mul(pt[1], pt[1], t1);
    return eq2(t1, t2);
}

void doubl(const uint256 pt[2][2], uint256 r[2][2]) {
  uint256 t[2];
  uint256 newx[2];
  uint256 newy[2];
  uint256 l[2];

  uint256 u3[2] = {3, 3};
  uint256 u2[2] = {2, 2};

  fq2_mul(pt[0], pt[0], t);
  fq2_mul(u3, t, t);
  fq2_mul(u2, pt[1], l);
  fq2_div(t, l, l);

  fq2_mul(l, l, t);
  fq2_mul(u2, pt[0], newx);
  fq2_sub(t, newx, newx);

  fq2_mul(l, pt[0], t);
  fq2_sub(t, pt[1], t);
  fq2_neg(l, newy);
  fq2_mul(newy, newx, newy);
  fq2_add(newy, t, newy);

  r[0][0] = newx[0];
  r[0][1] = newx[1];
  r[1][0] = newy[0];
  r[1][1] = newy[1];
}

}

}
