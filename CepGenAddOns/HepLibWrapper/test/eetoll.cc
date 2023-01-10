#include <HepLib.h>

int main() {
  Symbol me(" me "), mm(" mm "), e(" e ");
  Index mu(" mu "), nu(" nu ");
  Vector p(" p "), P(" P "), k(" k "), K(" K "), q(" q ");
  letSP(p) = me * me;
  letSP(P) = me * me;
  letSP(k) = mm * mm;
  letSP(K) = mm * mm;
#define gpm(p, m) (GAS(p) + m * GAS(1))
  ex tr1 = TR(gpm(P, -me) * GAS(mu) * gpm(p, me) * GAS(nu));
  ex tr2 = TR(gpm(k, mm) * GAS(mu) * gpm(K, -mm) * GAS(nu));
  ex res = pow(e, 4) / (4 * pow(SP(q), 2)) * tr1 * tr2;
  f or m_ us i ng _d im 4 = true;  // using 4 dimension
  res = form(res);                 // using form ( res , n ) with n >1 to print FORM script
  res = exfactor(res);
  hout << res.subs(me == 0) << endl;  // note hout is used
}
