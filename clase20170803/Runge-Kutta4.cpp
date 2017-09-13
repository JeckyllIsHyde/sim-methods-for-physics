#include <iostream>
#include <cmath>

using namespace std;

const double ERR = 1e-7;

double f(double x, double t) {
  return x;
}

double fExacta(double t) {
  return exp(t);
}

int main(void) {
  double x=1,t=0;
  double dt = 0.1;
  for (t=0;t<10;) {
    cout << t << " " << x << " " << fExacta(t) << endl;
    UnPasoDeRungeKutta4(x,t,dt);
  }
  return 0;
}
