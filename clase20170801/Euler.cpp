/* 

 Euler.cpp : output are three columns [t,f,f_exacta]

 */

#include <iostream>
#include <cmath>

using namespace std;

double f(double x, double t) {
  return x;
}

double fExacta(double t) {
  return exp(t);
}

void UnPasoDeEuler(double &x, double &t, double dt) {
  double dx;
  dx = dt*f(x,t);
  t+=dt;x+=dx;
}

int main(void) {
  double x=1,t=0;
  double dt = 0.01;
  for (t=0;t<10;){
    cout << t << " " << x << " " << fExacta(t) << endl;
    UnPasoDeEuler(x,t,dt);
  }
  return 0;
}
