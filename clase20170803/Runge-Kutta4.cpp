/* 

 Runge-Kutta4.cpp : output are three columns [t,f,f_exacta]

 */

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

void UnPasoDeRungeKutta4(double &x, double &t, double dt) {
  double dx1,dx2,dx3,dx4;
  dx1 = dt*f(x,t);
  dx2 = dt*f(x+0.5*dx1, t+0.5*dt);
  dx3 = dt*f(x+0.5*dx2, t+0.5*dt);
  dx4 = dt*f(x+dx3,t+dt);
  t+=dt;x+=((dx1+2*dx2+2*dx3+dx4)/6);
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
