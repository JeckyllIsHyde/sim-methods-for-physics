/* 

 Runge-Kutta4-2doOrden.cpp : 
 Output are three columns [t,f1,f2,f1_exacta,f2_exacta]

 */

#include <iostream>
#include <cmath>

using namespace std;

const double ERR = 1e-7;
const double Omega0 = 1;

double f1(double x1, double x2, double t) { return x2; }
double f2(double x1, double x2, double t) { return -Omega0*Omega0*x1;}
double fExacta(double t) { return sin(t);}

void UnPasoDeRK4_2doOrden(double &x1, double &x2, double &t, double dt) {
  double dx11,dx21,dx31,dx41,                        dx12,dx22,dx32,dx42;
  dx11 = dt*f1(x1,x2,t);                             dx12 = dt*f2(x1,x2,t);
  dx21 = dt*f1(x1+0.5*dx11,x2+0.5*dx12, t+0.5*dt);   dx22 = dt*f2(x1+0.5*dx11,x2+0.5*dx12, t+0.5*dt);
  dx31 = dt*f1(x1+0.5*dx21,x2+0.5*dx22, t+0.5*dt);   dx32 = dt*f2(x1+0.5*dx21,x2+0.5*dx22, t+0.5*dt);
  dx41 = dt*f1(x1+dx31,x2+dx32, t+dt);               dx42 = dt*f2(x1+dx31,x2+dx32, t+dt);
  t+=dt;
  x1+=((dx11+2*dx21+2*dx31+dx41)/6);                 x2+=((dx12+2*dx22+2*dx32+dx42)/6);
}

int main(void) {
  double x1=1, x2=1, t=0, dt=0.1;
  for (t=0;t<10; UnPasoDeRK4_2doOrden(x1,x2,t,dt))
    cout << t << " " << x1 << " " << x2 << " " << fExacta(t) << endl;
  return 0;
}
