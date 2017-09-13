/* 
   Solucion de una ecuacion valores de frontera  por el método
   de la lanzadera (shooting method), utilizando RK4
*/
#include <iostream>
#include <cmath>

using namespace std;

const double ERR = 1e-7;

double f1(double x1, double x2, double t) { return x2; }
double f2(double Omega0, double x1, double x2, double t) { return -Omega0*Omega0*x1;}
double fExacta(double t) { return sin(t);}

void UnPasoDeRK4_2doOrden(double w0, double &x1, double &x2, double &t, double dt) {
  double dx11,dx21,dx31,dx41,                        dx12,dx22,dx32,dx42;
  dx11 = dt*f1(x1,x2,t);                             dx12 = dt*f2(w0, x1,x2,t);
  dx21 = dt*f1(x1+0.5*dx11,x2+0.5*dx12, t+0.5*dt);   dx22 = dt*f2(w0, x1+0.5*dx11,x2+0.5*dx12, t+0.5*dt);
  dx31 = dt*f1(x1+0.5*dx21,x2+0.5*dx22, t+0.5*dt);   dx32 = dt*f2(w0, x1+0.5*dx21,x2+0.5*dx22, t+0.5*dt);
  dx41 = dt*f1(x1+dx31,x2+dx32, t+dt);               dx42 = dt*f2(w0, x1+dx31,x2+dx32, t+dt);
  t+=dt;
  x1+=((dx11+2*dx21+2*dx31+dx41)/6);                 x2+=((dx12+2*dx22+2*dx32+dx42)/6);
}

double f(double w0) {
  double x1=1, x2=1, t=0, dt=0.01;
  for (t=0;t<1; ) {
    UnPasoDeRK4_2doOrden(w0,x1,x2,t,dt);
  }
  return x1;
}

double ceroPorBiseccion(double a, double b) {
  double m, fa, fb, fm;
  fa = f(a); fb = f(b);
  do {
    m = (a+b)/2; fm = f(m);
    if (fa*fm<0)
      {b=m;fb = fm;}
    else
      {a=m;fa = fa;}
  } while (b-a > ERR);

  return (a+b)/2;
}

int main(void) {
  double w0;
  for (w0=0.1;w0<20; w0+=0.1)
    cout << f(w0) << endl;
  return 0;
}
