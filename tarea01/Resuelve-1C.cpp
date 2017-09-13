#include <iostream>
#include <cmath>

#include "Vector.h"

using namespace std;

double beta = 0.3;
double gama = 0.05;

vector3D f(vector3D &x, double t) {
  vector3D fx;
  fx.cargue(-beta*x.x()*x.y(),
	    beta*x.x()*x.y()-gama*x.y(),
	    gama*x.y());
  return fx;
}

void UnPasoDeRungeKutta4(vector3D &x, double &t, double dt) {
  vector3D dx1,dx2,dx3,dx4,xtmp;
  dx1 = dt*f(x,t);
  xtmp = x+0.5*dx1;
  dx2 = dt*f(xtmp, t+0.5*dt);
  xtmp = x+0.5*dx2;
  dx3 = dt*f(xtmp, t+0.5*dt);
  xtmp = x+dx3;
  dx4 = dt*f(xtmp,t+dt);
  t+=dt;x+=((dx1+2*dx2+2*dx3+dx4)/6);
}

int main(void) {
  double s0 = 0.999, // susceptibles
    i0 = 0.001,      // infectados
    r0 = 0.0;        // rertirados

  vector3D x, fx;
  double t, dt = 0.1;
  double R0;
  gama = 0.5;
  for (R0=0.5;R0<=10;R0+=0.1) {
    x.cargue(s0,i0,r0);
    for (t=0;t<2000;) {
      beta = R0*gama;
      UnPasoDeRungeKutta4(x,t,dt);
    }
    cout << R0 << " " << x.x() << endl;
    t=0.0;
  }

  return 0;
}
