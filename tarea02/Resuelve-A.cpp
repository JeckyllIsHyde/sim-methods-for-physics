/*
  Se toma la masa de Jupiter como 1newKg y la del sol como 1047newKg
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"

using namespace std;

const double G = 1.0;
const int N = 2;

const double Zeta = +0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi = -0.06626458266981849;

class Cuerpo;
class Colisionador;

void InicieAnimacion(void) {
  cout << "set terminal gif animate" << endl; 
  cout << "set output 'SolYJupiter.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [-1200:1200]" << endl;
  cout << "set yrange [-1200:1200]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set trange [0:7]" << endl;
  cout << "set isosamples 12" << endl;
}

void InicioCuadro(void) {
  cout << "plot 0,0 ";
}

void TermineCuadro(void) {
  cout << endl;
}

int main(void) {
  double t, dt; // 1000
  double tdibujo; int Ndibujos;
  Cuerpo planetas[N];
  int i;
  Colisionador newton;

  double m0=1047, m1=1, r=1000;
  double R0=100, R1=10;

  double M = m0+m1;
  double x0 = -m1/M*r, x1 = x0+r;
  
  double a = r, omega=sqrt(G*M/(a*a*a)),
    vy0=omega*x0, vy1=1.0*omega*x1,
    T=2*M_PI/omega, tmax=20.*T;

  dt = 20*T/(370); // M*T = n*dt
  
  //InicieAnimacion(); Ndibujos=500;
  //                ( x0, y0, z0,Vx0,Vy0,Vz0, m0, R0 )
  planetas[0].Inicio( x0,0.0,0.0,0.0,vy0,0.0, m0, R0 );
  planetas[1].Inicio( x1,0.0,0.0,0.0,vy1,0.0, m1, R1 );
  
  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt) {
    /*if (tdibujo>tmax/Ndibujos) {
      InicioCuadro();
      for (i=0;i<N;i++)
	planetas[i].Dibujese();
      TermineCuadro();
      tdibujo = 0;
    }*/
    cout << planetas[0].Getx() << " " << planetas[0].Gety()
	 << planetas[0].Getx() << " " << planetas[0].Gety() << endl;
    // Muevase con Omelyan PEFRL
    for (i=0;i<N;i++) planetas[i].Mueva_r(dt,Zeta);
    newton.CalculeTodasLasFuerzas(planetas); for (i=0;i<N;i++) planetas[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) planetas[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(planetas); for (i=0;i<N;i++) planetas[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) planetas[i].Mueva_r(dt,1-2*(Xi+Zeta));
    newton.CalculeTodasLasFuerzas(planetas); for (i=0;i<N;i++) planetas[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) planetas[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(planetas); for (i=0;i<N;i++) planetas[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) planetas[i].Mueva_r(dt,Zeta);
  }

  return 0;
}
