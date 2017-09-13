#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"

using namespace std;

const double GM = 1.0;

class Cuerpo;

class Cuerpo {
 private:
  vector3D r,v,F;
  double m,R;
 public:
  void Inicio(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void CalculeFuerza(void);
  void Arranque(double dt);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void) {return r.x();};
  double Gety(void) {return r.y();};
};

void Cuerpo::Inicio(double x0,double y0,double z0,
		    double Vx0,double Vy0,double Vz0,
		    double m0,double R0) {
  v.cargue(Vx0,Vy0,Vz0);
  r.cargue(x0,y0,z0);
  m = m0; R = R0;
}

void Cuerpo::CalculeFuerza(void) {
  double aux = -GM*m*pow(norma2(r),-1.5);
  F = aux*r;
}
				  
void Cuerpo::Arranque(double dt) {
  v-=F*(0.5*dt/m);
}
				  
void Cuerpo::Muevase(double dt) {
  v+=F*dt/m;
  r+=v*dt;
}

void Cuerpo::Dibujese(void) {
  cout << ", " << r.x() << "+"<< R << "*cos(t)," << r.y() << "+" << R << "*sin(t)";
}

void InicieAnimacion(void) {
  cout << "set terminal gif animate" << endl; 
  cout << "set output 'MiPlanetaConVectores.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [-120:120]" << endl;
  cout << "set yrange [-120:120]" << endl;
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
  double t, tmax, dt=10.0;
  double tdibujo; int Ndibujos;
  double m=1, R=5;
  double r=100, omega=sqrt(GM/(r*r*r)),
    v=omega*r, T=2*M_PI/omega; tmax=1.1*T;
  Cuerpo planeta;


  //InicieAnimacion(); Ndibujos=500;
  //          ( x0, y0, z0,Vx0,Vy0,Vz0, m0, R0)
  planeta.Inicio(r,0.0,0.0,
		 0.0,0.5*v,0.0,
		 m,R);
  planeta.CalculeFuerza();
  planeta.Arranque(dt);
  
  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt) {
    /*if (tdibujo>tmax/Ndibujos) {
      InicioCuadro();
      planeta.Dibujese();
      TermineCuadro();
      tdibujo = 0;
    }*/
    cout<<planeta.Getx() <<" " << planeta.Gety() << endl;
    planeta.CalculeFuerza();
    planeta.Muevase(dt);
  }
  
  return 0;
}
