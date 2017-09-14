/*
  Se toma la masa de Jupiter como 1newKg y la del sol como 1047newKg
  La distancia entre ambos es 1000 newm. 
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

class Cuerpo {
 private:
  vector3D r,v,F;
  double m,R;
 public:
  friend class Colisionador;
  void Inicio(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void BorreFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Mueva_r(double dt, double Constante);
  void Mueva_v(double dt, double Constante);
  void DibujeseRelativaA(double x0, double y0);
  double Getx(void) {return r.x();};
  double Gety(void) {return r.y();};
};

void Cuerpo::Inicio(double x0,double y0,double z0,
		    double Vx0,double Vy0,double Vz0,
		    double m0,double R0) {
  r.cargue(x0,y0,z0);
  v.cargue(Vx0,Vy0,Vz0);
  m = m0; R = R0;
}

void Cuerpo::BorreFuerza(void) {
  F.cargue(0.0,0.0,0.0);
}
				  
void Cuerpo::AgregueFuerza(vector3D F0) {
  F+=F0;
}

void Cuerpo::Mueva_r(double dt, double Constante) {
  r+=v*(Constante*dt);
}

void Cuerpo::Mueva_v(double dt, double Constante) {
  v+=F*(Constante*dt/m);
}

void Cuerpo::DibujeseRelativaA(double x0, double y0) {
  cout << ", " << r.x()-x0 << "+"<< R << "*cos(t)," << r.y()-y0 << "+" << R << "*sin(t)";
}

class Colisionador {
private:
  
public:
  void CalculeTodasLasFuerzas(Cuerpo* cuerpos);
  void CalculeLaFuerzaEntre(Cuerpo& cuerpo1, Cuerpo& cuerpo2);
};

void Colisionador::CalculeTodasLasFuerzas(Cuerpo* cuerpos) {
  int i,j;
  for(i=0;i<N;i++)
    cuerpos[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeLaFuerzaEntre(cuerpos[i],cuerpos[j]);
}

void Colisionador::CalculeLaFuerzaEntre(Cuerpo& cuerpo1,
					Cuerpo& cuerpo2) {
  vector3D F,dr;
  dr= cuerpo2.r-cuerpo1.r;
  double aux = G*cuerpo1.m*cuerpo2.m*pow(norma2(dr),-1.5);
  F = dr*aux;
  cuerpo1.AgregueFuerza(F);
  cuerpo2.AgregueFuerza(F*(-1.0));
}

void InicieAnimacion(void) {
  cout << "set terminal pngcairo" << endl; 
  cout << "set output 'SolYJupiter.png'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [-1200:1200]" << endl;
  cout << "set yrange [-1200:1200]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set trange [0:7]" << endl;
  //  cout << "set isosamples 12" << endl;
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

  dt = 20*T/(150); // M*T = n*dt
  
  InicieAnimacion(); Ndibujos=500;
  //                ( x0, y0, z0,Vx0,Vy0,Vz0, m0, R0 )
  planetas[0].Inicio( x0,0.0,0.0,0.0,vy0,0.0, m0, R0 );
  planetas[1].Inicio( x1,0.0,0.0,0.0,vy1,0.0, m1, R1 );
   
  cout << " set grid back ls 12 " << endl
       << " set title 'T = " << T
       << "[newS] , dt = " << dt << "' " << endl;
  InicioCuadro();
  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt) {
    if (tdibujo>tmax/Ndibujos) {
      planetas[0].DibujeseRelativaA( (m0*planetas[0].Getx()+m1*planetas[1].Getx())/M, 
				     (m0*planetas[0].Gety()+m0*planetas[0].Gety())/M );
      cout << " lc rgb 'red' ";
      planetas[1].DibujeseRelativaA( (m0*planetas[0].Getx()+m1*planetas[1].Getx())/M, 
				     (m0*planetas[0].Gety()+m0*planetas[0].Gety())/M );
      cout << " lc rgb 'blue' ";
      tdibujo = 0;
    }
    // cout << planetas[0].Getx() << " " << planetas[0].Gety()
    //	 << planetas[0].Getx() << " " << planetas[0].Gety() << endl;
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
  TermineCuadro();

  return 0;
}
