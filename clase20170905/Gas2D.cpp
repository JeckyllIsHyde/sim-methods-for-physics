#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"

using namespace std;

const double K = 1.e4;
const double Lx = 100.0, Ly = 100.0;
const int Nx = 7; int Ny = 8;
const int N = Nx*Ny;

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
  void Dibujese(void);
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

void Cuerpo::Dibujese(void) {
  cout << ", " << r.x() << "+"<< R << "*cos(t)," << r.y() << "+" << R << "*sin(t) ";
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
    for(j=i+1;j<N+4;j++)
      CalculeLaFuerzaEntre(cuerpos[i],cuerpos[j]);
}

void Colisionador::CalculeLaFuerzaEntre(Cuerpo& cuerpo1,
					Cuerpo& cuerpo2) {
  vector3D F2, dr, dr_u; double norma_dr, s;
  dr= cuerpo2.r-cuerpo1.r;
  norma_dr = norma(dr); dr_u = dr/norma_dr;
  s = (cuerpo1.R+cuerpo2.R)-norma_dr;
  if (s>0) {
    F2 = dr_u*K*pow(s,1.5);
    cuerpo2.AgregueFuerza(F2);
    cuerpo1.AgregueFuerza(F2*(-1.0));
  }
}

void InicieAnimacion(void) {
  cout << "set terminal gif animate" << endl; 
  cout << "set output 'Particulas.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [-10:110]" << endl;
  cout << "set yrange [-10:110]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set trange [0:7]" << endl;
  cout << "set isosamples 12" << endl;
}

void InicioCuadro(void) {
  cout << "plot 0,0 ";
  cout << " , " << 100./7 << "*t, 0";
  cout << " , " << 100./7 << "*t, 100";
  cout << " , 0," << 100./7 << "*t";
  cout << " , 100," << 100./7 << "*t";
}

void TermineCuadro(void) {
  cout << " title 'particulas' ";
  cout << endl;
}

int main(void) {
  double t, dt=1e-3;
  double tdibujo; int Ndibujos;
  Cuerpo granos[N+4];
  int i,j;
  Colisionador newton;

  double m0=1, R0=3, v=10;
  double Rpared=10000, Mpared=1000*m0 ;

  double T=Lx/v, tmax=5*T;
  
  InicieAnimacion(); Ndibujos=500;
  //                ( x0, y0, z0,Vx0,Vy0,Vz0, m0, R0 )
  // pared arriba
  granos[N+0].Inicio(Lx/2,Ly+Rpared,0.0,0.0,0.0,0.0, Mpared, Rpared );
  // pared abajo
  granos[N+1].Inicio(Lx/2,  -Rpared,0.0,0.0,0.0,0.0, Mpared, Rpared );
  // pared derercha
  granos[N+2].Inicio(Lx+Rpared,Ly/2,0.0,0.0,0.0,0.0, Mpared, Rpared );
  // pared izquierda
  granos[N+3].Inicio(  -Rpared,Ly/2,0.0,0.0,0.0,0.0, Mpared, Rpared );

  double  row, col;
  for (i=0; i<N; i++) {
    row = (int)i/(int)Nx; col = i-row*Nx;
    granos[i].Inicio((1+col)*Lx/(Nx+2),(1+row)*Ly/(Ny+2),0.0,0.8*v,0.6*v,0.0, m0, R0 );
  }

  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt) {
    if (tdibujo>tmax/Ndibujos) {
      InicioCuadro();
      for (i=0; i<N; i++)
	granos[i].Dibujese();
      TermineCuadro();
      tdibujo = 0;
    }
    // cout << granos[1].Getx() <<" " << granos[1].Gety() << endl; 
    // Muevase con Omelyan PEFRL
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Zeta);
    newton.CalculeTodasLasFuerzas(granos); for (i=0;i<N;i++) granos[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(granos); for (i=0;i<N;i++) granos[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,1-2*(Xi+Zeta));
    newton.CalculeTodasLasFuerzas(granos); for (i=0;i<N;i++) granos[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(granos); for (i=0;i<N;i++) granos[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Zeta);
  }
  
  return 0;
}
