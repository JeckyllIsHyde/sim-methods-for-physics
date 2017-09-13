#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"

using namespace std;

const double G = 9.81;
const double K = 200e9;
const int N = 2;

const double Zeta = +0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi = -0.06626458266981849;

class Cuerpo;
class Colisionador;

class Cuerpo {
 private:
  double tau, omega, th,
    m, R, L, I, xcorrido;
 public:
  friend class Colisionador;
  void Inicio(double th0,
	      double omega0,
	      double m0, double R0, double L0,
	      double x0corrido);
  void BorreFuerza(void);
  void InicieFuerza(void);
  void AgregueFuerza(double tau0);
  void Mueva_r(double dt, double Constante);
  void Mueva_v(double dt, double Constante);
  void Dibujese(void);
  double x(void) {return xcorrido+L*sin(th);};
  double y(void) {return -L*cos(th);};
};

void Cuerpo::Inicio(double th0,
		    double omega0,
		    double m0,double R0, double L0,
		    double x0corrido) {
  th = th0;
  omega = omega0;
  R = R0; L = L0; 
  m = m0; I = m*L*L;
  xcorrido = x0corrido;
}

void Cuerpo::BorreFuerza(void) {
  tau=0.0;
}
				  
void Cuerpo::InicieFuerza(void) {
  tau=-m*G*L*sin(th);
}

void Cuerpo::AgregueFuerza(double tau0) {
  tau+=tau0;
}

void Cuerpo::Mueva_r(double dt, double Constante) {
  th+=omega*(Constante*dt);
}

void Cuerpo::Mueva_v(double dt, double Constante) {
  omega+=tau*(Constante*dt/I);
}

void Cuerpo::Dibujese(void) {
  cout << ", " << x() << "+"<< R << "*cos(t),"
       << y() << "+" << R << "*sin(t), "
       << xcorrido << "+" << (x()-xcorrido)/7.0 << "*t,"
       << y()/7.0 << "*t";
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
    cuerpos[i].AgregueFuerza(-G*cuerpos[i].m*sin(cuerpos[i].th));
  for(i=0;i<N-1;i++)
      CalculeLaFuerzaEntre(cuerpos[i],cuerpos[i+1]);
};

void Colisionador::CalculeLaFuerzaEntre(Cuerpo& cuerpo1,
					Cuerpo& cuerpo2) {
  double s,d21,F2,tau2,tau1;
  d21 = cuerpo2.x() - cuerpo1.x();
  s = cuerpo2.R + cuerpo1.R - d21;
  if (s>0) {
    F2 = K*pow(s,1.5)*d21/fabs(d21);
    tau2 = F2*cuerpo2.L;
    tau1 = -F2*cuerpo1.L;
    cuerpo2.AgregueFuerza(tau2);
    cuerpo1.AgregueFuerza(tau1);
  }
}

void InicieAnimacion(void) {
  cout << "set terminal gif animate" << endl; 
  cout << "set output 'CunaDeNewton.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [-12:22]" << endl;
  cout << "set yrange [-12:0]" << endl;
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
  double t, dt=1e-4, tmax;
  double tdibujo; int Ndibujos;
  Cuerpo pendulos[N];
  int i;
  Colisionador newton;

  double m0=1, L0=10;
  double R0=1;

  double x0corrido = 0;
  
  double th0=-15*M_PI/180;

  double T = 2*M_PI*sqrt(L0/G);
  tmax = 1.5*T;

  InicieAnimacion(); Ndibujos=500;
  //                ( th0, w0, m0, R0, L0, x0corrido )
  pendulos[0].Inicio( th0,0.0, m0, R0,L0, x0corrido );
  for (i=1;i<N;i++)
    pendulos[i].Inicio( 0.0,0.0, m0, R0,L0, x0corrido+2*R0*i );
  
  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt) {
    if (tdibujo>tmax/Ndibujos) {
      InicioCuadro();
      for (i=0;i<N;i++)
	pendulos[i].Dibujese();
      TermineCuadro();
      tdibujo = 0;
    }
    // cout << pendulos[1].Getx() <<" " << pendulos[1].Gety() << endl; 
    // Muevase con Omelyan PEFRL
    for (i=0;i<N;i++) pendulos[i].Mueva_r(dt,Zeta);
    newton.CalculeTodasLasFuerzas(pendulos); for (i=0;i<N;i++) pendulos[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) pendulos[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(pendulos); for (i=0;i<N;i++) pendulos[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) pendulos[i].Mueva_r(dt,1-2*(Xi+Zeta));
    newton.CalculeTodasLasFuerzas(pendulos); for (i=0;i<N;i++) pendulos[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) pendulos[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(pendulos); for (i=0;i<N;i++) pendulos[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) pendulos[i].Mueva_r(dt,Zeta);
  }

  return 0;
}
