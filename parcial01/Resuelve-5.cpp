#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double K = 1000;
const double Kcentral = 1;
const double GAMMA = 0.5;

const double MASA = 1;
const double RADIO = 1;
const double CARGA = 1;
const double Lx = 200, Ly = 200;
const double Ranillo = 40;

const double Zeta = +0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi = -0.06626458266981849;

class Cuerpo {
private:
  vector3D r,v,F;
  double m,R;
public:
  friend class SimuladorDeReglas;
  void Inicio(double  x0, double y0,
	      double Vx0, double Vy0,
	      double m0,double R0);
  void BorreFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Mueva_r(double dt, double Constante);
  void Mueva_v(double dt, double Constante);
  void Dibujese(void);
  double Getx(void) {return r.x();};
  double Gety(void) {return r.y();};
};

void Cuerpo::Inicio(double x0,double y0,
		    double Vx0,double Vy0,
		    double m0, double R0) {
  r.cargue(x0,y0,0.0);
  v.cargue(Vx0,Vy0,0.0);
  F.cargue(0.0,0.0,0.0);
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
  cout << ", " << r.x() << "+"<< R << "*cos(t),"
       << r.y() << "+" << R << "*sin(t) ";
}

class SimuladorDeReglas {
private:
  int N;
public:
  SimuladorDeReglas(int n) : N(n) {} ;
  void CalculeTodasLasFuerzas(Cuerpo* cuerpos);
  void CalculeLaFuerzaEntre(Cuerpo& cuerpo1, Cuerpo& cuerpo2);
  double CalculeRmaxConfig(Cuerpo* cuerpos);
};

void SimuladorDeReglas::CalculeTodasLasFuerzas(Cuerpo* cuerpos) {
  int i, j;
  for (i=0; i<N; i++) {
    cuerpos[i].BorreFuerza();
    vector3D Fcentral = -Kcentral*cuerpos[i].r;
    vector3D Fviscosa = -GAMMA*cuerpos[i].v;
    cuerpos[i].AgregueFuerza(Fcentral+Fviscosa);
  }
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeLaFuerzaEntre(cuerpos[i],cuerpos[j]);
}

void SimuladorDeReglas::CalculeLaFuerzaEntre(Cuerpo& cuerpoi, Cuerpo& cuerpoj) {
  vector3D rij = cuerpoj.r-cuerpoi.r;
  double IrijI3 = pow(norma2(rij),1.5);
  vector3D Frepulsion = (K*CARGA*CARGA/IrijI3)*rij;
  cuerpoj.AgregueFuerza(Frepulsion);
  cuerpoi.AgregueFuerza((-1)*Frepulsion);
}

double SimuladorDeReglas::CalculeRmaxConfig(Cuerpo* cuerpos) {
  int i, j;
  double rmax = 0;
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      if (norma(cuerpos[j].r-cuerpos[i].r) > rmax)
	rmax = norma(cuerpos[j].r-cuerpos[i].r);
  return rmax;
	
}

int main(void) {
  double t, dt;
  
  int N, Nmax=201;
  for (N=200;N<Nmax;N+=10) {
    Cuerpo cargas[N];
    SimuladorDeReglas newton(N);
    Crandom rand64(1);
    int i,j;
    double th=0.0,thv;

    double T=2*M_PI*sqrt(MASA/Kcentral), tmax=5*T;
    dt = tmax/5000;

    //          (     x0,   y0,Vx0,Vy0,  m0,   R0 )
    for (i=0;i<N;i++) {
      th += M_PI/N;
      thv = 2*M_PI*rand64.r();
      cargas[i].Inicio(Ranillo*cos(th),Ranillo*sin(th),
		       -sin(thv),cos(thv),MASA,RADIO );
    }
    for (t=0;t<tmax;t+=dt) {
      // Muevase con Omelyan PEFRL
      for (i=0;i<N;i++) cargas[i].Mueva_r(dt,Zeta);
      newton.CalculeTodasLasFuerzas(cargas);
      for (i=0;i<N;i++) cargas[i].Mueva_v(dt,(1-2*Lambda)/2);
      for (i=0;i<N;i++) cargas[i].Mueva_r(dt,Xi);
      newton.CalculeTodasLasFuerzas(cargas);
      for (i=0;i<N;i++) cargas[i].Mueva_v(dt,Lambda);
      for (i=0;i<N;i++) cargas[i].Mueva_r(dt,1-2*(Xi+Zeta));
      newton.CalculeTodasLasFuerzas(cargas);
      for (i=0;i<N;i++) cargas[i].Mueva_v(dt,Lambda);
      for (i=0;i<N;i++) cargas[i].Mueva_r(dt,Xi);
      newton.CalculeTodasLasFuerzas(cargas);
      for (i=0;i<N;i++) cargas[i].Mueva_v(dt,(1-2*Lambda)/2);
      for (i=0;i<N;i++) cargas[i].Mueva_r(dt,Zeta);
    }

    // Calcule rmax de la configuracion
    cout << N << " " << newton.CalculeRmaxConfig(cargas) << endl;
  }
	 
  return 0;
}
