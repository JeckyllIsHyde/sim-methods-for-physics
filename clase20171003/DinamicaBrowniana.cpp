#include <iostream>
#include <cmath>
#include <Random64.h>

using namespace std;

const double tTotal = 400.0;
const double dt = 0.01;
const double L = 100;
const double T = 100;
const double MASA = 22.8916;
const double e=1;
const double D=0.132;
const double kb = 0.826;
const double Gamma = kb*T/(MASA*D);
const double sigma=sqrt(2*D*dt);
const double dtU2mGamma = dt/(2*MASA*Gamma);
  
class Cuerpo;

class Cuerpo {
 private:
  double x,F,Fold,m,R;
 public:
  void Inicio(double x0,double m0,double R0);
  void CalculeFuerza(double f);
  void Muevase(Crandom& ran64);
  void Dibujese(void) const;
  double Getx(void) const {return x;};
};

void Cuerpo::Inicio(double x0,double m0,double R0) {
  x = x0; m = m0; R = R0; F = 0; Fold = 0;
}

void Cuerpo::CalculeFuerza(double f) {
  Fold = F;  F = e*f;
}

void Cuerpo::Dibujese(void) const {
  cout << ", " << x << "+"<< R << "*cos(t)," << 0.0 << "+" << R << "*sin(t)";
}

void Cuerpo::Muevase( Crandom& ran64 ) {
  x+=dtU2mGamma*(3*F-Fold) + ran64.gauss(0.0,sigma);//Vx*dt;
}

void InicieAnimacion(void) {
  //cout << "set terminal gif animate" << endl; 
  //cout << "set output 'DifusionSodio.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [-10:110]" << endl;
  cout << "set yrange [-10:10]" << endl;
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
  Cuerpo Na;
  Crandom ran64(1);
  double t, tmax=40;
  double E = 0;

  
  InicieAnimacion();
  //          ( x0,   m0,  R0)
  Na.Inicio(L/2,MASA,4);
  Na.CalculeFuerza(E);
  
  for (t=0;t<tmax;t+=dt) {
    InicioCuadro();
    Na.Dibujese();
    TermineCuadro();
    Na.CalculeFuerza(E);
    Na.Muevase(ran64);
  }
  
  return 0;
}
