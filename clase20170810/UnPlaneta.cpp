#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double GM = 1.0;

class Cuerpo;

class Cuerpo {
 private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
 public:
  void Inicio(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  void Dibujese(void) const;
  double Getx(void) const {return x;};
  double Gety(void) const {return y;};
};

void Cuerpo::Inicio(double x0,double y0,double Vx0,double Vy0,double m0,double R0) {
  x = x0; y = y0; Vx = Vx0; Vy = Vy0; m = m0; R = R0;
}

void Cuerpo::CalculeFuerza(void) {
  double aux = -GM*m*pow(x*x+y*y,-1.5);
  Fx = aux*x;Fy = aux*y;
}
                                  
void Cuerpo::Muevase(double dt) {
  x+=Vx*dt;y+=Vy*dt;
  Vx+=Fx*dt/m;Vy+=Fy*dt/m;
}

void Cuerpo::Dibujese(void) const {
  cout << ", " << x << "+"<< R << "*cos(t)," << y << "+" << R << "*sin(t)";
}

void InicieAnimacion(void) {
  cout << "set terminal gif animate" << endl; 
  cout << "set output 'MiPlaneta.gif'" << endl; 
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
