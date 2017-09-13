#include <iostream>
#include <cmath>

using namespace std;

const double g = 9.81;

class Cuerpo;

class Cuerpo {
 private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
 public:
  void Inicio(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  double Getx(void) const {return x;};
  double Gety(void) const {return y;};
};

void Cuerpo::Inicio(double x0,double y0,double Vx0,double Vy0,double m0,double R0) {
  x = x0; y = y0; Vx = Vx0; Vy = Vy0; m = m0; R0 = R0;
}

void Cuerpo::CalculeFuerza(void) {
  Fx = 0;Fy = -m*g;
}
                                  
void Cuerpo::Muevase(double dt) {
  x+=Vx*dt;y+=Vy*dt;
  Vx+=Fx*dt/m;Vy+=Fy*dt/m;
}

int main(void) {
  double t, dt=0.01;
  Cuerpo balon;

  //          ( x0, y0,Vx0,Vy0,   m0,  R0)
  balon.Inicio(0.0,0.0,12.0,16.0,0.457,0.15);
  
  for (t=0;t<3.5;t+=dt) {
    cout<<balon.Getx() <<" " << balon.Gety() << endl;
    balon.CalculeFuerza();
    balon.Muevase(dt);
  }
  
  return 0;
}
