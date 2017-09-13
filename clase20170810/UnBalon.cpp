#include <iostream>
#include <cmath>

using namespace std;

const double g = 9.81;

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
