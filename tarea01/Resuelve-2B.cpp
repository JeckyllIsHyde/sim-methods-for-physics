#include <iostream>
#include <cmath>

#include <vector>

using namespace std;

class vector2D {
  double _x,_y;
public:
  vector2D(double x=0.0,double y=0.0) :_x(x), _y(y) {};
  vector2D(const vector2D &v) :_x(v.x()), _y(v.y()) {};
  void show(void) {cout << _x << " " << _y << endl;};
  double x(void) const {return _x;}
  double y(void) const {return _y;}
  vector2D &operator=(const vector2D& v) {
    this->_x=v.x();    this->_y=v.y();    return *this;};
  vector2D operator+(const vector2D& v) {
    return vector2D(this->x()+v.x(),this->y()+v.y());};
  vector2D operator+=(const vector2D& v) {
    *this = *this + v;    return *this;};
  friend vector2D operator*(double a,const vector2D& v) {
    return vector2D(a*v.x(),a*v.y());};
};

double lambda=1.0;
const double r0 = 1.0;
const double dr0 = 0.0;

vector2D f(vector2D &x, double t) {
  double r=x.x(), dr=x.y();
  return vector2D(dr,-(t*dr+(lambda*lambda*t*t)*r)/(t*t));
}

void UnPasoDeRungeKutta4(vector2D &x, double &t, double dt) {
  vector2D dx1,dx2,dx3,dx4,xtmp;
  dx1 = dt*f(x,t);  xtmp = x+0.5*dx1;
  dx2 = dt*f(xtmp, t+0.5*dt);  xtmp = x+0.5*dx2;
  dx3 = dt*f(xtmp, t+0.5*dt);  xtmp = x+dx3;
  dx4 = dt*f(xtmp,t+dt);
  t+=dt;x+=(1.0/6)*(dx1+2*dx2+2*dx3+dx4);
}

double f_lambda(double a, double l) {
  double t=0.0,dt=0.01,Rf=1.0;
  vector2D x(r0,dr0);
  for (t=0.01;t<a;)
    UnPasoDeRungeKutta4(x,t,dt);
  return x.x();
}

int main(void) {
  double a=1.0;

  for (lambda=0.01;lambda<80.0;lambda+=0.01) {
    cout << lambda << " " << f_lambda(a,lambda) << endl;
  }
  
  return 0;
}
