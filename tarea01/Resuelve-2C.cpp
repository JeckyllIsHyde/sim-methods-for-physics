#include <iostream>
#include <cmath>

#include <vector>

using namespace std;

class vector2D {
  double _x,_y;
public:
  vector2D(double x=0.0,double y=0.0) :_x(x), _y(y) {};
  vector2D(const vector2D &v) :_x(v.x()), _y(v.y()) {};
  void show(void) const {cout << this->x() << " "
			      << this->y() << endl;};
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

const double ERR = 1e-7;

double lambda=1.0;
const double r0 = 1.0;
const double dr0 = 0.0;

vector2D f(vector2D &x, double t, double lambda) {
  double r=x.x(), dr=x.y();
  return vector2D(dr,-(t*dr+(lambda*lambda*t*t)*r)/(t*t));
}

void UnPasoDeRungeKutta4(vector2D &x, double &t,
			 double dt, double l) {
  vector2D dx1,dx2,dx3,dx4,xtmp;
  dx1 = dt*f(x,t,l);  xtmp = x+0.5*dx1;
  dx2 = dt*f(xtmp, t+0.5*dt,l);  xtmp = x+0.5*dx2;
  dx3 = dt*f(xtmp, t+0.5*dt,l);  xtmp = x+dx3;
  dx4 = dt*f(xtmp,t+dt,l);
  x+=(1.0/6)*(dx1+2*dx2+2*dx3+dx4);
}

double f_lambda(double a, double l) {
  double t=0.0,dt=0.001;
  vector2D x(r0,dr0);
  for (t=0.01;t<=a;t+=dt)
    UnPasoDeRungeKutta4(x,t,dt,l);
  return x.x();
}

double ceroPorBiseccion(double a, double b, double r) {
  double m, fa, fb, fm;
  fa = f_lambda(r,a); fb = f_lambda(r,b);
  do {
    m = (a+b)/2; fm = f_lambda(r,m);
    if (fa*fm<0)
      {b=m;fb = fm;}
    else
      {a=m;fa = fa;}
  } while (b-a > ERR);

  return (a+b)/2;
}

int main(void) {
  double a=5.0;
  int NMODOS = 25;
  double lmodos[NMODOS]; 
  double sem[] = { 1, 3, 7,10,13,
		   16,20,23,26,29,
		   33,35,38,41,44,
		   47,51,54,57,61,
		   64,67,70,74,77,80}; 
  
  for (int i=0;i<NMODOS;i++) {
    lmodos[i] = ceroPorBiseccion(sem[i],sem[i+1],a);
    //    cout << i+1 << " ";
    //    cout << lmodos[i] << " ";
    //    cout << f_lambda(a,lmodos[i]) << endl;
  }

  double t=0.0,dt=0.001;
  vector<vector2D> xs;
  for (int i=0;i<NMODOS;i++)
       xs.push_back(vector2D(r0,dr0));
  for (t=0.01;t<=a;t+=dt) {
    cout << t << " ";
    for (int j=0;j<NMODOS;j++) {
      UnPasoDeRungeKutta4(xs[j],t,dt,lmodos[j]);
      cout << xs[j].x() << " ";
    }
    cout << endl;
  }

  return 0;
}
