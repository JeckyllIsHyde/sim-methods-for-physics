/*

  integralSimson.cpp : 
  Is an example of how to find zeros in bessel functions, the 
  bessel function is calculated with Simpson's Rule. Zeros are
  found by bisections.

 */

#include <iostream>
#include <cmath>

using namespace std;

const double ERR = 1e-7;

double f(double alfa, double x, double t) {
  return cos(alfa*t-x*sin(t));
}

double IntegralPorSimpson(double alfa, double x, double a, double b, double N) {
  double h = (b-a)/N, suma,t;
  int i;
  for ( suma=0, i=0; i<=N; i++ ) {
    t = a+i*h;
    if ( i==0 || i==N )
      suma += f(alfa,x,t);
    else if ( i%2==0 )
      suma += 2*f(alfa,x,t);
    else
      suma += 4*f(alfa,x,t);    
  }
  return suma*h/3;
}

double Bessel(double alfa, double x) {
  return IntegralPorSimpson(alfa,x,0,M_PI,50)/M_PI;
}

double ceroPorBiseccion(double alfa, double a, double b) {
  double m, fa, fb, fm;
  fa = Bessel(alfa,a); fb = Bessel(alfa,b);
  do {
    m = (a+b)/2; fm = Bessel(alfa,m);
    if (fa*fm<0)
      {b=m;fb = fm;}
    else
      {a=m;fa = fa;}
  } while (b-a > ERR);

  return (a+b)/2;
}

int main(void) {

  double a, b, alfa=0;
  cout << ceroPorBiseccion(alfa, 2, 3) << endl;
  cout << ceroPorBiseccion(alfa, 5, 7) << endl;
  cout << ceroPorBiseccion(alfa, 7, 9) << endl;
  cout << ceroPorBiseccion(alfa, 11, 13) << endl;
    
  return 0;
}
