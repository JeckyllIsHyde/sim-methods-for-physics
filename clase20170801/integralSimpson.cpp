#include <iostream>
#include <cmath>

using namespace std;

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

int main(void) {

  double a, b, alfa=0;
  cout << ceroPorBiseccion(alfa, 2, 3) << endl;
  cout << ceroPorBiseccion(alfa, 5, 7) << endl;
  cout << ceroPorBiseccion(alfa, 7, 9) << endl;
  cout << ceroPorBiseccion(alfa, 11, 13) << endl;
    
  return 0;
}
