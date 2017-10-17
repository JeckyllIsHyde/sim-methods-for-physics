#include <iostream>
#include <cmath>

using namespace std;

const int Lx=200;
const double p=0.5;

class LatticeGas {
private:
  int v[2];    // V[i] con i=0, derecha y i=1 izquierda
  double n[Lx][2], nnew[Lx][2]; // n[ix][i]: la probabilidad de q se encuentre una bolita en una posicion
public:
  LatticeGas(void);
  void inicie(double mu, double sigma);
  void muestre(void);
  void colisione(void);
  void adveccione(void);
  double GetRho(int ix) {return n[ix][0]+n[ix][1]; };
  double GetSigma2(void);
};

LatticeGas::LatticeGas(void) {
  v[0] = 1; v[1] = -1;
}
void LatticeGas::inicie(double mu, double sigma) {
  for (int ix=0;ix<Lx;ix++)
    n[ix][0] = n[ix][1] = 1.0/(2.*sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2));
}
void LatticeGas::muestre(void) {
  for (int ix=0;ix<Lx;ix++)
    cout << GetRho(ix) << endl;
}
void LatticeGas::colisione(void) {
  // de n a nnew
  for (int ix=0;ix<Lx;ix++) {
    nnew[ix][0]=p*n[ix][0]+(1-p)*n[ix][1];// dejelo quieto
    nnew[ix][1]=p*n[ix][1]+(1-p)*n[ix][0];// dejelo quieto
  }
}
void LatticeGas::adveccione(void) {
  // de nnew a n
  for (int ix=0;ix<Lx;ix++)
    for (int i=0;i<2;i++)
      n[(ix+v[i]+Lx)%Lx][i]=nnew[ix][i];
}

// ------ funciones globales -----
const int N=1000;

double LatticeGas::GetSigma2(void) {
  double N, xprom, xsigma2;
  int ix;
  // calculo N
  for (N=0,ix=0;ix<Lx;ix++)
    N+=GetRho(ix);
  int c;
  for (xprom=0,ix=0;ix<N;ix++)
    xprom+=ix*GetRho(ix);
  xprom/=N;
  for (xsigma2=0,ix=0;ix<Lx;ix++)
    xsigma2+=pow(ix-xprom,2)*GetRho(ix);
  xsigma2/=N;
  return xsigma2;
}

int main() {
  LatticeGas difusion;
  int t,tmax=400;

  difusion.inicie(Lx/2,Lx/12);
  for (t=0;t<tmax;t++) {
    cout << t << " " << difusion.GetSigma2() << endl;
    difusion.colisione();
    difusion.adveccione();
  }
  // Mostrar resultado
  // difusion.muestre();
  return 0;
}
