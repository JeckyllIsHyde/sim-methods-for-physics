#include <iostream>
#include <cmath>
#include <Random64.h>

using namespace std;

const int Lx=200;
const double p=0.5;

class LatticeGas {
private:
  int v[2];    // V[i] con i=0, derecha y i=1 izquierda
  int n[Lx][2], nnew[Lx][2]; // n[ix][i]
public:
  LatticeGas(void);
  void inicie(Crandom& ran64);
  void muestre(int copia);
  void colisione(Crandom& ran64);
  void adveccione(void);
  int dondeEstaLaBolita(void);
};

LatticeGas::LatticeGas(void) {
  v[0] = 1; v[1] = -1;
}
void LatticeGas::inicie(Crandom& ran64) {
  for (int ix=0;ix<Lx;ix++)
    for (int i=0;i<2;i++)
      n[ix][i] = nnew[ix][i] = 0;
  if (ran64.r()<0.5)
    n[Lx/2][0] = 1;
  else
    n[Lx/2][1] = 1;
}
void LatticeGas::muestre(int copia) {
  for (int i=0;i<2;i++) {
    for (int ix=0;ix<Lx;ix++)
      if (copia==0)
	std::cout << n[ix][i] << " ";
      else if (copia==1)
        cout << nnew[ix][i] << " ";
    cout << endl;
  }
  std::cout << std::endl;
}
void LatticeGas::colisione(Crandom& ran64) {
  // de n a nnew
  for (int ix=0;ix<Lx;ix++) {
    if (ran64.r()<p)
      for (int i=0;i<2;i++)
	nnew[ix][i]=n[ix][i];// dejelo quieto
    else
      for (int i=0;i<2;i++)
	nnew[ix][i]=n[ix][(i+1)%2];// volteelo
  }
}
void LatticeGas::adveccione(void) {
  // de nnew a n
  for (int ix=0;ix<Lx;ix++)
    for (int i=0;i<2;i++)
      n[(ix+v[i]+Lx)%Lx][i]=nnew[ix][i];
}
int LatticeGas::dondeEstaLaBolita(void) {
  int ix=0;
  while (n[ix][0]+n[ix][1]==0)
    ix++;
  return ix;
}

// ------ funciones globales -----
const int N=1000;

double sigma2(LatticeGas* automs) {
  double xprom, xsigma2;
  int c;
  for (xprom=0,c=0;c<N;c++)
    xprom+=automs[c].dondeEstaLaBolita();
  xprom/=N;
  for (xsigma2=0,c=0;c<N;c++)
    xsigma2+=pow(automs[c].dondeEstaLaBolita()-xprom,2);
  xsigma2/=(N-1);
  return xsigma2;
}

int main() {
  LatticeGas difusion[N];
  Crandom ran64(1012);
  int c,t,tmax=400;

  for (c=0;c<N;c++)
    difusion[c].inicie(ran64);
  for (t=0;t<tmax;t++) {
    cout << t << " " << sigma2(difusion) << endl;
    for (c=0;c<N;c++)
      difusion[c].colisione(ran64);
    for (c=0;c<N;c++)
      difusion[c].adveccione();
  }
  

  return 0;
}
