#include <iostream>
#include <fstream>

#include <stdlib.h>

#include <Random64.h>

using namespace std;

const int V=240; // vertical
const int H=240; // horizontal
const int P=9; // profundidad
const int Nerf=801; // muestras de la funcion de error
const int Nr=28; // numero de receptores

class LatticeGas {
private:
  char MATRIX[V][H][P],MATRIXnew[V][H][P];
  int gaba_total;
  int bound[Nr], state[Nr], ReceptorX[Nr], ReceptorY[Nr];
  double erf[Nerf]; // Funcion de error
  double p_bound[Nr];
public:
  void lea(void);
  void gaba_inicial(void) {gaba_total=0;}
  void poner_receptores( Crandom& ran );
  void colisione(Crandom& dado,double p0, double p1);
  void avance(void) {}
};

const double tmax = 320000;

int main(void) {

  int t=0;

  Crandom num(8);

  LatticeGas sinapsis;

  sinapsis.lea();
  sinapsis.gaba_inicial();
  sinapsis.poner_receptores(num);
  //
  for ( t=0;t<=tmax;t++ ) {
    sinapsis.colisione(num,0.1,0.18);
    sinapsis.avance();
  }
  
  return 0;
}

// ****** Funciones miscelanias
char desplaza(char& numero, int rmax) {
  int r; char l, mask5, a=0, mask0;
  l=numero; mask5=1<<5; mask0=1<<0;
  if(rmax==0);
  else {
    for (r=1;r<rmax;r++)
      if((l&mask5)==mask5) {l<<=1;l=(l|mask0)&63;}
      else {l<<1;l=l&63;}
  }
  numero=1;
  return numero;
}

// ****** Implementacion de Lattice Gas
void LatticeGas::lea(void) {
  int j=0;
  double dato;

  ifstream entra("erf2.dat");
  ofstream sale("opC.dat");

  if (!entra.is_open()) {
    cerr << "No se puede abrir el archivo \"erf2.dat\"" << endl;
    exit(1);
  }
  while (!entra.eof()) {
    entra >> dato;
    erf[j] = dato;
    sale << float(j)/100-4 << "\t" << dato << "\t" << erf[j] << endl;
    j++;
  }
  entra.close();
  sale.close();
}

void LatticeGas::poner_receptores( Crandom& ran ) {
  int i,k,dx,dy;
  float hx=H/2,ky=V/2,rv=10;
  for (i=0;i<Nr;) {
    bound[i]=state[i]=0.0;
    p_bound[i]=1.0*0;
    ReceptorX[i] = int(ran.gauss(hx,rv));
    ReceptorY[i] = int(ran.gauss(ky,rv));
    if (i>0)
      for (k=0;k<i;) {
	  dx=ReceptorX[i]-ReceptorX[k];
	  dy=ReceptorY[i]-ReceptorY[k];
	if(dx!=0 || dy!=0) k++;
	else {
	  if(dx==0) ReceptorX[i] = int(ran.gauss(hx,rv));
	  if(dy==0) ReceptorY[i] = int(ran.gauss(ky,rv));
	}
      }
    i++;
  }
}

void LatticeGas::colisione(Crandom& dado,double p0, double p1) {
  int v,h,k,RMAX;
  double p; // p0: avanza, p1: rota

  for (v=1;v<V-1;v++)
    for (h=1;h<H-1;h++)
      for (k=1;k<P-1;k++)
	if ((1<=MATRIX[v][h][k])&&(MATRIX[v][h][k]<=63)) {
	  if (k==1||k==P-2/**/)
	    desplaza(MATRIX[v][h][k],3);
	  else {
	    p=dado.r();
	    if (p<=p0) RMAX=0; // no rota
	    else if (p0<p && p<=(p0+p1)) RMAX=1; // desplaza 1
	    else if ((p0+p1)<p && p<=(p0+2*p1)) RMAX=2; // desplaza 2
	    else if ((p0+2*p1)<p && p<=(p0+3*p1)) RMAX=4; // desplaza 4
	    else if ((p0+3*p1)<p && p<=(p0+4*p1)) RMAX=5; // desplaza 5
	    else if ((p0+4*p1)<p && p<=1) RMAX=3; // desplaza 3
	    desplaza(MATRIX[v][h][k],RMAX);
	  }
	}
}
