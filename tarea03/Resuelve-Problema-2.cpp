#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/*  Lattice-Boltzmann para la ecuacion de onda
 *     2
 *     |0
 * 3 --*--1  con w0
 *     |
 *     4
 */
const int Lx=400;
const int Ly=200;

const int Q=5;
const double W0 = 1./3;

const double C = 0.5; // C < 0.707 celdas por click
const double _3C2 = 3*C*C;
const double AUX0 = 1-_3C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

// parametros discrimiante espejo
const double icx = 50;
const double icy = 100;
const double ipx = 100;
const double R = 100;
const double R2 = R*R; 

class LatticeBoltzmann {
private:
  double w[Q];
  int v[2][Q];    // V[alpha][i] con alpha=0 como x y alpha=1 como y
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]: probabilidad
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool calculeConLosNew);
  double Jx(int ix, int iy);
  double Jy(int ix, int iy);
  double feq(int i, double rho0, double Jx0, double Jy0);
  bool esEspejo(int ix,int iy) {
    double x=(double)ix, y=(double)iy;
    return ( ((x-icx)*(x-icx)+(y-icy)*(y-icy)>R2)
	     && (x>ipx) )? true: false;
  }
  void inicie(double rho0, double Jx0, double Jy0);
  void imponerFrontera(int ix, int iy,
		       double& rho0, double& Jx0, double& Jy0,
 		       int t);
  void colisione(int t);
  void adveccione(void);
  void imprimase(const char* nombreArchivo, int t);
};

LatticeBoltzmann::LatticeBoltzmann(void) {
  w[0] = W0;
  w[1] = w[2] = w[3] = w[4] = W0/2;

  v[0][0] = 0;
  v[1][0] = 0;

  v[0][1] = 1;  v[0][2] = 0;  v[0][3] = -1;  v[0][4] =  0;
  v[1][1] = 0;  v[1][2] = 1;  v[1][3] =  0;  v[1][4] = -1;

}
double LatticeBoltzmann::rho(int ix, int iy, bool calculeConLosNew) {
  int i; double suma;
  for (i=0, suma=0;i<Q;i++)
    if (calculeConLosNew)
      suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix, int iy) {
  int i; double suma;
  for (i=0, suma=0;i<Q;i++)
      suma+=v[0][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix, int iy) {
  int i; double suma;
  for (i=0, suma=0;i<Q;i++)
      suma+=v[1][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::feq(int i, double rho0,
			     double Jx0, double Jy0) {
  if (i==0)
    return AUX0*rho0;
  else
    return w[i]*( _3C2*rho0+3*( v[0][i]*Jx0+v[1][i]*Jy0 ) );
}
void LatticeBoltzmann::inicie(double rho0, double Jx0, double Jy0) {
  // de n a nnew
  //  double
  for (int ix=0;ix<Lx;ix++)
    for (int iy=0;iy<Ly;iy++)
      for (int i=0;i<Q;i++)
	f[ix][iy][i] = feq(i, rho0, Jx0, Jy0 );
}
void LatticeBoltzmann::imponerFrontera(int ix, int iy,
				       double& rho0, double& Jx0, double& Jy0,
				       int t) {
  double A=10, lambda=10, omega=2*M_PI*C/lambda;
  if (ix==0)
    rho0=A*sin(omega*t);
}
void LatticeBoltzmann::colisione(int t) {
  // de n a nnew
  double rho0, Jx0, Jy0; 
  for (int ix=0;ix<Lx;ix++)
    for (int iy=0;iy<Ly;iy++) {
      rho0 = rho(ix,iy, false);
      Jx0 = Jx(ix,iy), Jy0=Jy(ix,iy);
      imponerFrontera(ix, iy, rho0, Jx0, Jy0, t);
      if (esEspejo(ix,iy)) {
	fnew[ix][iy][0]=f[ix][iy][0];
	fnew[ix][iy][1]=f[ix][iy][3];
	fnew[ix][iy][2]=f[ix][iy][4];
	fnew[ix][iy][3]=f[ix][iy][1];
	fnew[ix][iy][4]=f[ix][iy][2];
      }
      else
	for (int i=0;i<Q;i++)
	  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]
	    +Utau*feq(i, rho0, Jx0, Jy0 );
    }
}
void LatticeBoltzmann::adveccione(void) {
  // de nnew a n
  for (int ix=0;ix<Lx;ix++)
    for (int iy=0;iy<Ly;iy++)
      for (int i=0;i<Q;i++)
	f[(ix+v[0][i]+Lx)%Lx][(iy+v[1][i]+Ly)%Ly][i] =
	  fnew[ix][iy][i];
}
void LatticeBoltzmann::imprimase(const char *nombreArchivo, int t) {
  ofstream MiArchivo( nombreArchivo );
  double rho0, Jx0, Jy0;
  for (int ix=0;ix<Lx/2;ix++) {
    for (int iy=0;iy<Ly;iy++) {
      rho0=rho(ix,iy,true);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
      imponerFrontera(ix,iy,rho0,Jx0,Jy0,t);
      MiArchivo << ix << " "
		<< iy  << " "
		<< rho0 << endl;
    }
    MiArchivo << endl;    
  }
  MiArchivo.close();
}

int main() {
  LatticeBoltzmann ondas;
  int t,tmax=500;

  double rho0=0, Jx0=0, Jy0=0;
  
  ondas.inicie(rho0 ,Jx0, Jy0);
  for (t=0;t<tmax;t++) {
    ondas.colisione(t);
    ondas.adveccione();
  }
  // Mostrar resultado
  ondas.imprimase("ondas_espejo.dat",t);

  return 0;
}
