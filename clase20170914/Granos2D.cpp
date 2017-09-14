#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double g = 9.8;
const double K = 1.e4;
const double Gamma = 50, Kcundall = 10, MU = 0.4;
const double Lx = 100.0, Ly = 100.0;

const int Nx = 1, Ny = 1;
const int N = Nx*Ny;

const double Zeta = +0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi = -0.06626458266981849;

const double ERFF = 1.e-8;

class Cuerpo;
class Colisionador;

class Cuerpo {
private:
  vector3D r,v,F,w,tau;
  double m,R,th,I;
public:
  friend class Colisionador;
  void Inicio(double  x0, double y0, double z0,
	      double Vx0, double Vy0, double Vz0,
	      double th0, double w0,
	      double m0,double R0);
  void BorreFuerzaYTorque(void);
  void AgregueFuerza(vector3D F0);
  void AgregueTorque(vector3D tau0);
  void Mueva_r(double dt, double Constante);
  void Mueva_v(double dt, double Constante);
  void Dibujese(void);
  double Getx(void) {return r.x();};
  double Gety(void) {return r.y();};
};

void Cuerpo::Inicio(double x0,double y0,double z0,
		    double Vx0,double Vy0,double Vz0,
		    double th0, double w0,
		    double m0, double R0) {
  r.cargue(x0,y0,z0);
  v.cargue(Vx0,Vy0,Vz0);
  w.cargue(0.0,0.0, w0);
  m = m0; R = R0; th = th0; I = 2.0/5*m*R*R;
}

void Cuerpo::BorreFuerzaYTorque(void) {
  F.cargue(0.0,0.0,0.0); tau.cargue(0.0,0.0,0.0);
}
				  
void Cuerpo::AgregueFuerza(vector3D F0) {
  F+=F0;
}

void Cuerpo::AgregueTorque(vector3D tau0) {
  tau+=tau0;
}

void Cuerpo::Mueva_r(double dt, double Constante) {
  r+=v*(Constante*dt); th+=w.z()*(Constante*dt);
}

void Cuerpo::Mueva_v(double dt, double Constante) {
  v+=F*(Constante*dt/m);  w+=tau*(Constante*dt/I);
}

void Cuerpo::Dibujese(void) {
  cout << ", " << r.x() << "+"<< R << "*cos(t),"
       << r.y() << "+" << R << "*sin(t), "
       << r.x() << "+" << R*cos(th)/7.0 << "*t, "
       << r.y() << "+" << R*sin(th)/7.0 << "*t ";
}

class Colisionador {
private:
  vector3D L[N+4][N+4];  bool EnColision[N+4][N+4];
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo* cuerpos, double dt);
  void CalculeLaFuerzaEntre(Cuerpo& cuerpo1, Cuerpo& cuerpo2,
			    vector3D& ele, bool& enColision,
			    double dt);
};

void Colisionador::Inicie(void) {
  int i, j;
  for(i=0; i<N; i++)
    for(j=i+1; j<N+4; j++) {
      L[i][j].cargue( 0.0,0.0,0.0 );
      EnColision[i][j] = false;
    }
}

void Colisionador::CalculeTodasLasFuerzas(Cuerpo* cuerpos, double dt) {
  int i,j;
  vector3D g_vector; g_vector.cargue(0,-g,0);
  for(i=0;i<N+4;i++)
    cuerpos[i].BorreFuerzaYTorque();
  for(i=0;i<N;i++)
    cuerpos[i].AgregueFuerza(cuerpos[i].m*g_vector);
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++)
      CalculeLaFuerzaEntre(cuerpos[i],cuerpos[j],
			   L[i][j],EnColision[i][j], dt);
}

void Colisionador::CalculeLaFuerzaEntre(Cuerpo& cuerpo1,
					Cuerpo& cuerpo2,
					vector3D& ele,
					bool& enColision,
					double dt) {
  vector3D F2, Fn, Ft, r21, n, t, vc, vcn, vct;
  double d21, s, m1, m2, m12, R1, R2, vcn_n, vct_t, Fn_n, Ft_tmax;
  r21= cuerpo2.r-cuerpo1.r;
  d21 = norma(r21);
  s = (cuerpo1.R+cuerpo2.R)-d21;
  if (s>0) {
    // geometria y dinamica del contacto
    m1 = cuerpo1.m; m2 = cuerpo2.m; m12 = m1*m2/(m1+m2);
    R1 = cuerpo1.R; R2 = cuerpo2.R;
    n = r21/d21;
    // calcular velocidad de contacto y el vector tangente
    vc = (cuerpo2.v-cuerpo1.v)-(cuerpo2.w^n)*R2-(cuerpo1.w^n)*R1;
    vcn_n = vc*n;    vcn = n*vcn_n;
    vct = vc-vcn;    vct_t = norma(vct);
    if (vct_t < ERFF)
      t.cargue(0,0,0);
    else
      t = vct/vct_t;

    // fuerzas normales
    // fuerza de Hertz
    Fn_n = K*pow(s,1.5);
    // disipacion plastica
    Fn_n-= m12*sqrt(s)*Gamma*vcn_n;
    if (Fn_n<0)
      Fn_n = 0;
    Fn = n*Fn_n;

    // fuerzas tangenciales
    // fuerza estatica
    ele+=(vct*dt);
    Ft = ele*(-Kcundall);
    // fuerza cinetica
    Ft_tmax = MU*Fn_n;
    if ( norma(Ft)>Ft_tmax )
      Ft = ele*(-Ft_tmax/norma(ele));

    // construir fuerza total
    F2 = Fn+Ft;
    cuerpo2.AgregueFuerza(F2);           cuerpo2.AgregueTorque((n*(-R2))^Ft);
    cuerpo1.AgregueFuerza(F2*(-1.0));    cuerpo1.AgregueTorque((n*R1)^(Ft*(-1.0)));

    enColision = true;
  }
  else if (enColision == true) {
    ele.cargue(0,0,0); enColision=false;
  }
}

void InicieAnimacion(void) {
  // cout << "set terminal gif animate" << endl; 
  // cout << "set output 'Granos.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [-10:110]" << endl;
  cout << "set yrange [-10:110]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set trange [0:7]" << endl;
  cout << "set isosamples 12" << endl;
}

void InicioCuadro(void) {
  cout << "set title 'Pariculas'\n plot 0,0 ";
  cout << " , " << 100./7 << "*t, 0";
  cout << " , " << 100./7 << "*t, 100";
  cout << " , 0," << 100./7 << "*t";
  cout << " , 100," << 100./7 << "*t";
}

void TermineCuadro(void) {
  cout << endl;
}

int main(void) {
  double t, dt=1e-3;
  double tdibujo; int Ndibujos;
  Cuerpo granos[N+4];
  int i,j;
  Colisionador newton; //newton.Inicie();
  Crandom ran64(1); double theta;

  double m0=1, R0=6, v=10, w = 10;
  double Rpared=10000, Mpared=1000*m0 ;

  double T=Lx/v, tmax=2*T;

  InicieAnimacion(); Ndibujos=2000;
  //                (  x0,       y0, z0,Vx0,Vy0,Vz0,th0,W0,     m0,     R0 )
  // pared arriba
  granos[N+0].Inicio(Lx/2,Ly+Rpared,0.0,0.0,0.0,0.0,  0, 0, Mpared, Rpared );
  // pared abajo
  granos[N+1].Inicio(Lx/2,  -Rpared,0.0,0.0,0.0,0.0,  0, 0, Mpared, Rpared );
  // pared derercha
  granos[N+2].Inicio(Lx+Rpared,Ly/2,0.0,0.0,0.0,0.0,  0, 0,  Mpared, Rpared );
  // pared izquierda
  granos[N+3].Inicio(  -Rpared,Ly/2,0.0,0.0,0.0,0.0,  0, 0,  Mpared, Rpared );

  double  row, col;
  for (i=0; i<N; i++) {
    row = (int)i/(int)Nx; col = i-row*Nx;
    theta = 2*M_PI*ran64.r();
    granos[i].Inicio( (1+col)*Lx/(Nx+2), (1+row)*Ly/(Ny+2), 0.0,
		      0*v*cos(theta),0*v*sin(theta),0.0, 0.0, w,
		      m0, R0 );
  }

  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt) {
    if (tdibujo>tmax/Ndibujos) {
      InicioCuadro();
      for (i=0; i<N; i++)
	granos[i].Dibujese();
      TermineCuadro();
      tdibujo = 0;
    }
    // cout << granos[1].Getx() <<" " << granos[1].Gety() << endl; 
    // Muevase con Omelyan PEFRL
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Zeta);
    newton.CalculeTodasLasFuerzas(granos, dt); for (i=0;i<N;i++) granos[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(granos, dt); for (i=0;i<N;i++) granos[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,1-2*(Xi+Zeta));
    newton.CalculeTodasLasFuerzas(granos, dt); for (i=0;i<N;i++) granos[i].Mueva_v(dt,Lambda);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Xi);
    newton.CalculeTodasLasFuerzas(granos, dt); for (i=0;i<N;i++) granos[i].Mueva_v(dt,(1-2*Lambda)/2);
    for (i=0;i<N;i++) granos[i].Mueva_r(dt,Zeta);
  }

  return 0;
}
