#include <iostream>
#include <cmath>

using namespace std;

const double deltat = 0.0001;
const double g=1;
const double l=1;
const double m=1;

class CuerpoRigido {
private:
  double I1,I2,I3,q0,q1,q2,q3,wx,wy,wz,Nx,Ny,Nz;
public:
  void Inicie(double Ix,double Iy,double Iz,
	      double th0,double phi0,double psi0,
	      double wx0,double wy0,double wz0);
  void Rote(double dt);
  void CalculeTorque(void);
  void CalcularXYZ(double& x,double& y,double& z);
  double getWx(void) {return wx;};
  double getWy(void) {return wy;};
  double getWz(void) {return wz;};
  double getTh(void) { return acos((q0*q0-q1*q1-q2*q2+q3*q3)/(q0*q0+q1*q1+q2*q2+q3*q3));};
};

void CuerpoRigido::Inicie(double Ix,double Iy,double Iz,
			  double th0,double phi0,double psi0,
			  double wx0,double wy0,double wz0) {
  I1 = Ix;I2 = Iy;I3 = Iz;
  wx = wx0;wy = wy0;wz = wz0;
  q0 = cos(0.5*th0)*cos(0.5*(phi0+psi0));
  q1 = sin(0.5*th0)*cos(0.5*(phi0-psi0));
  q2 = sin(0.5*th0)*sin(0.5*(phi0-psi0));
  q3 = cos(0.5*th0)*sin(0.5*(phi0+psi0));
}

void CuerpoRigido::CalculeTorque(void) {
  Nx = 2*(q2*q3+q0*q1)*m*g*l;
  Ny = -2*(q1*q3-q0*q2)*m*g*l;
  Nz = 0; 
}

void CuerpoRigido::Rote(double dt) {
  double q0old=q0,q1old=q1,q2old=q2,q3old=q3;
  double wxold=wx,wyold=wy,wzold=wz;

  q0+=dt*0.5*(-q1old*wx-q2old*wy-q3old*wz);
  q1+=dt*0.5*( q0old*wx-q3old*wy+q2old*wz);
  q2+=dt*0.5*( q3old*wx+q0old*wy-q1old*wz);
  q3+=dt*0.5*(-q2old*wx+q1old*wy+q0old*wz);

  wx+=dt*((Nx/I1)+wyold*wzold*(I2-I3)/I1);
  wy+=dt*((Ny/I2)+wxold*wzold*(I3-I1)/I2);
  wz+=dt*((Nz/I3)+wxold*wyold*(I1-I2)/I3);
}

void CuerpoRigido::CalcularXYZ(double& x,double& y,double& z) {
  double xi=x,yi=y,zi=z;

  x=(q0*q0+q1*q1-q2*q2-q3*q3)*xi+2*(q1*q2-q0*q3)*yi+2*(q1*q3-q0*q2)*zi;
  y=2*(q1*q2+q0*q3)*xi+(q0*q0-q1*q1+q2*q2-q3*q3)*yi+2*(q2*q3-q0*q1)*zi;
  z=2*(q1*q3+q0*q2)*xi+2*(q2*q3+q0*q1)*yi+(q0*q0-q1*q1-q2*q2+q3*q3)*zi;
}

int main () {
  CuerpoRigido trompo;
  double t;
  
  //           (Ix,  Iy, Iz,theta0,phi0,psi0   ,wx0,wy0,wz0);
  trompo.Inicie(0.2,0.2,1.0,M_PI/7,0.0,-M_PI/2,0,1.0,9);

  for (t=0;t<=0.5;t+=deltat) {
    cout << t << " " << trompo.getTh()*180/M_PI << endl;
    trompo.CalculeTorque();
    trompo.Rote(deltat);
  }
}
