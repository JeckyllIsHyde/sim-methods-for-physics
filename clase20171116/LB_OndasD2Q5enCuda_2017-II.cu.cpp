// Programa que hace c=a+b para a,b,c vectores en CUDA
#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 128
#define Ly 128
#define Nx 8
#define Ny 8
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Nx;

const double C=0.5;     // OJO:  C< sqrt(1/2)=0.707
const double A=10;
const double Lambda=10;
const double Omega=2*M_PI*C/Lambda;

//--------------------KERNELS----------------
__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];

__global__ void d_Adveccione(float *d_f0,size_t pitchf0,
				  float *d_f1,size_t pitchf1,
				  float *d_f2,size_t pitchf2,
				  float *d_f3,size_t pitchf3,
				  float *d_f4,size_t pitchf4,
				  float *d_f0new,size_t pitchf0new,
				  float *d_f1new,size_t pitchf1new,
				  float *d_f2new,size_t pitchf2new,
				  float *d_f3new,size_t pitchf3new,
				  float *d_f4new,size_t pitchf4new){
  //Definir variables
  int ix,iy; float *f0,*f1,*f2,*f3,*f4; float *f0new,*f1new,*f2new,*f3new,*f4new;

  //Determinar la posicion que atiende el hilo (thread)
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;

  //LEE EN fnew Y ESCRIBE EN f
  //Leer los datos
  f0new=d_f0new+(ix*pitchf0new)/sizeof(float)+iy;
  f1new=d_f1new+(ix*pitchf0new)/sizeof(float)+iy;
  f2new=d_f2new+(ix*pitchf0new)/sizeof(float)+iy;
  f3new=d_f3new+(ix*pitchf0new)/sizeof(float)+iy;
  f4new=d_f4new+(ix*pitchf0new)/sizeof(float)+iy;
  
  //Determinar las posiciones de escritura
  f0=d_f0+(((ix+d_Vx[0]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[0]+Ly)%Ly); 
  f1=d_f1+(((ix+d_Vx[1]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[1]+Ly)%Ly); 
  f2=d_f2+(((ix+d_Vx[2]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[2]+Ly)%Ly); 
  f3=d_f3+(((ix+d_Vx[3]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[3]+Ly)%Ly); 
  f4=d_f4+(((ix+d_Vx[4]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[4]+Ly)%Ly); 

  //Escribirlos en las nuevas posiciones
  (*f0)=(*f0new);
  (*f1)=(*f1new);
  (*f2)=(*f2new);
  (*f3)=(*f3new);
  (*f4)=(*f4new);
}
__device__ float d_C0(int ix,int iy){
  return 0.5;
}
__device__ float d_f0eq(float rho0,float Jx0,float Jy0,float C0){
  return rho0*(1-2*C0*C0);
}
__device__ float d_feq(float rho0,float Jx0,float Jy0,float C0,int i){
  return 3*d_w[i]*(C0*C0*rho0+d_Vx[i]*Jx0+d_Vy[i]*Jy0);
}
__device__ float d_rho(int ix,int iy,float f0,float f1,float f2,float f3, float f4, float rhoFuente){
  if(ix==Lx/2 && iy==Ly/2)
    return rhoFuente;
  else
    return f0+f1+f2+f3+f4;
}
__global__ void d_Colisione(float *d_f0,size_t pitchf0,
			  float *d_f1,size_t pitchf1,
			  float *d_f2,size_t pitchf2,
			  float *d_f3,size_t pitchf3,
			  float *d_f4,size_t pitchf4,
			  float *d_f0new,size_t pitchf0new,
			  float *d_f1new,size_t pitchf1new,
			  float *d_f2new,size_t pitchf2new,
			  float *d_f3new,size_t pitchf3new,
			  float *d_f4new,size_t pitchf4new,
			  float rhoFuente){
  //Definir variables
  int ix,iy; float *f0,*f1,*f2,*f3,*f4; float *f0new,*f1new,*f2new,*f3new,*f4new;
  float rho0,Jx0,Jy0,C0;

  //Determinar la posicion que atiende el hilo (thread)
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;

  //Leer los datos en f
  f0=d_f0+(ix*pitchf0)/sizeof(float)+iy; 
  f1=d_f1+(ix*pitchf1)/sizeof(float)+iy; 
  f2=d_f2+(ix*pitchf2)/sizeof(float)+iy; 
  f3=d_f3+(ix*pitchf3)/sizeof(float)+iy; 
  f4=d_f4+(ix*pitchf4)/sizeof(float)+iy; 

  f0new=d_f0new+(ix*pitchf0new)/sizeof(float)+iy; 
  f1new=d_f1new+(ix*pitchf1new)/sizeof(float)+iy; 
  f2new=d_f2new+(ix*pitchf2new)/sizeof(float)+iy; 
  f3new=d_f3new+(ix*pitchf3new)/sizeof(float)+iy; 
  f4new=d_f4new+(ix*pitchf4new)/sizeof(float)+iy; 
  
  //Calcular las cantidades macroscopicas
  //rho
  rho0=d_rho(ix,iy,(*f0),(*f1),(*f2),(*f3),(*f4),rhoFuente);
  //Jx0
  Jx0=(*f1)-(*f3); //Jx0=d_Vx[1]*(*f1)+d_Vx[2]*(*f2)+d_Vx[3]*(*f3)+d_Vx[4]*(*f4);
  //Jy0
  Jy0=(*f2)-(*f4); //Jy0=d_Vy[1]*(*f1)+d_Vy[2]*(*f2)+d_Vy[3]*(*f3)+d_Vy[4]*(*f4);
  //C0
  C0=d_C0(ix,iy);
  
  //Colisionar y escribe en fnew
  (*f0new)=(*f0)-2*((*f0)-d_f0eq(rho0,Jx0,Jy0,C0));
  (*f1new)=(*f1)-2*((*f1)-d_feq(rho0,Jx0,Jy0,C0,1));
  (*f2new)=(*f2)-2*((*f2)-d_feq(rho0,Jx0,Jy0,C0,2));
  (*f3new)=(*f3)-2*((*f3)-d_feq(rho0,Jx0,Jy0,C0,3));
  (*f4new)=(*f4)-2*((*f4)-d_feq(rho0,Jx0,Jy0,C0,4));
}
//------------------- CLASES ----------------
class LatticeBoltzmann{
private:
  float h_w[5]; // w[i]
  int h_Vx[5],h_Vy[5]; // Vx[i],Vy[i]

  float h_f0[Lx][Ly]; float *d_f0; size_t pitchf0;
  float h_f1[Lx][Ly]; float *d_f1; size_t pitchf1;
  float h_f2[Lx][Ly]; float *d_f2; size_t pitchf2;
  float h_f3[Lx][Ly]; float *d_f3; size_t pitchf3;
  float h_f4[Lx][Ly]; float *d_f4; size_t pitchf4;

  float h_f0new[Lx][Ly]; float *d_f0new; size_t pitchf0new;
  float h_f1new[Lx][Ly]; float *d_f1new; size_t pitchf1new;
  float h_f2new[Lx][Ly]; float *d_f2new; size_t pitchf2new;
  float h_f3new[Lx][Ly]; float *d_f3new; size_t pitchf3new;
  float h_f4new[Lx][Ly]; float *d_f4new; size_t pitchf4new;
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  void Inicie(void);
  void Adveccione(void);
  void Colisione(int t);
  void Muestre(int Who);
  double h_rho(int ix,int iy,int t);
  void Imprimase(char * NombreArchivo,int t);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Construir las matrices virtuales en el Device
  cudaMallocPitch((void**) &d_f0,&pitchf0,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f1,&pitchf1,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f2,&pitchf2,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f3,&pitchf3,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f4,&pitchf4,Ly*sizeof(float),Lx);

  cudaMallocPitch((void**) &d_f0new,&pitchf0new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f1new,&pitchf1new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f2new,&pitchf2new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f3new,&pitchf3new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f4new,&pitchf4new,Ly*sizeof(float),Lx);
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  cudaFree(d_f0);  cudaFree(d_f1);  cudaFree(d_f2);  cudaFree(d_f3);  cudaFree(d_f4);
  cudaFree(d_f0new);  cudaFree(d_f1new);  cudaFree(d_f2new);  cudaFree(d_f3new);  cudaFree(d_f4new);
}
void LatticeBoltzmann::Inicie(void){
  //CONSTANTES
  //---Cargar las constantes en el Host-----------------
  //Cargar los pesos
  h_w[0]=1.0/3;    h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  //Cargar los vectores
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;  
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  h_w[0]=1.0/3;    h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  //------Enviarlas al Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,5*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,5*sizeof(int),0,cudaMemcpyHostToDevice);
  //FUNCIONES DE DISTRIBUCION
  int ix,iy;
  //Cargar valores en el Host
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      h_f0[ix][iy]=h_f1[ix][iy]=h_f2[ix][iy]=h_f3[ix][iy]=h_f4[ix][iy]=0;
      h_f0new[ix][iy]=h_f1new[ix][iy]=h_f2new[ix][iy]=h_f3new[ix][iy]=h_f4new[ix][iy]=0;
    }
  //Llevar al Device
  cudaMemcpy2D(d_f0,pitchf0,h_f0,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f1,pitchf1,h_f1,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f2,pitchf2,h_f2,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f3,pitchf3,h_f3,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f4,pitchf4,h_f4,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);

  cudaMemcpy2D(d_f0new,pitchf0new,h_f0new,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f1new,pitchf1new,h_f1new,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f2new,pitchf2new,h_f2new,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f3new,pitchf3new,h_f3new,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f4new,pitchf4new,h_f4new,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
}
void LatticeBoltzmann::Adveccione(void){
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  d_Adveccione<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
						d_f1,pitchf1,
						d_f2,pitchf2,
						d_f3,pitchf3,
						d_f4,pitchf4,
						d_f0new,pitchf0new,
						d_f1new,pitchf1new,
						d_f2new,pitchf2new,
						d_f3new,pitchf3new,
						d_f4new,pitchf4new);
}
void LatticeBoltzmann::Colisione(int t){
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  d_Colisione<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
					       d_f1,pitchf1,
					       d_f2,pitchf2,
					       d_f3,pitchf3,
					       d_f4,pitchf4,
					       d_f0new,pitchf0new,
					       d_f1new,pitchf1new,
					       d_f2new,pitchf2new,
					       d_f3new,pitchf3new,
					       d_f4new,pitchf4new,
					       A*sin(Omega*t));
}
void LatticeBoltzmann::Muestre(int Who){
  int ix,iy;
  //Devolver los datos al Host
  cudaMemcpy2D(h_f0,Ly*sizeof(float),d_f0,pitchf0,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f1,Ly*sizeof(float),d_f1,pitchf1,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f2,Ly*sizeof(float),d_f2,pitchf2,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f3,Ly*sizeof(float),d_f3,pitchf3,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f4,Ly*sizeof(float),d_f4,pitchf4,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);

  cudaMemcpy2D(h_f0new,Ly*sizeof(float),d_f0new,pitchf0new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f1new,Ly*sizeof(float),d_f1new,pitchf1new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f2new,Ly*sizeof(float),d_f2new,pitchf2new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f3new,Ly*sizeof(float),d_f3new,pitchf3new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f4new,Ly*sizeof(float),d_f4new,pitchf4new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);

  //Mostrar Who=0 muestra a f, Who=1 muestra a fnew
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      if(Who==0) cout<<h_f0[ix][iy]<<" "; else  cout<<h_f0new[ix][iy]<<" "; 
    cout<<endl;
  }
  cout<<endl;
}
double LatticeBoltzmann::h_rho(int ix,int iy,int t){
  if(ix==Lx/2 && iy==Ly/2)
    return A*sin(Omega*t);
  else{
    double suma;
    suma=h_f0new[ix][iy]+h_f1new[ix][iy]+h_f2new[ix][iy]+h_f3new[ix][iy]+h_f4new[ix][iy];
    return suma;
  }
}
void LatticeBoltzmann::Imprimase(char * NombreArchivo,int t){
  //Devolver los datos al Host
  cudaMemcpy2D(h_f0new,Ly*sizeof(float),d_f0new,pitchf0new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f1new,Ly*sizeof(float),d_f1new,pitchf1new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f2new,Ly*sizeof(float),d_f2new,pitchf2new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f3new,Ly*sizeof(float),d_f3new,pitchf3new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy2D(h_f4new,Ly*sizeof(float),d_f4new,pitchf4new,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);

  //Imprimir en un archivo
  ofstream MiArchivo(NombreArchivo);
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++)
      MiArchivo<<ix<<" "<<iy<<" "<<h_rho(ix,iy,t)<<endl;
    MiArchivo<<endl;
  }
  MiArchivo.close();
}
// ----------------- FUNCIONES GLOBALES ----------
int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=100;

  Ondas.Inicie();

  for(t=0;t<tmax;t++){
    Ondas.Adveccione();
    Ondas.Colisione(t);
  }

  //  Ondas.Muestre(1);
  Ondas.Imprimase("Ondas_CUDA.dat",t);

  return 0;
}
