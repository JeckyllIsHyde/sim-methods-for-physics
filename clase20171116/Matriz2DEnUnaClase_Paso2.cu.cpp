// Programa que hace c=a+b para a,b,c vectores en CUDA
#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 16
#define Ly 8
#define Nx 8
#define Ny 8
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Nx;
#define Q 5

//--------------------KERNELS----------------
__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];

__global__ void IncrementarMatriz(float *d_f0,size_t pitchf0,
				  float *d_f1,size_t pitchf1,
				  float *d_f2,size_t pitchf2,
				  float *d_f3,size_t pitchf3,
				  float *d_f4,size_t pitchf4){
  int ix,iy; float *a;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;

  a=d_f2+(ix*pitchf2)/sizeof(float)+iy;
  
  (*a)++;
}
//------------------- CLASES ----------------
class LatticeBoltzmann{
private:
  float h_w[Q];
  int h_Vx[Lx],h_Vy[Ly]; 
  float h_f0[Lx][Ly]; float*d_f0; size_t pitchf0; 
  float h_f1[Lx][Ly]; float*d_f1; size_t pitchf1; 
  float h_f2[Lx][Ly]; float*d_f2; size_t pitchf2; 
  float h_f3[Lx][Ly]; float*d_f3; size_t pitchf3; 
  float h_f4[Lx][Ly]; float*d_f4; size_t pitchf4; 
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  void Inicie(void);
  void Incremente(void);
  void Imprimase(void);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Construir las matrices virtuales en el Device
  cudaMallocPitch((void**) &d_f0,&pitchf0,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f1,&pitchf1,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f2,&pitchf2,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f3,&pitchf3,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f4,&pitchf4,Ly*sizeof(float),Lx);
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  cudaFree(d_f0);
  cudaFree(d_f1);
  cudaFree(d_f2);
  cudaFree(d_f3);
  cudaFree(d_f4);
}
void LatticeBoltzmann::Inicie(void){
  // COnstantes
  h_w[0] = 1./3; h_w[1] = h_w[2] = h_w[3] = h_w[4] = 1./6;
  h_Vx[0] = 0;
  h_Vy[0] = 0;
  // cargar vectores
  h_Vx[0] = 0; h_Vx[1] = 1; h_Vx[2] = 0; h_Vx[3] = -1; h_Vx[4] =  0;
  h_Vy[0] = 0; h_Vy[1] = 0; h_Vy[2] = 1; h_Vy[3] =  0; h_Vy[4] = -1;
  // enviarlos al device
  cudaMemcpyToSymbol(d_w, h_w, Q*sizeof(float),0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx, h_Vx, Q*sizeof(int),0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy, h_Vy, Q*sizeof(int),0, cudaMemcpyHostToDevice);
  // Funciones de distribucion
  int ix,iy;
  //Cargar valores en el Host
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      h_f0[ix][iy]=h_f1[ix][iy]=h_f2[ix][iy]=h_f3[ix][iy]=h_f4[ix][iy]=ix*Ly+iy;
  //Llevar al Device
  cudaMemcpy2D(d_f0,pitchf0,h_f0,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f1,pitchf1,h_f1,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f2,pitchf2,h_f2,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f3,pitchf3,h_f3,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  cudaMemcpy2D(d_f4,pitchf4,h_f4,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
}
void LatticeBoltzmann::Imprimase(void){
  int ix,iy;
  //Devolver los datos al Host
  cudaMemcpy2D(h_f0,Ly*sizeof(float),d_f0,pitchf0,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_f2[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;
}
void LatticeBoltzmann::Incremente(void){
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  IncrementarMatriz<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
						       d_f1,pitchf1,
						       d_f2,pitchf2,
						       d_f3,pitchf3,
						       d_f4,pitchf4);
}
// ----------------- FUNCIONES GLOBALES ----------
int main(void){
  LatticeBoltzmann Ondas;

  Ondas.Inicie();
  Ondas.Imprimase();
  Ondas.Incremente();
  Ondas.Imprimase();

  return 0;
}
