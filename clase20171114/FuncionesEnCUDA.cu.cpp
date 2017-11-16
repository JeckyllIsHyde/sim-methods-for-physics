// Suma de dos vectores c=a+b en CUDA
#include <iostream>
#include <cuda_runtime_api.h>
using namespace std;

#define Lx 128
#define Nx 64
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];

__device__ float SumeleUno(float x){
  return x+1;
}
__global__ void SumarConstante(float * d_a){
  int ix;
  ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_a[ix]=SumeleUno(d_a[ix]);
}

int main(){
  float h_w[5]; // w[i]
  int h_Vx[5],h_Vy[5]; // Vx[i],Vy[i]
  int ix;

  //DECLARAR LAS MATRICES
  //CONSTANTES
  //---Cargar las constantes en el Host-----------------
  //Cargar los pesos
  h_w[0]=1.0/3;    h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  //Cargar los vectores
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;  
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  //------Enviarlas al Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,5*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,5*sizeof(int),0,cudaMemcpyHostToDevice);

  //Datos
  float h_a[Lx];
  //En el Device
  float *d_a;  cudaMalloc((void**) &d_a,Lx*sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargar los datos en el Host
  for(ix=0;ix<Lx;ix++)
    h_a[ix]=ix;
  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumarConstante<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_a,d_a,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++)
    cout<<ix<<" "<<h_a[ix]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
