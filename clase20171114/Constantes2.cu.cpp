// Difusion 1D 
#include <iostream>
#include <cuda_runtime_api.h>
using namespace std;

//--------------------KERNELS----------------
__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];
__global__ void Test(float *d_test){
  (*d_test)=d_Vx[1];
}

int main(void){
  float h_w[5]; // w[i]
  int h_Vx[5],h_Vy[5]; // Vx[i],Vy[i]
  /*Test*/ float h_test[1]; float*d_test; cudaMalloc((void**) &d_test,sizeof(float));

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

  //TEST
  //Llevar al Device
  /*Test*/ h_test[0]=0; cudaMemcpy(d_test,h_test,sizeof(float),cudaMemcpyHostToDevice);
  
  //Cambiar en el Device
  dim3 ThreadsPerBlock(1,1,1);
  dim3 BlocksPerGrid(1,1,1);
  Test<<<BlocksPerGrid,ThreadsPerBlock>>>(d_test);

  //Devolver al Host
  cudaMemcpy(h_test,d_test,sizeof(float),cudaMemcpyDeviceToHost);

  cout<<h_test[0]<<endl;

  cudaFree(d_test);

  return 0;
}
