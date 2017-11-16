#include <iostream>
#include <fstream>
#include <cmath>

#include <GL/glew.h>
#include <GL/glut.h>

#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

using namespace std;

/* PROGRAMAR EN CUDA:
   Se trabaja con dos copias de datos
     1. en el Host
     2. en el Divice
   Se realiza 6 pasos
     1. Declarar datos en el Host y el Device
     2. Inicializar y copiar datos en el Host y el Device
     3. Organizar Threads, Blocks and Grids
     4. Devolver datos al Host
     5. Post-Procesar en el Host
     6. Liberar memoria dinámica
*/

// CUDA extiende el lenguaje de C/C++ con directivas
__constant__ float d_w[5]; // en memoria de constates
__device float SumeleUno(float x ){
  return x+1;
}
__global__ void SumeleAlgo(float *d_test) {
  ( *d_test )+=SumeleUno( d_w[0] );
}

int main(void ) {
  // Declarar los datos MATRICES
  float h_test[1]; // memoria en el host
  float h_w[5]; // mermoria para constantes
  h_w[0] = 1.0/3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6
  float* d_test; // debe ser un apuntador para el device
  cudaMalloc((void**) &d_test, sizeof(float)); // CUDA asigna

  // Inicializar los datos
  //  cargan en el host
  h_test[0] = 10;
  //  enviar al device (a donde, de donde, tamaño, dirección)
  cudaMemcpy(d_test,h_test,sizeof(float),cudaMemcpyHostToDevice);
  cudaMemToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);

  // Organizar y ejecutar el proceso en el device
  // Dividir datos por blocks y grids, para especificar tareas
  // Los blocks están formados por hilos (threads)
  dim3 ThreadsPerBlock(1,1,1);
  dim3 BlockPerGrid(1,1,1);
  SumeleAlgo<<<BlocksPerGrid,ThreadsPerBlock>>>(d_test);

  // Devolver los datos al host
  cudaMemcpy(h_test,d_test,sizeof(float),cudaMemcpyDeviceToHost);

  // Procesar en el host
  cout << h_test[0] << endl;

  // Liberar memoria dinámica
  cudaFree(d_test);

  return 0;
}
