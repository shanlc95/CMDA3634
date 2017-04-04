/*

To compile:
nvcc -O3  -o cudaMandelbrot cudaMandelbrot.cu png_util.c -I. -lm -lpng

To run:
./cudaMandelbrot

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cuda.h>

extern "C"
{
#include "png_util.h"
}

__global__ void mandelbrot(const int NRe, 
			   const int NIm,
			   const float minRe,
			   const float minIm,
			   const float dRe,
			   const float dIm,
			   const float cRe,
			   const float cIm,
			   float* __restrict__ c_count ) {

  int n,m;
  
  n = threadIdx.x + blockIdx.x*blockDim.x;
  m = threadIdx.y + blockIdx.y*blockDim.y;

  float zRe = minRe + n*dRe;
  float zIm = minIm + m*dIm;
 
  int Nt = 200;
  int t, cnt=0;
  for(t=0;t<Nt;++t){
    
    // z = z^2 + c
    //   = (zRe + i*zIm)*(zRe + i*zIm) + (cRe + i*cIm)
    //   = zRe^2 - zIm^2 + 2*i*zIm*zRe + cRe + i*cIm
    float zReTmp = zRe*zRe - zIm*zIm + zRe;
    zIm = 2.f*zIm*zRe + cIm;
    zRe = zReTmp;
    
    cnt += (zRe*zRe+zIm*zIm<4.f);
  }
  
  c_count[n + m*NRe] = cnt;

}


int main(int argc, char **argv){

  const int NRe = 2048;
  const int NIm = 2048;

  const float cRe = 0.285;
  const float cIm = 0.01;

  /* box containing sample points */
  const float centRe = -.759856, centIm= .125547;
  const float diam  = 0.151579;
  const float minRe = centRe-0.5*diam;
  const float remax = centRe+0.5*diam;
  const float minIm = centIm-0.5*diam;
  const float immax = centIm+0.5*diam;

  const float dRe = (remax-minRe)/(NRe-1.f);
  const float dIm = (immax-minIm)/(NIm-1.f);

  float *h_count = (float*) calloc(NRe*NIm, sizeof(float));

  float *c_count;
  
  cudaMalloc(&c_count, NRe*NIm*sizeof(float));
  
  // specify two dimensional array of threads
  int TRe = 16;
  int TIm = 16;
  int BRe = NRe/TRe;
  int BIm = NIm/TIm;
  
  dim3 B(BRe, BIm);
  dim3 T(TRe, TIm);
  
  double tic = clock();


  // call mandelbrot from here
  mandelbrot <<< B, T >>> (NRe, NIm, minRe, minIm, dRe, dIm, cRe, cIm, c_count);

  cudaDeviceSynchronize();

  double toc = clock();

  double elapsed = (toc-tic)/CLOCKS_PER_SEC;
  
  printf("elapsed time %f\n", elapsed);

  cudaMemcpy(h_count, c_count, NRe*NIm*sizeof(float), cudaMemcpyDeviceToHost);

  FILE *png = fopen("cudaMandelbrot.png", "w");
  write_hot_png(png, NRe, NIm, h_count, 0, 80);
  fclose(png);

}
