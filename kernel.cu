#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/async/copy.h>
#include <thrust/async/reduce.h>
#include <thrust/functional.h>
#include <thrust/complex.h>
#include <thrust/random.h>
#include <numeric>
#include <cuComplex.h>

typedef thrust::complex<double> complex;

__global__ void myKernel(   int dim,
                            int iterations,
                            int nSamples,
                            complex* Ep,
                            complex* Em,
                            complex* Pp,
                            complex* Pm,
                            complex* dEpdz,
                            complex* dEmdz,
                            complex* dEp2dz2,
                            complex* dEm2dz2,
                            complex* dEpt0,
                            complex* dEmt0,
                            complex* dEpdt,
                            complex* dEmdt,
                            complex* tmpdz2,
                            complex* tmpdz,
                            double* power
                        )
{
    double J  = 1100;
    double R1 = 0.09;
    double R2 = 1;
    double c = 299792458;
    double n = 3.3;
    double kpp = -2.000000000000000e-24;
    double aw = 400;
    double gammaK = 0;
    double g0 = 819.4299;
    double gc = 1;
    double T1 = 4.0000e-13;
    double T2 = 5.0000e-14;
    double dz =  8.0000e-07;
    double dt = 8.8061e-15;
    double Psat = 8.189848267404146e+12;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    //Initialize Ep and Em
    if(tid < dim) {
        Ep[tid] = cos(2*3.14*2*tid*tid/dim/dim);
        Em[tid] = 0;
    }

    for(int i = 0; i < iterations ; i++)
    {
        // Calculate Power
        if(tid < dim)
        {
            Pp[tid] = abs(Ep[tid]) * abs(Ep[tid]);
            Em[tid] = abs(Em[tid]) * abs(Em[tid]);
        }
        // Calculate dEpdz
        if(tid == 0)
        {
            tmpdz[tid] = Em[tid+1]*sqrt(R1);

        }
        else if(tid < dim + 1)
        {
            tmpdz[tid]  = Ep[tid-1];
        }
        __syncthreads();
        dEpdz[tid] = (tmpdz[tid+1] - tmpdz[tid])/ dz;
        //calculate dEmdz
        
    }
    
}

std::vector<double> wrapper(int dim)
{
    int threadsPerBlock = 128;
    int blocksPerGrid = 32; 
    dim = threadsPerBlock*blocksPerGrid - 2;
    complex* Ep;
    complex* Em;
    complex* Pp;
    complex* Pm;
    complex* dEpdz;
    complex* dEmdz;
    complex* dEp2dz2;
    complex* dEm2dz2;
    complex* dEpt0;
    complex* dEmt0;
    complex* dEpdt;
    complex* dEmdt;
    complex* tmpdz2;
    complex* tmpdz;
    double* power;
    size_t size = dim * sizeof(complex);
    cudaMalloc(&Ep, size);
    cudaMalloc(&Em, size);
    cudaMalloc(&Pp, size);
    cudaMalloc(&Pm, size);
    cudaMalloc(&dEpdz, size);
    cudaMalloc(&dEmdz, size);
    cudaMalloc(&dEp2dz2, size);
    cudaMalloc(&dEm2dz2, size);
    cudaMalloc(&dEpt0, size);
    cudaMalloc(&dEmt0, size);
    cudaMalloc(&dEpdt, size);
    cudaMalloc(&dEmdt, size);
    cudaMalloc(&tmpdz2, size);
    cudaMalloc(&tmpdz, size);
    cudaMalloc(&power, 5'000*sizeof(double));
    myKernel<<<blocksPerGrid, threadsPerBlock>>>(Ep, Em, Pp, Pm, dEpdz, dEmdz, dEp2dz2, dEm2dz2, dEpt0, dEmt0, dEpdt, dEmdt, tmpdz2, tmpdz, power);
    std::vector<double> power_host(5000);
    cudaError_t err = cudaMemcpy(power_host.data(), power, 5000 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(Ep);
    cudaFree(Em);
    cudaFree(Pp);
    cudaFree(Pm);
    cudaFree(dEpdz);
    cudaFree(dEmdz);
    cudaFree(dEp2dz2);
    cudaFree(dEm2dz2);
    cudaFree(dEpt0);
    cudaFree(dEmt0);
    cudaFree(dEpdt);
    cudaFree(dEmdt);
    cudaFree(tmpdz2);
    cudaFree(tmpdz);
    cudaFree(power);
    return power_host;
}