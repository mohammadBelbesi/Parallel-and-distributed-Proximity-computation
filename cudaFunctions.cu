#include "cuda_runtime.h"
#include "helper_cuda.h"
#include "myProto.h"



__global__ void ComputeXYKernel(Point* points, double* x_vals, double* y_vals, int N, double t) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < N) {
        x_vals[i] = ((points[i].x2 - points[i].x1) / 2.0) * sin(t * 3.14 / 2.0) + (points[i].x2 + points[i].x1) / 2.0;
        y_vals[i] = points[i].a * x_vals[i] + points[i].b;
    }
}

void ComputeXYValues(Point* points, double* x_vals, double* y_vals, int N, double t) {
    Point* d_points = nullptr;
    double* d_x_vals = nullptr;
    double* d_y_vals = nullptr;

    cudaMalloc((void**)&d_points, N * sizeof(Point));
    cudaMalloc((void**)&d_x_vals, N * sizeof(double));
    cudaMalloc((void**)&d_y_vals, N * sizeof(double));

    cudaMemcpy(d_points, points, N * sizeof(Point), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int numBlocks = (N + threadsPerBlock - 1) / threadsPerBlock;

    ComputeXYKernel<<<numBlocks, threadsPerBlock>>>(d_points, d_x_vals, d_y_vals, N, t);

    cudaMemcpy(x_vals, d_x_vals, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(y_vals, d_y_vals, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_points);
    cudaFree(d_x_vals);
    cudaFree(d_y_vals);
}