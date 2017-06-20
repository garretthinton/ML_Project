#ifndef DEVICE_H
#define DEVICE_H

#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>

__global__ void g_Initialize_Array(	float *dev_frame,
									float *dev_frame_long,
									unsigned int real_dim_x,
									unsigned int real_dim_y);

__global__ void g_Contracting_Max(	float *dev_frame_long,
									unsigned int *max, 									
									unsigned int *center_x,
									unsigned int *center_y,
									unsigned int real_dim_x,
									unsigned int real_dim_y,
									unsigned int iteration_num, 
									unsigned int iteration_curr,
									float shrinkFactor);

__global__ void g_Zero_Array(		float *dev_frame, 
									float *dev_frame_long,
									unsigned int dim_x,
									unsigned int dim_y);
									
									
#endif