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
									
__global__ void g_GetValues(		float *dev_frame,	// This is the full 2D array
									int dim,			// 0 is x, 1 is y
									float *dev_x);		// This is the output array
									
__global__ void g_GetCenters(		float *dev_x,
									unsigned int *dev_broad_x);
									
#endif