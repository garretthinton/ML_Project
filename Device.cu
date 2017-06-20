
#include "Device.h"
#include <stdio.h>
#include <math.h>

#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))


/*
Title:  d_Initialize_Array
Description:  This kernal function takes the frame that is given and expands it so that it can be used as part of the Contracting Grid algorithm.  
Inputs:  
	-float *dev_frame:
	-unsigned int real_dim_x: 
	-unsigned int real_dim_y:
Outputs: 
	-float *dev_frame_long:
*/
__global__ void g_Zero_Array(float *dev_frame, float *dev_frame_long, unsigned int dim_x, unsigned int dim_y)
{
	int idx_x = threadIdx.x;
	int idx_y = blockIdx.x;
	//int real_dim_x = blockDim.x;
	//int real_dim_y = gridDim.x;
	
	dev_frame_long[MAP_2D(blockDim.x, gridDim.x, idx_x, idx_y)] = 0;
}

__global__ void g_Initialize_Array(float *dev_frame, float *dev_frame_long, unsigned int real_dim_x, unsigned int real_dim_y)
{ 
		int idx_x = threadIdx.x;
		int idx_y = blockIdx.x;
		int dim_x = blockDim.x;
		int dim_y = gridDim.x;
		int extend_x = (real_dim_x - dim_x)/2;
		int extend_y = (real_dim_y - dim_y)/2;
		
		dev_frame_long[MAP_2D(real_dim_x, real_dim_y, idx_x + extend_x, idx_y + extend_y)] = dev_frame[MAP_2D(dim_x,dim_y,idx_x, idx_y)];
}

/*
Title:  d_Contracting_Max
Description:  This kernal function takes the big frame that is made in the d_Initialize_Array() function and chooses 
Inputs:  
	-float *dev_frame:
	-unsigned int real_dim_x: 
	-unsigned int real_dim_y:
Outputs: 
	-float *dev_frame_long:
*/
__global__ void g_Contracting_Max(	float *dev_frame_long,
									unsigned int *max, 									
									unsigned int *center_x,	// X-dimension center of the enlarged frame
									unsigned int *center_y, // Y-dimension center of the enlarged frame
									unsigned int real_dim_x,
									unsigned int real_dim_y,
									unsigned int iteration_num, 
									unsigned int iteration_curr,
									float shrinkFactor)
{	
	// Should be 16 threads happening here.
	__shared__ float pixel[16];

	unsigned int idx_x = threadIdx.x % 4;
	unsigned int idx_y = threadIdx.x / 4;
	unsigned int dim_x = real_dim_x / (1 + (1 / shrinkFactor));
	unsigned int dim_y = real_dim_y / (1 + (1 / shrinkFactor));
	
	// Get the width of the current grid being applied
	unsigned int width_x = dim_x / (pow((double)shrinkFactor, (double)iteration_curr));
	unsigned int width_y = dim_y / (pow((double)shrinkFactor, (double)iteration_curr));
	
	// This should mark the top left spot for the square grid being used
	unsigned int start_x = (*center_x) - width_x / 2; 
	unsigned int start_y = (*center_y) - width_y / 2;	
	
	// Each thread will have a pix_x and pix_y in the frame which will be tested
	unsigned int pix_x = start_x + idx_x * width_x / 3;
	unsigned int pix_y = start_y + idx_y * width_y / 3;
	
	pixel[MAP_2D(4, 4, idx_x, idx_y)] = dev_frame_long[MAP_2D(real_dim_x, real_dim_y, pix_x, pix_y)];
	//printf("pixel[x,y]: %d \n", pixel[MAP_2D(4,4,idx_x,idx_y)]);
	//printf("dev_frame_long[x,y]: %d \n", dev_frame_long[MAP_2D(real_dim_x,real_dim_y, pix_x, pix_y)]);
	
	
	__syncthreads();
	
	//find max, recursively call the next iteration
	if(threadIdx.x == 0)
	{
		*max = 0;
		for(int i = 1; i < 16; i++)
		{
			if(pixel[i] > pixel[*max])
			{
				*max = i;
			}
		}
		
		//Assign max as the new center:  
		*center_x = (start_x + (*max % 4) * width_x / 3);
		*center_y = (start_y + (*max / 4) * width_y / 3);
		
		
		iteration_curr++;
		if(iteration_curr < iteration_num){			
			g_Contracting_Max<<<1,16>>>(dev_frame_long, max, center_x, center_y, real_dim_x, real_dim_y, iteration_num, iteration_curr, shrinkFactor);
		}
		else
		{
			return;
		}
	}
	
	//__syncthreads();
	
	}
