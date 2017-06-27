
#include "Device.h"



#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))

// Referred to in Device.h
__global__ void g_Zero_Array(float *dev_frame_long)
{
	int idx_x = threadIdx.x;
	int idx_y = blockIdx.x;
	
	dev_frame_long[MAP_2D(blockDim.x, gridDim.x, idx_x, idx_y)] = 0;
}

// Referred to in Device.h
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

// Referred to in Device.h
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
	
}

	// Referred to in Device.h
	__global__ void g_GetValues(	float *dev_frame,	// This is the full 2D array
									int dim,			// 0 is x, 1 is y
									float *dev_x)
	{
		// How wide is each block
		unsigned int dim_curr = blockIdx.x;
		
		// Which thread is this
		unsigned int idx = threadIdx.x;
		
		__shared__ float *array_dim;		
		array_dim = (float*)malloc(sizeof(float) * blockDim.x);
		
		__syncthreads();
		
		// sum each element in the 1D array
		if(dim == 0)
		{
		array_dim[idx] = 	dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx, dim_curr)] 		+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 1, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 2, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 3, dim_curr)]	+
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 4, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 5, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 6, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 7, dim_curr)]	+
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 8, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 9, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 10, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 11, dim_curr)]	+
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 12, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 13, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 14, dim_curr)] 	+ 
							dev_frame[MAP_2D(blockDim.x * 16, gridDim.x, 16 * idx + 15, dim_curr)];
		}
		else
		{
			// blockDim and griddim need to be switched
			array_dim[idx] = 	dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx)] 		+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 1)] 	+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 2)] 	+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 3)]	+
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 4)] 	+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 5)] 	+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 6)] 	+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 7)]	+
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 8)] 	+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 9)] 	+ 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 10)] + 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 11)]	+
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 12)] + 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 13)] + 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 14)] + 
								dev_frame[MAP_2D( gridDim.x, blockDim.x * 16, dim_curr, 16 * idx + 15)];
		}
		
		//dev_x[dim_curr] = 0;
		
		__syncthreads();
		
		if(idx == 0)
		{
			dev_x[dim_curr] = 0;
			for(int i = 0; i < blockDim.x; i++)
			{
				dev_x[dim_curr] += array_dim[i];
			}
		}		
		__syncthreads();
	}
	
	// Referred to in Device.h
	__global__ void g_GetCenters(	float *dev_x, 
									unsigned int *dev_broad_x)
	{
		unsigned int idx = threadIdx.x;
		
		__shared__ float *max;
		max = (float*)malloc(sizeof(float) * blockDim.x);
		
		//get the 21 rows with the max combined value.
		max[idx] = 0;
		
		//if(idx < 10){ printf("idx: %d", idx);}
		
		__syncthreads();
		//write a device function for this
		for(int i = 0;i<21;i++)
		{
			max[idx] = max[idx] + dev_x[idx + i];
		}
		
		__syncthreads();
		
		// Find the max value of max
		// Write a Device function for this
		
		if(idx ==0)
		{
			*dev_broad_x = 0;
			for(int i = 0; i < blockDim.x; i++)
			{
				if(max[*dev_broad_x] < max[i])
				{
					*dev_broad_x = i;
				}
			}		
		// add 11 to the max and return it to dev_broad_x
		dev_broad_x = dev_broad_x + 11;
		}
		
		__syncthreads();
			
	}
	
	// Referred to in Device.h
	__global__ void g_RemoveSharpCenter(	float *dev_frame,
											unsigned int dim_x,
											unsigned int dim_y,
											unsigned int sharp_x,
											unsigned int sharp_y)
	{
		int idx_x = threadIdx.x - 5;
		int idx_y = blockIdx.x - 5;
		
		if(idx_x == 0)
		{
			idx_x = 6;
		}
		if(idx_y == 0)
		{
			idx_y = 6;
		}
		
		// Take the frame and replace the sharpCenter with the pixels around it.  In this case, it replaces the pixels with those that are 
		//	 diagonal to it.
		dev_frame[MAP_2D(dim_x,dim_y, sharp_x + idx_x, sharp_x + idx_y)] = 
		dev_frame[MAP_2D(dim_x,dim_y, (sharp_x + idx_x) +5*((idx_x > 0) - (idx_x < 0)), (sharp_y + idx_y) + 5*((idx_y > 0) - (idx_y < 0)))];		
	}
	
	// Referred to in Device.h
	__global__ void g_Contracting_Max_Broad_Thrust(	float *dev_frame_long,
								unsigned int *max, 									
								unsigned int *center_x,
								unsigned int *center_y,
								unsigned int real_dim_x,
								unsigned int real_dim_y,
								unsigned int iteration_num, 
								unsigned int iteration_curr,
								float shrinkFactor)
{
	/*
	// Should be 16 threads happening here.
	__shared__ float pixel[16];
	float dev_array[21 * 21];
	
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
	
	// Call a global function which takes the pix_x and pix_y values, sums their surrounding values, and puts them in the pixel[] array.
	g_PutInArray<<<21,21>>>(dev_frame_long, dev_array, pix_x, pix_y);

	auto sum = thrust::reduce(thrust::seq, dev_array, dev_array + (21 * 21), (float) 0, thrust::plus<float>());

	pixel[MAP_2D(4, 4, idx_x, idx_y)] = sum;
	
	__syncthreads();
	
	//find max, recursively call the next iteration
	if(threadIdx.x == 0)
	{
		//float *maxValue;
		
		// Somehow, if this function can return the right index, this will work.
		auto maxValue = thrust::max_element(thrust::seq, pixel, pixel + 16);
		
		*max = maxValue - pixel; 
		
		printf("*max:  %d\n",*max);
		
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
	
	*/
	//__syncthreads();
}
	
__global__ void g_PutInArray(	float *dev_frame_long,
								float *dev_array,
								unsigned int center_x,
								unsigned int center_y)
{

	dev_array[threadIdx.x + blockDim.x * blockIdx.x] = dev_frame_long[MAP_2D(blockDim.x, gridDim.x, center_x + threadIdx.x - 11, center_y + blockIdx.x - 11)];
	
}
	
__global__ void g_PutInArray_1(	float *dev_frame_long,
								float *dev_array,
								unsigned int center_x,
								unsigned int center_y,
								float* storedData)
{	
	// Assign the array
	dev_array[threadIdx.x + blockDim.x * blockIdx.x] = dev_frame_long[MAP_2D(blockDim.x, gridDim.x, center_x + threadIdx.x - 11, center_y + blockIdx.x - 11)];
	
	//sync everything.
	__syncthreads();
	
	//store the data in an array to be returned
	if(threadIdx.x == 0)
	{
		for(int i = 0; i < blockDim.x; i++)
		{
			storedData[blockIdx.x] = storedData[blockIdx.x] + dev_array[blockIdx.x * blockDim.x + i];
		}
	}
	
	
}	
	
	// Referred to in Device.h
	__global__ void g_Contracting_Max_Broad(	float *dev_frame_long,
												unsigned int *max, 									
												unsigned int *center_x,
												unsigned int *center_y,
												unsigned int real_dim_x,
												unsigned int real_dim_y,
												unsigned int iteration_num, 
												unsigned int iteration_curr,
												float shrinkFactor)
{
	/*
	// Should be 16 threads happening here.
	__shared__ float pixel[16];
	float *dev_array;
	unsigned int sum = 0;
	float *dev_dataStored;	
	dev_dataStored = (float*)malloc(21 * sizeof(float));
	dev_array = (float*)malloc(21 * 21 * sizeof(float));
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
	
	// Call a global function which takes the pix_x and pix_y values, sums their surrounding values, and puts them in the pixel[] array.
	//g_PutInArray<<<21,21>>>(dev_frame_long, dev_array, pix_x, pix_y);
	
	g_PutInArray_1<<<21,21>>>(dev_frame_long, dev_array, pix_x, pix_y, dev_dataStored);
	
	for(int i = 0; i < 21; i++)
	{
		sum = sum + dev_dataStored[i];
	}
	
	
	//auto sum = thrust::reduce(thrust::seq, dev_array, dev_array + (21 * 21), (float) 0, thrust::plus<float>());

	pixel[MAP_2D(4, 4, idx_x, idx_y)] = sum;
	
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
		
		printf("*max:  %d\n",*max);
		
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
*/
}
		
	
	
	
	
	