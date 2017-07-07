
#include "Device.h"



#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))

__device__ float pix_C;

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
		
//Adjust this for not only sum, but also max/min:  do a float 2 to store both the max and the min
template <unsigned int blockSize>
__global__ void g_reduce6_sum(int *g_idata, int *g_odata, unsigned int n)
{
	extern __shared__ int sdata[];
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockSize*2) + tid;
	unsigned int gridSize = blockSize*2*gridDim.x;
	sdata[tid] = 0;
	while (i < n) 
	{ 
		sdata[tid] += g_idata[i] + g_idata[i+blockSize]; 
		i += gridSize; 
	}
	__syncthreads();
	
	if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
	if (tid < 32) {
	if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
	if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
	if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
	if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
	if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
	if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
	}
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

/* //Adjust this for not only sum, but also max/min:  do a float 2 to store both the max and the min
template <unsigned int blockSize>
__global__ void g_reduce6_MinMax(float *g_idata, float *g_odata, unsigned int n)
{
	extern __shared__ float sdata1[];
	__shared__ float maxPix;
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockSize*2) + tid;
	unsigned int gridSize = blockSize*2*gridDim.x;
	sdata1[tid] = 0;
	while (i < n) 
	{ 
		maxPix = max(g_idata[i],g_idata[i + blockSize]);
		sdata1[tid] = max(sdata1[tid],maxPix); 
		i += gridSize; 
	}
	__syncthreads();
	
	if (blockSize >= 512) { if (tid < 256) { sdata1[tid] = max(sdata1[tid], sdata1[tid + 256]); } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata1[tid] = max(sdata1[tid], sdata1[tid + 128]); } __syncthreads(); }
	if (blockSize >= 128) { if (tid < 64) { sdata1[tid] = max(sdata1[tid], sdata1[tid + 64]);  } __syncthreads(); }
	if (tid < 32) {
	if (blockSize >= 64) sdata1[tid] = max(sdata1[tid], sdata1[tid + 32]);
	if (blockSize >= 32) sdata1[tid] = max(sdata1[tid], sdata1[tid + 16]);
	if (blockSize >= 16) sdata1[tid] = max(sdata1[tid], sdata1[tid + 8]);
	if (blockSize >= 8) sdata1[tid] = max(sdata1[tid], sdata1[tid + 4]);
	if (blockSize >= 4) sdata1[tid] = max(sdata1[tid], sdata1[tid + 2]);
	if (blockSize >= 2) sdata1[tid] = max(sdata1[tid], sdata1[tid + 1]);
	}
	// This shows that the final values are placed in an array of output data that needs to have a max found with it.
	if (tid == 0) g_odata[blockIdx.x] = sdata1[0];
} */


__global__ void g_getPeaks(	float *dev_frame_long,
							float *dev_gauss_Sharp,			
							float *pix_Compare,
							unsigned int *center_x,	// X-dimension center of the enlarged frame
							unsigned int *center_y, // Y-dimension center of the enlarged frame
							unsigned int real_dim_x,
							unsigned int real_dim_y,
							unsigned int iteration_num, 
							unsigned int iteration_curr,
							unsigned int grid,
							float shrinkFactor)
{
	
	// Should be 16 threads happening here.
	extern __shared__ float pixel[];// = (float*)malloc(grid * grid * sizeof(float));
	float *pix_Compare1 = (float*)malloc(sizeof(float));
	//*pix_Compare = 0;
	//float pix_Compare[1];
	//extern __shared__ float pix_Compare[];
	unsigned int idx_x = threadIdx.x % grid;
	unsigned int idx_y = threadIdx.x / grid;
	unsigned int dim_x = (real_dim_x - 20) / (1 + (1 / shrinkFactor));
	unsigned int dim_y = (real_dim_y - 20) / (1 + (1 / shrinkFactor));
	
	// Get the width of the current grid being applied
	unsigned int width_x = dim_x / (pow((double)shrinkFactor, (double)iteration_curr));
	unsigned int width_y = dim_y / (pow((double)shrinkFactor, (double)iteration_curr));
	
	// This should mark the top left spot for the square grid being used
	unsigned int start_x = (*center_x) - width_x / 2; 
	unsigned int start_y = (*center_y) - width_y / 2;	
	
	// Each thread will have a pix_x and pix_y in the frame which will be tested
	unsigned int pix_x = start_x + idx_x * width_x / (grid - 1);
	unsigned int pix_y = start_y + idx_y * width_y / (grid - 1);
	
	//__shared__ float pix_C[16] ;//= (float*)malloc(sizeof(float) * grid * grid);
	
	//pixel[MAP_2D(grid, grid, idx_x, idx_y)] = dev_frame_long[MAP_2D(real_dim_x, real_dim_y, pix_x, pix_y)];
	
	//printf("pixel[x,y]: %d \n", pixel[MAP_2D(4,4,idx_x,idx_y)]);
	//printf("dev_frame_long[x,y]: %d \n", dev_frame_long[MAP_2D(real_dim_x,real_dim_y, pix_x, pix_y)]);
	
	// This function takes the pixels to be compared, compares them with a gaussian distribution, and returns a best fit number.
	//Remember, shared memory cannot be put into a parameter
	if(threadIdx.x == 0)
	{
		printf("iteration_curr: %d \n",iteration_curr);
	}
	g_compareGauss<<<1,25, 25 * sizeof(float) * sizeof(float)>>>(	dev_frame_long,
																	dev_gauss_Sharp,
																	pix_Compare1, 
																	5,				// This is for the Gaussian grid, not the Contracting grid.
																	grid,
																	real_dim_x,
																	real_dim_y,
																	pix_x,
																	pix_y,
																	idx_x,
																	idx_y);
												
	
	
	/*
	float sum[25];
	float subtraction = 0;
	for(int i = 0; i < 5;i++)
	{
		for(int j = 0; j < 5; j++)
		{
			subtraction = dev_frame_long[MAP_2D(real_dim_x, real_dim_y, pix_x + i - 2,pix_y + j - 2)]  - dev_gauss_Sharp[MAP_2D(5,5,i,j)];
		
			// The pix_x is 1, and is causing the above to crash
			sum[idx_x + idx_y * 5] = subtraction * subtraction;
		}
	}
	
	for(int i = 1; i < 5 * 5;i++)
	{				
		sum[0] += sum[i];				
	}
	pixel[idx_x + idx_y * 4] = sum[0];
	
	*/
	
	pixel[idx_x + (idx_y * 4)] = *pix_Compare1;

	__syncthreads();

	if(threadIdx.x == 0)
	{
		int bestFit_Sharp = 0;
		
		// This needs to be replaced with the gaussian to compare against
 		for(int i = 1; i < grid * grid; i++)
		{
			if(pixel[i] == 0)
			{
				printf("iteration, i: %d %d", iteration_curr, i);
			}
			//printf("pixel[%d]: %d ", i, pixel[i]);
			//if(pixel[i] < pixel[bestFit_Sharp])
			if(pixel[i] < pixel[bestFit_Sharp])
			{
				bestFit_Sharp = i;
			}
		
		} 		
		printf("\t bestFit_Sharp: %d \n \n", bestFit_Sharp);
		//Assign bestFit_Sharp as the new center:  
		*center_x = (start_x + (bestFit_Sharp % grid) * width_x / (grid - 1));
		*center_y = (start_y + (bestFit_Sharp / grid) * width_y / (grid - 1));
		
		
		iteration_curr++;
		
		if(iteration_curr < iteration_num ){		// iteration_num
			g_getPeaks<<<1, grid * grid, grid * grid * sizeof(float)>>>(	dev_frame_long, 
																			dev_gauss_Sharp,
																			pix_Compare,
																			center_x,
																			center_y,
																			real_dim_x,
																			real_dim_y,
																			iteration_num,
																			iteration_curr,
																			grid,
																			shrinkFactor);
		}
		else
		{
			return;
		}
		
	}
	__syncthreads();
}
	
	__global__ void g_compareGauss(	float* dev_frame_long,
									float *dev_gauss_Sharp,
									float *pixel, 
									unsigned int gauss_grid,
									unsigned int grid,
									unsigned int real_dim_x,
									unsigned int real_dim_y,
									unsigned int pix_x,
									unsigned int pix_y,
									unsigned int i_x,
									unsigned int i_y)
	{
		extern __shared__ float sum[];
		unsigned int idx_x = threadIdx.x % gauss_grid;
		unsigned int idx_y = threadIdx.x / gauss_grid;
		float subtraction = 0;
		subtraction = dev_frame_long[MAP_2D(real_dim_x, real_dim_y, pix_x + idx_x - 2,pix_y + idx_y - 2)]  - dev_gauss_Sharp[MAP_2D(5,5,idx_x,idx_y)];
		
		// The pix_x is 1, and is causing the above to crash
		sum[idx_x + (idx_y * gauss_grid)] = subtraction * subtraction;
		
		// This may not work because there is more than one block
		__syncthreads();
		
		// replace this with reduce6 algorithm
		if(threadIdx.x  == 0)
		{			
			for(int i = 1;i<gauss_grid * gauss_grid;i++)
			{				
				sum[0] += sum[i];				
			}
			*pixel = sum[0];
		}
		
		
		//__syncthreads();		
	}
	
	

	
	
	
	
	