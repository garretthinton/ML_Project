//	University of Arizona
//	Center for Gamma-Ray Imaging
//	
//
//
//
//
//
//
//  
// Author and Revision Date 
// Garrett Hinton		Jun 13 2017
//
//
//
// Title:	****Contracting_Grid****
// Description:  This class takes in a frame, as well as some parameters for how the user would like the algorithm to behave.  // It has the capabilities to output the center points for the first and second screen in the imaging process.  It also has the 
// ability to output the direction of propagation for beta particle being imaged.  The idea of this class is that using the
// Cine class, the user could read in the data and analyze it in real time.
		
#include "Contracting_Grid.h"

using namespace std;
		
#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))

#define CHECK(x) do {\
	cudaError_t err = (x);\
	if ( err!= cudaSuccess)\
	{\
		fprintf(stderr, "API error %s:%d Returned:%d\n", \
		__FILE__, __LINE__, err);\
		printf(cudaGetErrorString(err));\
		exit(1);\
	}\
	}while(0)
		
	//#define min(a, b) ((a) > (b))? (b): (a)
	//#define max(a, b) ((a) > (b))? (a): (b)
	
	// could be the max function
//Adjust this for not only sum, but also max/min:  do a float 2 to store both the max and the min
template <unsigned int blockSize>
__global__ void g_reduce6_MinMax(float *g_idata, float *g_odata, unsigned int n)
{
	// Initialize these variables
	extern __shared__ float sdata1[];
	__shared__ float maxPix[512];
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*(blockSize*2) + tid;
	unsigned int gridSize = blockSize*2*gridDim.x;
	sdata1[tid] = 0;
	
	while (i < (n - blockSize)) 
	{ 
		maxPix[tid] = max(g_idata[i],g_idata[i + blockSize]);
		sdata1[tid] = max(sdata1[tid],maxPix[tid]); 
		i += gridSize; 
	}
	__syncthreads();
	
	if (blockSize >= 512) { if (tid < 256) { sdata1[tid] = max(sdata1[tid], sdata1[tid + 256]); } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata1[tid] = max(sdata1[tid], sdata1[tid + 128]); } __syncthreads(); }
	if (blockSize >= 128) { if (tid < 64) { sdata1[tid] = max(sdata1[tid], sdata1[tid + 64]);  } __syncthreads(); }
	// These are all on the same warp, and therefore do not need a __syncthreads() command
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
	
}
	
// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid()
{
	// standard assumptions for algorithms
	iterations = 7;
	shrinkFactor = 2;
	distance = .075;
	grid = 4;
	sigma_Sharp = 3;
	sigma_Broad = 10;
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y  = real_dim_x = real_dim_y = dim_x = dim_y = 0;
	
	// May need to instantiate 'frame'
	frame = NULL;
	gauss_Sharp = NULL;
	gauss_Broad = NULL;
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c)
{
	// standard assumptions for algorithms
	shrinkFactor = 2;
	distance = .075;
	grid = 4;
	sigma_Sharp = 5;
	sigma_Broad = 20;
	
	// Initialize the sharp gaussian array
	gauss_Sharp = new float[25];
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{ 
			gauss_Sharp[MAP_2D(5,5,i,j)] = exp((-1.0/2.0 * (pow((2.0 - (double)i ),2.0) + pow((2.0 - (double)j),2.0))) / pow(sigma_Sharp,2.0));
		}
	}
	
	gauss_Broad = new float[21*21];
	for(int i = 0; i < 21; i++)
	{
		for(int j = 0; j < 21; j++)
		{ 
			gauss_Broad[MAP_2D(21,21,i,j)] = exp((-1.0/2.0 * (pow((10.0 - (double)i ),2.0) + pow((10.0 - (double)j),2.0))) / pow(sigma_Broad,2.0));
		}
	}
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	// Get the biggest dimension of the frame and put it into dim_Max
	unsigned int dim_Max = dim_x;
	
	if(dim_y > dim_x)
	{
		dim_Max = dim_y;
	}
	
	// Using dim_Max, calculate the iterations needed for the algorithm to succeed
	iterations = log(dim_Max)/log(shrinkFactor);
	
	// The +20 is to account for the extra 10 spaces that the broad algorithm will go outside of its chosen point.
	real_dim_x = (dim_x + (dim_x / shrinkFactor)) + 20;
	real_dim_y = (dim_y + (dim_y / shrinkFactor)) + 20;
	
	//
	CHECK(cudaMalloc((void**)&dev_frame_long, real_dim_x * real_dim_y * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_frame, dim_x * dim_y * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_a, 128 * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_a_final, sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_center_x, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_y, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&pix_Compare, grid * grid * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_gauss_Sharp, 25 * sizeof(float)));
	
	// Copy the data to the GPU
	//CHECK(cudaMemcpy(dev_frame, frame, dim_x * dim_y * sizeof(float), cudaMemcpyHostToDevice));
	center_x = real_dim_x/2;
	center_y = real_dim_y/2;
	iteration_curr = 0;
	A = 0;
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c, float shrinkFactor_in)
{	
	// Checks that the value for the shrinkFactor is acceptable
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}
	else
	{
		cout << "shrinkFactor value must be a floating point number bigger than 1. " << endl;
		exit(EXIT_FAILURE);
	}
	// standard assumption for algorithms
	distance = .075;
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	// Get the biggest dimension of the frame and put it into dim_Max
	unsigned int dim_Max = dim_x;
	
	if(dim_y > dim_x)
	{
		dim_Max = dim_y;
	}
	
	// Using dim_Max, calculate the iterations needed for the algorithm to succeed
	iterations = log(dim_Max)/log(shrinkFactor);
	
	real_dim_x = dim_x + (dim_x / shrinkFactor);
	real_dim_y = dim_y + (dim_y / shrinkFactor);
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c, unsigned int iterations_in, float shrinkFactor_in)
{
	// standard assumptions for algorithms			
	distance = .075;
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	//Checks that the value for iterations is acceptable
	if(iterations_in > 0)
	{
		iterations = iterations_in;
	}
	else
	{
		cout << "Iterations value must be an integer bigger than 0. " << endl;
		exit(EXIT_FAILURE);
	}
	
	//Checks that the value for the shrinkFactor is acceptable
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}
	else
	{
		cout << "shrinkFactor value must be a floating point number bigger than 1. " << endl;
		exit(EXIT_FAILURE);
	}
	
	real_dim_x = dim_x + (dim_x / shrinkFactor);
	real_dim_y = dim_y + (dim_y / shrinkFactor);
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c, float shrinkFactor_in, float distance_in)
{
	shrinkFactor = 2;
	distance = .075;
	
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}
	else
	{
		cout << "shrinkFactor value must be a floating point number bigger than 1. " << endl;
		exit(EXIT_FAILURE);
	}
	
	
	if(distance_in > 0)
	{
		distance = distance_in;
		
	}
	else
	{
		cout << "Distance value must be an integer bigger than 0. " << endl;
		exit(EXIT_FAILURE);
	}
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	// Get the biggest dimension of the frame and put it into dim_Max
	unsigned int dim_Max = dim_x;
	
	if(dim_y > dim_x)
	{
		dim_Max = dim_y;
	}
	
	// Using dim_Max, calculate the iterations needed for the algorithm to succeed
	iterations = log(dim_Max)/log(shrinkFactor);
	
	real_dim_x = dim_x + (dim_x / shrinkFactor);
	real_dim_y = dim_y + (dim_y / shrinkFactor);
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c, unsigned int iterations_in, float shrinkFactor_in, float distance_in)
{
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	//Checks that the value for iterations is acceptable
	if(iterations_in > 0)
	{
		iterations = iterations_in;
	}
	else
	{
		cout << "Iterations value must be an integer bigger than 0. " << endl;
		exit(EXIT_FAILURE);
	}
	
	//Checks that the value for the shrinkFactor is acceptable
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}
	else
	{
		cout << "shrinkFactor value must be a floating point number bigger than 1. " << endl;
		exit(EXIT_FAILURE);
	}

	//Checks that the value for the distance is acceptable
	if(distance_in > 0)
	{
		distance = distance_in;
	}
	else
	{
		cout << "Distance value must be an integer bigger than 0. " << endl;
		exit(EXIT_FAILURE);
	}
	
	real_dim_x = dim_x + (dim_x / shrinkFactor);
	real_dim_y = dim_y + (dim_y / shrinkFactor);
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(unsigned int iterations_in, float shrinkFactor_in, float distance_in)
{
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y = real_dim_x = real_dim_y = 0;
	
	// The real dimensions of the frame will be received later if this is the contructor used
	dim_x = dim_y = 0;
	
	frame = NULL;
	//Checks that the value for iterations is acceptable
	if(iterations_in > 0)
	{
		iterations = iterations_in;
	}
	else
	{
		cout << "Iterations value must be an integer bigger than 0. " << endl;
		exit(EXIT_FAILURE);
	}
	
	//Checks that the value for the shrinkFactor is acceptable
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}
	else
	{
		cout << "shrinkFactor value must be a floating point number bigger than 1. " << endl;
		exit(EXIT_FAILURE);
	}

	//Checks that the value for the distance is acceptable
	if(distance_in > 0)
	{
		distance = distance_in;
	}
	else
	{
		cout << "Distance value must be an integer bigger than 0. " << endl;
		exit(EXIT_FAILURE);
	}
}

//Commented
// Essence of the contracting grid algorithm and the purpose of this class.  Referred to in Contracting_Grid.h
void Contracting_Grid::findSharpCenter_Rec()
{
	/*
	float *dev_frame_long;
	
	// This is the frame data that will be going into the GPU
	float *dev_frame;
	unsigned int *max;
	unsigned int *dev_center_x;
	unsigned int *dev_center_y;
	unsigned int iteration_curr = 0;
	unsigned int center_x = real_dim_x / 2;
	unsigned int center_y = real_dim_y / 2;
	unsigned int totalSamples = real_dim_x * real_dim_y;
	unsigned int dim_Samples = dim_x * dim_y;
	
	// Allocate memory in the GPU
	CHECK(cudaMalloc((void**)&dev_frame, dim_Samples * sizeof(float)));
	CHECK(cudaMalloc((float**)&dev_frame_long, totalSamples * sizeof(float)));
	
	CHECK(cudaMalloc((void**)&max, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_x, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_y, sizeof(unsigned int)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, dim_Samples * sizeof(float), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(dev_center_x, &center_x, sizeof(unsigned int), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(dev_center_y, &center_y, sizeof(unsigned int), cudaMemcpyHostToDevice));
	
	g_Zero_Array<<<real_dim_x, real_dim_y>>>(	dev_frame_long);
	
	g_Initialize_Array<<<dim_x, dim_y>>>(	dev_frame,
											dev_frame_long,
											real_dim_x,
											real_dim_y);
	
	// may need a cudaDeviceSynchronize here
	cudaDeviceSynchronize();
	CHECK(cudaGetLastError());
	
	g_Contracting_Max<<<1,16>>>(dev_frame_long,
								max, 									
								dev_center_x,
								dev_center_y,
								real_dim_x,
								real_dim_y,
								iterations, 
								iteration_curr,
								shrinkFactor);
	
	CHECK(cudaGetLastError());
	// may need a cudaDeviceSynchronize here
	cudaDeviceSynchronize();
	
	// Copy the data back to the host
	CHECK(cudaMemcpy(&center_x, (void*)dev_center_x, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(&center_y, (void*)dev_center_y, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	
	sharp_x = center_x - ((real_dim_x - dim_x) / 2);
	sharp_y = center_y - ((real_dim_y - dim_y) / 2);
	
	float temp1 = frame[MAP_2D(dim_x,dim_y,0,0)];
	int tempi = 0;
	int tempj = 0;
	for(int i = 0; i< dim_x;i++){
		for(int j = 0; j < dim_y; j++){
			if(temp1 < frame[MAP_2D(dim_x,dim_y,i,j)]){
				temp1 = frame[MAP_2D(dim_x,dim_y,i,j)];
				tempi = i;
				tempj = j;
			}
		}
	}

	
	cout<<"True max: "<< temp1 <<endl;
	cout<< "i , j: \t" << tempi << " , " << tempj <<endl;
	
	
	// Free the data
	CHECK(cudaFree(dev_frame));
	CHECK(cudaFree(dev_frame_long));
	CHECK(cudaFree(max));
	CHECK(cudaFree(dev_center_x));
	CHECK(cudaFree(dev_center_y));	
	*/
}

//Commented
// referred to in Contracting_Grid.h
void Contracting_Grid::findBroadCenter_brute()
{	
	/*
	// This is the frame data that will be going into the GPU
	float *dev_frame;
	float *dev_x;
	float *dev_y;
	unsigned int *dev_broad_x;
	unsigned int *dev_broad_y;
	unsigned int dim_Samples = dim_x * dim_y;
	
	// Allocate memory in the GPU
	CHECK(cudaMalloc((void**)&dev_frame, dim_Samples * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_x, dim_x * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_y, dim_y * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_broad_x, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_broad_y, sizeof(unsigned int)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, dim_Samples * sizeof(float), cudaMemcpyHostToDevice));
	
	// Take the sharpCenter(x,y) and delete those rows and columns + or - 10 lines
	g_RemoveSharpCenter<<<11,11>>>(dev_frame, dim_x, dim_y, sharp_x, sharp_y);
	
	
	//0 or 1 depending on x or y dimension collection
	if(dim_x % 16 == 0 && dim_y % 16 == 0)
	{
		g_GetValues<<<dim_x, dim_y/16>>>(dev_frame, 0, dev_x);
		g_GetValues<<<dim_y, dim_x/16>>>(dev_frame, 1, dev_y);
		
		g_GetCenters<<<1, dim_x - 21>>>(dev_x, dev_broad_x);
		g_GetCenters<<<1, dim_y - 21>>>(dev_y, dev_broad_y);
	}
	
	// Add up the sums for every 21 lines of data up to the sharpCenter +- 10.  If the lines are increasing, 
	//   interpolate the data to adjust the data for the lines that were ignored.  Choose the center of the highest 21 lines 
	CHECK(cudaMemcpy(&broad_x, (void*)dev_broad_x, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(&broad_y, (void*)dev_broad_y, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	
	cout<< "broad_x: " << broad_x << endl;
	cout<< "broad_y: " << broad_y << endl;
	
	cudaFree(dev_frame);
	cudaFree(dev_x);
	cudaFree(dev_y);
	cudaFree(dev_broad_x);
	cudaFree(dev_broad_y);
	*/
}

// referred to in Contracting_Grid.h
/*void Contracting_Grid::findBroadCenter_Rec_Thrust()
{
	// This is the frame data that will be going into the GPU
	float *dev_frame;
	float *dev_frame_long;
	unsigned int *max;
	unsigned int *dev_center_x;
	unsigned int *dev_center_y;
	unsigned int iteration_curr = 0;
	unsigned int real_dim_x_broad = real_dim_x + 22;
	unsigned int real_dim_y_broad = real_dim_y + 22;
	unsigned int center_x = real_dim_x_broad / 2;
	unsigned int center_y = real_dim_y_broad / 2;
	unsigned int totalSamples = real_dim_x_broad * real_dim_y_broad;
	unsigned int dim_Samples = dim_x * dim_y;
	
	// Allocate memory in the GPU
	CHECK(cudaMalloc((void**)&dev_frame, dim_Samples * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_center_x, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_y, sizeof(unsigned int)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, dim_Samples * sizeof(float), cudaMemcpyHostToDevice));
	
	// Take the sharpCenter(x,y) and delete those rows and columns + or - 10 lines
	g_RemoveSharpCenter<<<11,11>>>(dev_frame, dim_x, dim_y, sharp_x, sharp_y);
	
	// Allocate memory in the GPU
	CHECK(cudaMalloc((void**)&dev_frame, dim_Samples * sizeof(float)));
	CHECK(cudaMalloc((float**)&dev_frame_long, totalSamples * sizeof(float)));
	
	CHECK(cudaMalloc((void**)&max, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_x, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_y, sizeof(unsigned int)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, dim_Samples * sizeof(float), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(dev_center_x, &center_x, sizeof(unsigned int), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(dev_center_y, &center_y, sizeof(unsigned int), cudaMemcpyHostToDevice));
	
	g_Zero_Array<<<real_dim_x, real_dim_y>>>(	dev_frame_long);
	
	g_Initialize_Array<<<dim_x, dim_y>>>(	dev_frame,
											dev_frame_long,
											real_dim_x_broad,
											real_dim_y_broad);
	
	// may need a cudaDeviceSynchronize here
	cudaDeviceSynchronize();
	CHECK(cudaGetLastError());
	
	g_Contracting_Max_Broad_Thrust<<<1,16>>>(dev_frame_long,
								max, 									
								dev_center_x,
								dev_center_y,
								real_dim_x_broad,
								real_dim_y_broad,
								iterations, 
								iteration_curr,
								shrinkFactor);
	
	CHECK(cudaGetLastError());
	// may need a cudaDeviceSynchronize here
	cudaDeviceSynchronize();
	
	// Copy the data back to the host
	CHECK(cudaMemcpy(&center_x, (void*)dev_center_x, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(&center_y, (void*)dev_center_y, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	
	broad_x = center_x - ((real_dim_x - dim_x) / 2);
	broad_y = center_y - ((real_dim_y - dim_y) / 2);
	
	// Free the data
	CHECK(cudaFree(dev_frame));
	CHECK(cudaFree(dev_frame_long));
	CHECK(cudaFree(max));
	CHECK(cudaFree(dev_center_x));
	CHECK(cudaFree(dev_center_y));	
}
*/



void Contracting_Grid::findBroadCenter_Rec()
{
	/*
	// This is the frame data that will be going into the GPU
	float *dev_frame;
	float *dev_frame_long;
	unsigned int *max;
	unsigned int *dev_center_x;
	unsigned int *dev_center_y;
	unsigned int iteration_curr = 0;
	unsigned int real_dim_x_broad = real_dim_x + 22;
	unsigned int real_dim_y_broad = real_dim_y + 22;
	unsigned int center_x = real_dim_x_broad / 2;
	unsigned int center_y = real_dim_y_broad / 2;
	unsigned int totalSamples = real_dim_x_broad * real_dim_y_broad;
	unsigned int dim_Samples = dim_x * dim_y;
	
	// Allocate memory in the GPU
	CHECK(cudaMalloc((void**)&dev_frame, dim_Samples * sizeof(float)));
	CHECK(cudaMalloc((void**)&dev_center_x, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_y, sizeof(unsigned int)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, dim_Samples * sizeof(float), cudaMemcpyHostToDevice));
	
	// Take the sharpCenter(x,y) and delete those rows and columns + or - 10 lines
	g_RemoveSharpCenter<<<11,11>>>(dev_frame, dim_x, dim_y, sharp_x, sharp_y);
	
	// Allocate memory in the GPU
	CHECK(cudaMalloc((void**)&dev_frame, dim_Samples * sizeof(float)));
	CHECK(cudaMalloc((float**)&dev_frame_long, totalSamples * sizeof(float)));
	
	CHECK(cudaMalloc((void**)&max, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_x, sizeof(unsigned int)));
	CHECK(cudaMalloc((void**)&dev_center_y, sizeof(unsigned int)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, dim_Samples * sizeof(float), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(dev_center_x, &center_x, sizeof(unsigned int), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(dev_center_y, &center_y, sizeof(unsigned int), cudaMemcpyHostToDevice));
	
	g_Zero_Array<<<real_dim_x, real_dim_y>>>(	dev_frame_long);
	
	g_Initialize_Array<<<dim_x, dim_y>>>(	dev_frame,
											dev_frame_long,
											real_dim_x_broad,
											real_dim_y_broad);
	
	// may need a cudaDeviceSynchronize here
	cudaDeviceSynchronize();
	CHECK(cudaGetLastError());
	
	
	//This part needs to change for this algorithm
	g_Contracting_Max_Broad_Thrust<<<1,16>>>(dev_frame_long,
								max, 									
								dev_center_x,
								dev_center_y,
								real_dim_x_broad,
								real_dim_y_broad,
								iterations, 
								iteration_curr,
								shrinkFactor);
	
	CHECK(cudaGetLastError());
	// may need a cudaDeviceSynchronize here
	cudaDeviceSynchronize();
	
	// Copy the data back to the host
	CHECK(cudaMemcpy(&center_x, (void*)dev_center_x, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(&center_y, (void*)dev_center_y, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	
	broad_x = center_x - ((real_dim_x - dim_x) / 2);
	broad_y = center_y - ((real_dim_y - dim_y) / 2);
	
	// Free the data
	CHECK(cudaFree(dev_frame));
	CHECK(cudaFree(dev_frame_long));
	CHECK(cudaFree(max));
	CHECK(cudaFree(dev_center_x));
	CHECK(cudaFree(dev_center_y));	
	*/
}


// referred to in Contracting_Grid.h
void Contracting_Grid::findDirection()
{
	direction_x = sharp_x - broad_x;
	direction_y = sharp_y - broad_y;
}

//Setter for the frame field
void Contracting_Grid::Frame(float *input_Frame)
{
	frame = new float[dim_x * dim_y];

	for(int i = 0; i<dim_x;i++)
	{
		for(int j = 0; j<dim_y;j++)
		{
			frame[MAP_2D(dim_x, dim_y, i, j)] =  input_Frame[MAP_2D(dim_x, dim_y, i, j)];
		}
	}
}

void Contracting_Grid::find_4D_Center()
{
	//find the sharpCenter for now, then add on the broadCenter later.
	
	clock_t t = clock();
	
	//t(x,y:x0,y0,x1,y1) = A*exp(-1/2[(x-x0)^2 + (y-y0)^2]/sigma^2) + A*exp(-1/2[(x-x1)^2 + (y-y1)^2]/sigma^2);
	//
	//real_dim_x += 4;
	//real_dim_y += 4;
	//make a new long frame
	float *frame_long = (float*)malloc(real_dim_x * real_dim_y *sizeof(float));
	
	for(int i = 0; i < real_dim_x * real_dim_y; i++)
	{
		frame_long[i] = 0.0f;
	}
	
	float A = 0;
	
	// Give an adequate sigma value for the Gaussian function
	// In the future, these sigma values should be evaluated in the constructor, based on the distance value given between surfaces.
	//float sigma_Sharp = 10;//3;
	//float sigma_Broad = 8;
	
	int extend_x = (real_dim_x - dim_x)/2;
	int extend_y = (real_dim_y - dim_y)/2;
	
	//int temp_max = 0;
	//search the frame and find the max value.  Set as A
	for(int i = 0; i < dim_x; i++)
	{
		for(int j = 0; j < dim_y; j++)
		{
			frame_long[MAP_2D(real_dim_x, real_dim_y, extend_x + i, extend_y + j)] = frame[MAP_2D(dim_x, dim_y, i, j)];
			if(A < frame[MAP_2D(dim_x,dim_y,i,j)])
			{
				A = frame[MAP_2D(dim_x,dim_y,i,j)];
			}
		}
	}
	
 	//float gauss[5][5];
	// Get the sharpCenter Gaussian array
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			gauss_Sharp[MAP_2D(5,5,i,j)] = A * gauss_Sharp[MAP_2D(5,5,i,j)];
		}
	} 
	
	cout<<"A: "<< A<< endl;
	
	// Get the original values for the contracting grid
	
	unsigned int center_x = real_dim_x/2;
	unsigned int center_y = real_dim_y/2;

	unsigned int *pix_x = new unsigned int[grid * grid];
	unsigned int *pix_y = new unsigned int[grid * grid];
	unsigned int idx_x = 0;
	unsigned int idx_y = 0;
	
	float *sum = new float[grid*grid];
	int max = 0;
	// Get the sharp peak values for the 
	for(int m = 0; m < iterations; m++)
	{

		
		//find new widths and 
				// Get the width of the current grid being applied
		unsigned int width_x = dim_x / (pow((double)shrinkFactor, (double)m));
		unsigned int width_y = dim_y / (pow((double)shrinkFactor, (double)m));
		
		// This should mark the top left spot for the square grid being used
		unsigned int start_x = (center_x) - width_x / 2; 
		unsigned int start_y = (center_y) - width_y / 2;	
		
 		/* for(int k = 0; k< grid * grid; k++){
			idx_x = k % 4;
			idx_y = k / 4;
			// Each thread will have a pix_x and pix_y in the frame which will be tested
			pix_x[k] = start_x + idx_x * width_x / 3;
			pix_y[k] = start_y + idx_y * width_y / 3;
			
			sum[k] = 0;
		}  */
	

		max = 0;
		for(int k = 0; k<grid*grid;k++)
		{	
			idx_x = k % 4;
			idx_y = k / 4;
			// Each thread will have a pix_x and pix_y in the frame which will be tested
			pix_x[k] = start_x + idx_x * width_x / (grid - 1);
			pix_y[k] = start_y + idx_y * width_y / (grid - 1);
			sum[k] = 0;
			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					
					// Given the values compare with the Gaussian curves
					sum[k] += pow(frame_long[MAP_2D(real_dim_x, real_dim_y, pix_x[k] + i - 2, pix_y[k] + j - 2)] - gauss_Sharp[MAP_2D(5,5,i,j)],2.0);
									
				}
			}

			// Choose the one where the magnitude is lowest.
			if(sum[k] < sum[max])
			{
				max = k;
				sharp_x = pix_x[k];
				sharp_y = pix_y[k];
				center_x = sharp_x;
				center_y = sharp_y;
			}	
				
		}
		
	}

	//cout<<"frame_long[x,y]: "<<frame_long[MAP_2D(real_dim_x,real_dim_y,sharp_x,sharp_y)]<<endl;
	
	sharp_x = sharp_x - ((real_dim_x - dim_x) / 2);
	sharp_y = sharp_y - ((real_dim_y - dim_y) / 2);
	
	cout<<"CPU: "<<endl;
	cout<<"sharp_x: "<< sharp_x << endl;
	cout<<"sharp_y: "<< sharp_y << endl;
	//Get time
	t = clock() - t;
	cout << "Time: " << (float)t/CLOCKS_PER_SEC << endl;
	
}

void Contracting_Grid::find_4D_Center_Cuda()
{
	//find the sharpCenter for now, then add on the broadCenter later.
	
	clock_t t = clock();
	clock_t t0 = clock();
	
	
	//int extend_x = (real_dim_x - dim_x)/2;
	//int extend_y = (real_dim_y - dim_y)/2;
	
	// Get the original values for the contracting grid
	
	//unsigned int center_x = real_dim_x/2;
	//unsigned int center_y = real_dim_y/2;
	//unsigned int iteration_curr = 0;
	//float *dev_frame_long;// = (float*)malloc(real_dim_x * real_dim_y *sizeof(float));
	//float *dev_frame;
	//float *dev_a;
	//float *dev_a_final;
	//unsigned int *dev_center_x;
	//unsigned int *dev_center_y;
	//float *pix_Compare;
	//float *dev_gauss_Sharp;
	
	//CHECK(cudaMalloc((void**)&dev_frame_long, real_dim_x * real_dim_y * sizeof(float)));
	//CHECK(cudaMalloc((void**)&dev_frame, dim_x * dim_y * sizeof(float)));
	//CHECK(cudaMalloc((void**)&dev_a, 128 * sizeof(float)));
	//CHECK(cudaMalloc((void**)&dev_a_final, sizeof(float)));
	//CHECK(cudaMalloc((void**)&dev_center_x, sizeof(unsigned int)));
	//CHECK(cudaMalloc((void**)&dev_center_y, sizeof(unsigned int)));
	//CHECK(cudaMalloc((void**)&pix_Compare, grid * grid * sizeof(float)));
	//CHECK(cudaMalloc((void**)&dev_gauss_Sharp, 25 * sizeof(float)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, dim_x * dim_y * sizeof(float), cudaMemcpyHostToDevice));
	
	t0 = clock() - t0;
	cout << "Time0: " << (float)t0/CLOCKS_PER_SEC << endl;
	clock_t t1 = clock();
	
	g_Zero_Array<<<1,128>>>(	dev_a);
	
	t1 = clock() - t1;
	cout << "Time1: " << (float)t1/CLOCKS_PER_SEC << endl;
	clock_t t2 = clock();
	
	// Set all dev_frame_long to zero
	g_Zero_Array<<<real_dim_x, real_dim_y>>>(	dev_frame_long);
	
	t2 = clock() - t2;
	cout << "Time2: " << (float)t2/CLOCKS_PER_SEC << endl;
	clock_t t3 = clock();
	
		// Put frame values into the frame_long array, see if a max value can be returned.
	g_Initialize_Array<<<dim_x, dim_y>>>(	dev_frame,
											dev_frame_long,
											real_dim_x,
											real_dim_y);
	
	t3 = clock() - t3;
	cout << "Time3: " << (float)t3/CLOCKS_PER_SEC << endl;
	clock_t t4 = clock();
	
	//search the frame and find the max value.  Set as A
	// Reduce: max
	unsigned int blocks = 128;
	const unsigned int blockSize = 256;//64 next
	const unsigned int blockSize1 = 64;
	g_reduce6_MinMax<blockSize><<<blocks, blockSize * 2, blockSize * 2 * sizeof(float)>>>(dev_frame, dev_a, dim_x * dim_y);
	
	//float *a = (float*)malloc(sizeof(float) * blocks);
	//CHECK(cudaMemcpy(a, dev_a, blocks * sizeof(float), cudaMemcpyDeviceToHost));
	
/* 	for(int i = 0; i < 128; i = i+4)
	{
		cout<< "a :"<< a[i] << "\t" << a[i+1] << "\t" << a[i+2] << "\t" << a[i+3] << endl;
	} */
	
	
	g_reduce6_MinMax<blockSize1><<<1, blocks, blocks * sizeof(float)>>>(dev_a, dev_a_final, blocks);
	
	t4 = clock() - t4;
	cout << "Time4: " << (float)t4/CLOCKS_PER_SEC << endl;
	clock_t t5 = clock();
	
	CHECK(cudaMemcpy(&A, dev_a_final, sizeof(float), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(dev_center_x, &center_x, sizeof(unsigned int), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(dev_center_y, &center_y, sizeof(unsigned int), cudaMemcpyHostToDevice));
	
	t5 = clock() - t5;
	cout << "Time5: " << (float)t5/CLOCKS_PER_SEC << endl;
	clock_t t6 = clock();
	
	//cout<<"got here: "<< 4 <<endl;
	//cout<<"A: " << A<< endl;
	
	
	
	/*
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			gauss_Sharp[MAP_2D(5,5,i,j)] = A * gauss_Sharp[MAP_2D(5,5,i,j)];
		}
	} 
	*/
	
	//cout<<"got here: "<< 5 <<endl;
	
	CHECK(cudaMemcpy(dev_gauss_Sharp, gauss_Sharp, 25 * sizeof(float), cudaMemcpyHostToDevice));
	
	float *sum = new float[grid*grid];
	//int max = 0;
	unsigned int grid_squared = grid * grid;
	unsigned int shared_Size =  sizeof(float) * grid_squared;
	
	t6 = clock() - t6;
	cout << "Time6: " << (float)t6/CLOCKS_PER_SEC << endl;
	
	clock_t t7 = clock();
	
	//cout<<"got here: "<< 6 <<endl;
	//cout<<"iterations:  "<< iterations<<endl;
	g_getPeaks<<<1, grid_squared, shared_Size>>>(	dev_frame_long,
													dev_gauss_Sharp, 	
													pix_Compare,
													dev_center_x,	// X-dimension center of the enlarged frame
													dev_center_y, 	// Y-dimension center of the enlarged frame
													real_dim_x,
													real_dim_y,
													iterations, 
													iteration_curr,
													grid,
													shrinkFactor);

	t7 = clock() - t7;
	cout << "Time7: " << (float)t7/CLOCKS_PER_SEC << endl;
	clock_t t8 = clock();	
	
	//cout<<"got here: "<< 7 <<endl;					
	//cout<<"frame_long[x,y]: "<<frame_long[MAP_2D(real_dim_x,real_dim_y,sharp_x,sharp_y)]<<endl;
	// Copy the data back to the host
	CHECK(cudaMemcpy(&center_x, (void*)dev_center_x, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(&center_y, (void*)dev_center_y, sizeof(unsigned int), cudaMemcpyDeviceToHost));
	
	sharp_x = center_x - ((real_dim_x - dim_x) / 2);
	sharp_y = center_y - ((real_dim_y - dim_y) / 2);
	

	//Get time

	t8 = clock() - t8;
	
	cout << "Time8: " << (float)t8/CLOCKS_PER_SEC << endl;
	
	t = clock() - t;
	cout << "Time: " << (float)t/CLOCKS_PER_SEC << endl;
	
	cout<<"CUDA: "<<endl;
	cout<<"sharp_x: "<< sharp_x << endl;
	cout<<"sharp_y: "<< sharp_y << endl;
	
	CHECK(cudaFree(dev_frame_long));
	CHECK(cudaFree(dev_frame));
	CHECK(cudaFree(dev_a));
	CHECK(cudaFree(dev_a_final));
	CHECK(cudaFree(dev_center_x));
	CHECK(cudaFree(dev_center_y));
	CHECK(cudaFree(pix_Compare));
	CHECK(cudaFree(dev_gauss_Sharp));
	
}





