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
		
// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid()
{
	// standard assumptions for algorithms
	iterations = 7;
	shrinkFactor = 2;
	distance = .075;
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = 0;
	direction_x = direction_y  = real_dim_x = real_dim_y = dim_x = dim_y = 0;
	
	// May need to instantiate 'frame'
	frame = NULL;
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c)
{
	// standard assumptions for algorithms
	shrinkFactor = 2;
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

// Essence of the contracting grid algorithm and the purpose of this class.  Referred to in Contracting_Grid.h
void Contracting_Grid::findSharpCenter_Rec()
{
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
	
	g_Zero_Array<<<real_dim_x, real_dim_y>>>(	dev_frame,
												dev_frame_long,
												dim_x,
												dim_y);
	
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
	
	
	/*
	cout<< "sharp_x: " << sharp_x << endl;
	cout<< "sharp_y: " << sharp_y << endl;
	
	cout<< "frame[sharp_x, sharp_y]: "		<<  frame[MAP_2D(dim_x,dim_y, sharp_x, sharp_y)]<<endl;
	cout<< "frame[sharp_x-1, sharp_y]: "	<<    frame[MAP_2D(dim_x,dim_y,sharp_x - 1, sharp_y)]<<endl;
	cout<< "frame[sharp_x+1, sharp_y]: "	<<    frame[MAP_2D(dim_x,dim_y,sharp_x + 1, sharp_y)]<<endl;
	cout<< "frame[sharp_x, sharp_y-1]: "	<<   frame[MAP_2D(dim_x,dim_y,sharp_x, sharp_y - 1)]<<endl;
	cout<< "frame[sharp_x-1, sharp_y-1]: "	<<   frame[MAP_2D(dim_x,dim_y,sharp_x - 1, sharp_y - 1)]<<endl;
	cout<< "frame[sharp_x+1, sharp_y-1]: "	<<    frame[MAP_2D(dim_x,dim_y,sharp_x + 1, sharp_y - 1)]<<endl;
	cout<< "frame[sharp_x, sharp_y+1]: "	<<   frame[MAP_2D(dim_x,dim_y,sharp_x, sharp_y + 1)]<<endl;
	cout<< "frame[sharp_x-1, sharp_y+1]: "	<<    frame[MAP_2D(dim_x,dim_y,sharp_x - 1, sharp_y + 1)]<<endl;
	cout<< "frame[sharp_x+1, sharp_y+1]: "	<<   frame[MAP_2D(dim_x,dim_y,sharp_x + 1, sharp_y + 1)]<<endl; 
	*/
	
	/*
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
	*/
	
	// Free the data
	CHECK(cudaFree(dev_frame));
	CHECK(cudaFree(dev_frame_long));
	CHECK(cudaFree(max));
	CHECK(cudaFree(dev_center_x));
	CHECK(cudaFree(dev_center_y));	
}

// referred to in Contracting_Grid.h
void Contracting_Grid::findBroadCenter_brute()
{	

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