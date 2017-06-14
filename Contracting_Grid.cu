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
		
#include <iostream>
#include "Contracting_Grid.h"
//#include "Cine.cu"
		
using namespace std;
		
#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))

#define CHECK(x) do {\
	cudaError_t err = (x);\
	if ( err!= cudaSuccess)\
	{\
		fprintf(stderr, "API error %s:%d Returned:%d\n",__FILE__, __LINE__, err);\
		exit(1);\
	}while(0)
		
// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid()
{
	// standard assumptions for algorithms
	iterations = 7;
	shrinkFactor = 2;
	distance = .075;
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = tempSquareCenter = 0;
	direction_x = direction_y  = real_dim_x = real_dim_y = dim_x = dim_y = 0;
	
	// May need to instantiate 'frame'
	frame = NULL;
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c)
{
	// standard assumptions for algorithms
	iterations = 7;
	shrinkFactor = 2;
	distance = .075;
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = tempSquareCenter = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c, unsigned int iterations_in, float shrinkFactor_in)
{
	// standard assumptions for algorithms			
	distance = .075;
	
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = tempSquareCenter = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	//Checks that the value for iterations is acceptable
	if(iterations_in > 0)
	{
		iterations = iterations_in;
	}
	
	//Checks that the value for the shrinkFactor is acceptable
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}
}

// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(Cine c, unsigned int iterations_in, float shrinkFactor_in, float distance_in)
{
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = tempSquareCenter = 0;
	direction_x = direction_y = 0;
	
	// Get dimensions of the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	//Checks that the value for iterations is acceptable
	if(iterations_in > 0)
	{
		iterations = iterations_in;
	}
	
	//Checks that the value for the shrinkFactor is acceptable
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}

	//Checks that the value for the distance is acceptable
	if(distance_in > 0)
	{
		distance = distance_in;
	}
}


// referred to in Contracting_Grid.h
Contracting_Grid::Contracting_Grid(unsigned int iterations_in, float shrinkFactor_in, float distance_in)
{
	// Set these numbers to 0.
	sharp_x = sharp_y = broad_x = broad_y = tempSquareCenter = 0;
	direction_x = direction_y = real_dim_x = real_dim_y = 0;
	
	// The real dimensions of the frame will be received later if this is the contructor used
	dim_x = dim_y = 0;
	
	frame = NULL;
	//Checks that the value for iterations is acceptable
	if(iterations_in > 0)
	{
		iterations = iterations_in;
	}
	
	//Checks that the value for the shrinkFactor is acceptable
	if(shrinkFactor_in > 1)
	{
		shrinkFactor = shrinkFactor_in;
	}

	//No check is needed yet
	distance = distance_in;	
}

 // referred to in Contracting_Grid.h
/* void Contracting_Grid::getFrameData(&Cine c)
{
	// Get the x and y dimensions to know how much memory to allocate for the frame
	dim_x = c.Dim_x();
	dim_y = c.Dim_y();
	
	real_dim_x = dim_x + (dim_x / shrinkFactor);
	real_dim_y = dim_y + (dim_y / shrinkFactor);
	
	// Allocate memory 
	//frame = new float[dim_x * dim_y];
	malloc(frame, real_dim_x * real_dim_y * sizeof(float));
	// Copy the values over
	for(i = 0; i < dim_x; ++i) {
	  for(j = 0; j < dim_y; ++j) {
		//frame[MAP_2D(dim_x, dim_y, i, j,)] = c.data_frame[MAP_2D(dim_x, dim_y, i, j)];
		frame[MAP_2D(dim_x, dim_y, dim_x/(2*shrinkFactor) + i, dim_y/(2*shrinkFactor) + j)] = c.data_frame[MAP_2D(dim_x, dim_y, i, j)];
	  }
	}
}  */

// Essence of the contracting grid algorithm and the purpose of this class.  Referred to in Contracting_Grid.h
void Contracting_Grid::findSharpCenter()
{
	/*
	// this is the frame data that will be going into the GPU
	float *dev_frame;
	
	// Allocate memory in the GPU
	CHECK(cudaMalloc((void**)&frame, sizeof(float)));
	
	// Copy the data to the GPU
	CHECK(cudaMemcpy(dev_frame, frame, real_dim_x * real_dim_y * sizeof(float), cudaMemcpyHostToDevice));
	
	// Copy the data back to the host
	CHECK(cudaMemcpy(frame, dev_frame, real_dim_x * real_dim_y * sizeof(float), cudaMemcpyDeviceToHost));
	
	// Free the data
	cudaFree(dev_frame);
	*/
}

// referred to in Contracting_Grid.h
void Contracting_Grid::findBroadCenter(){}

// referred to in Contracting_Grid.h
void Contracting_Grid::findDirection(){}

// Cine Contracting_Grid::operator=(const Cine &c)
// {
	// Cine cine();
	// cine.Filename(c.Filename());
	// cine.Dim_x(c.Dim_x());
	// cine.Dim_y(c.Dim_y());
	// cine.Dim_z(c.Dim_z());
	// cine.read_cine_3d();
// }

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
