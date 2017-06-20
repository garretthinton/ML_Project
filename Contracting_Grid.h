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

#ifndef CONTRACTING_GRID_H
#define CONTRACTING_GRID_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "Cine.h"
#include "Device.h"

#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))

//class Cine;

class Contracting_Grid
{
	private:
		// The amount of iterations of the Contracting Grid algorithm to be run.
		unsigned int iterations;
		
		// The factor to shrink the square in the Contracting Grid algorithm. 
		//   Or in other words: Size_next = Size / shrinkFactor
		float shrinkFactor;
		
		// The distance between surfaces in the system used 
		float distance;
		
		// The image frame to be analyzed
		float *frame;
				
		// The point with the max value in the current iteration of the Contracting Grid algorithm
		//unsigned int tempSquareCenter;
		
		// The sharpCenter pixel data points
		unsigned int sharp_x, sharp_y;
		
		// The broadCenter pixel data points
		unsigned int broad_x, broad_y;
		
		// The direction of propagation for the beta particle
		float direction_x, direction_y;
		
		// The expanded grid's x-dimension extents
		unsigned int real_dim_x;
		
		// The expanded grid's y-dimension extents
		unsigned int real_dim_y;
		
		// The number of pixels in the x-dimension of the frame
		unsigned int dim_x;
		
		// The number of pixels in the y-dimension of the frame
		unsigned int dim_y;
		
	public:
	
		/*
		Title: Contracting_Grid
		Description:  The Default Constructor.  This makes some assumptions about the imaging system, such as the distance
			in between the surfaces.  
		Inputs:		NONE
		Outputs:	NONE
		*/
		Contracting_Grid();
		
		/*
		Title: Contracting_Grid
		Description:  An overloaded Constructor.  This gets the dim_x, dim_y, and frame data from the Cine class and makes some assumptions about the rest of the class, such as the distance in between the surfaces and algorithm details.  
		Inputs:  
			-Cine c:  This is a friend member defined in Cine.h.
		Outputs:	NONE
		*/
		Contracting_Grid(	Cine c	);
		
		/*
		Title: Contracting_Grid
		Description:  An overloaded Constructor.  This gets the dim_x, dim_y, and frame data from the Cine class and makes some assumptions about the rest of the class, such as the distance in between the surfaces.  
		Inputs:  
			-Cine c:  This is a friend member defined in Cine.h.
			-float shrinkFactor:  How much to shrink the square grid with each iteration.
		Outputs:	NONE
		*/
		Contracting_Grid(	Cine c, 
							float shrinkFactor_in);
		
		/*
		Title: Contracting_Grid
		Description:  An overloaded Constructor.  This gets the dim_x, dim_y, and frame data from the Cine class and makes some assumptions about the rest of the class, such as the distance in between the surfaces.  
		Inputs:  
			-Cine c:  This is a friend member defined in Cine.h.
			-int iterations:	This is the amount of iterations to be run on the algorithm.
			-float shrinkFactor:  How much to shrink the square grid with each iteration.
		Outputs:	NONE
		*/
		Contracting_Grid(	Cine c, 
							unsigned int iterations,
							float shrinkFactor);
							
		/*
		Title: Contracting_Grid
		Description:  An overloaded Constructor.  This gets the dim_x, dim_y, and frame data from the Cine class and gets all other needed data from the constructor.  
		Inputs:  
			-Cine c:  This is a friend member defined in Cine.h.
			-float shrinkFactor:  How much to shrink the square grid with each iteration.
			-const int distance_in:	the distance between surfaces in the imaging system.
		Outputs:	NONE
		*/
		Contracting_Grid(	Cine c, 
							float shrinkFactor_in,
							float distance_in);
		
		/*
		Title: Contracting_Grid
		Description:  An overloaded Constructor.  This gets the dim_x, dim_y, and frame data from the Cine class and gets all other needed data from the constructor.  
		Inputs:  
			-Cine c:  This is a friend member defined in Cine.h.
			-int iterations:	This is the amount of iterations to be run on the algorithm.
			-float shrinkFactor:  How much to shrink the square grid with each iteration.
			-const int distance_in:	the distance between surfaces in the imaging system.
		Outputs:	NONE
		*/
		Contracting_Grid(	Cine c, 
							unsigned int iterations,
							float shrinkFactor,
							float distance_in);
		
		/*
		Title: Contracting_Grid
		Description:  An overloaded Constructor.  Gets all needed data from the constructor, except for the frame data and dimension data.  
		Inputs:  
			-int iterations:	This is the amount of iterations to be run on the algorithm.
			-float shrinkFactor:  How much to shrink the square grid with each iteration.
			-const int distance_in:	the distance between surfaces in the imaging system.
		Outputs:	NONE
		*/
		Contracting_Grid(	unsigned int iterations,
							float shrinkFactor,
							float distance_in);
		
		/*
		Title: findSharpCenter
		Description:  This finds the sharp peak that corresponds to a beta particle on the second surface.  This is accomplished using the Contracting_Grid algorithm and CUDA code.  GPU's are needed for this function. 
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findSharpCenter_Rec();
		
		/*
		Title: findBroadCenter
		Description:  This finds the smooth, broad peak that corresponds to a beta particle on the first surface.  This is accomplished using the Contracting_Grid algorithm and CUDA code.  GPU's are needed for this function. 
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findBroadCenter();
		
		/*
		Title: findDirection
		Description:  This finds the direction of the beta particle being imaged.  This is done by subtracting the sharpCenter point from the broadCenter point and using trig, 
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findDirection();
	
		/*
		Title: getFrameData
		Description:  This function takes a Cine class, and gets three private members of the class to read in:  data_frame, dim_x, and dim_y. 
		Inputs:  	
			-Cine c:  This class has the getFrameData function as a friend, and therefore this function can access the needed private variables.
		Outputs:	NONE
		*/
		//void getFrameData(	Cine c	);
		/* void getFrameData(Cine c)
		{
			// Get the x and y dimensions to know how much memory to allocate for the frame
			dim_x = c.Dim_x();
			dim_y = c.Dim_y();
			
			real_dim_x = dim_x + (dim_x / shrinkFactor);
			real_dim_y = dim_y + (dim_y / shrinkFactor);
			
			// Allocate memory 
			//frame = new float[dim_x * dim_y];
			frame = (float *)malloc( real_dim_x * real_dim_y * sizeof(float));
			// Copy the values over
			for(int i = 0; i < dim_x; ++i) {
			  for(int j = 0; j < dim_y; ++j) {
				//frame[MAP_2D(dim_x, dim_y, i, j,)] = c.data_frame[MAP_2D(dim_x, dim_y, i, j)];
				frame[MAP_2D(dim_x, dim_y, dim_x/(2*shrinkFactor) + i, dim_y/(2*shrinkFactor) + j)] = c.data_frame[MAP_2D(dim_x, dim_y, i, j)];
			  }
			}
		}  */
		
		/*
		Title:  Operator overload of assignment operator.
		Description:  This will be used in some of the constructors.
		Inputs:
			const Cine &c: This is the Cine that needs to be copied.
		Outputs:
			-This operator returns a Cine class.
		*/
		//Cine operator=(const Cine &c);
		
		void Frame(float* frame_in);
};



#endif