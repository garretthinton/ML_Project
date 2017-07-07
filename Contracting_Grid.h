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
#include <time.h>
#include "Cine.h"
#include "Device.h"
#include <algorithm>

#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))

//class Cine;

class Contracting_Grid
{
	private:
		// The amount of iterations of the Contracting Grid algorithm to be run.
		unsigned int iterations;
		
		// The factor to shrink the square in the Contracting Grid algorithm. 
		//   Or in other words: Size_next = Size / shrinkFactor^iteration
		float shrinkFactor;
		
		// The distance between surfaces in the system used 
		float distance;
		
		// The image frame to be analyzed
		float *frame;
		
		// stores the 2D gaussian functions to be used to check for the Sharp and Broad peaks.
		float *gauss_Sharp, *gauss_Broad;
		
		float *dev_frame, *dev_frame_long;
			
		float *dev_a;
		float *dev_a_final;
		unsigned int *dev_center_x;
		unsigned int *dev_center_y;
		float *pix_Compare;
		float *dev_gauss_Sharp;
		float A;
		
		
		unsigned int center_x; 
		unsigned int center_y; 
		unsigned int iteration_curr;	
		//the sigmas needed for the gaussian shapes.  
		float sigma_Sharp, sigma_Broad;
		
		// The sharpCenter pixel data points
		unsigned int sharp_x, sharp_y;
		
		// The broadCenter pixel data points
		unsigned int broad_x, broad_y;
		
		// The direction of propagation for the beta particle
		float direction_x, direction_y;
		
		// The expanded grid's x and y-dimension extents
		unsigned int real_dim_x, real_dim_y;
		
		// The number of pixels in the x and y-dimension of the frame
		unsigned int dim_x, dim_y;
		
		// The number of points on a side of the square grid.  The total number of points would be grid^2.
		unsigned int grid;
		
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
		Description:  This finds the sharp peak that corresponds to a beta particle on the second surface.  This is accomplished using the Contracting_Grid algorithm and CUDA code.  GPU's are needed for this function.  This function is recursive.
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findSharpCenter_Rec();
		
		/*
		Title: findBroadCenter_brute
		Description:  This uses a brute force method and finds the smooth, broad peak that corresponds to a beta particle on the first surface.  GPU's are needed for this function. 
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findBroadCenter_brute();
		
		/*
		Title: findBroadCenter_Rec
		Description:  This uses a modified Contracting Grid algorithm to find the broad peak in the imager.   GPU's are needed for this function.  This function uses only programmer code to maximize the code's effectiveness.
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findBroadCenter_Rec_Thrust();
		
				/*
		Title: findBroadCenter_Rec
		Description:  This uses a modified Contracting Grid algorithm to find the broad peak in the imager.   GPU's are needed for this function.  This function uses the Thrust library to simplify the code.
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findBroadCenter_Rec();
		
		
		/*
		Title: findDirection
		Description:  This finds the direction of the beta particle being imaged.  This is done by subtracting the sharpCenter point from the broadCenter point and using trig, 
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void findDirection();
		
		/*
		Title:  Frame
		Description:  Reads in any 2D image that is stored as a float*, and stores it as the float* frame in the Contracting_Grid class.
		Inputs:
			-float* frame_in:  The input frame that is to be turned into the frame field in the Contracting_Grid class.
		Outputs:	NONE
		*/
		void Frame(float* frame_in);
		
		/*
		Title: find_4D_center
		Description:  Uses a 4D contracting grid method to find both the broad and sharp peaks.  
		Inputs:  	NONE
		Outputs:	NONE
		*/
		void find_4D_Center();
		
		//Cuda version
		void find_4D_Center_Cuda();
};



#endif