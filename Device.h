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
// Title:	****Device.h****
// Description:  This file is for functions ran on the GPU, whether they be kernals, or normal __device__ functions.  

#ifndef DEVICE_H
#define DEVICE_H

#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/extrema.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
//#include <arrayfire.h>

/*
Title: g_Initialize_Array
Description:  Takes the array copies in the frame data to a larger 2D array.  The data is copied so that it is in the center of the bigger 2D array.
Inputs:		
	-float *dev_frame: 			This is the frame data that has been copied onto the GPU.
	-float *dev_frame_long:		This is the bigger array that the dev_frame data will be copied onto.  Therefore, this is an output.
	-unsigned int real_dim_x:	This is the x-dimension of the bigger array to be used as a reference.
	-unsigned int real_dim_y:	This is the y-dimension of the bigger array to be used as a reference.
Outputs:	
	-float *dev_frame_long:  	Mentioned in inputs.
*/
__global__ void g_Initialize_Array(	float *dev_frame,
									float *dev_frame_long,
									unsigned int real_dim_x,
									unsigned int real_dim_y);

/*
Title: g_Contracting_Max
Description: 	This is a recursive implementation of the Contracting Grid algorithm.  It takes in different parameters and outputs the 
				max value's coordinates on the dev_frame_long variable.
Inputs:		
	-float *dev_frame_long:			This is the bigger array that the dev_frame data will be copied onto.  Therefore, this is an output.
	-unsigned int *max:				This is the max value of each iteration of the Contracting Grid algorithm.  The greatest value will be returned after the last iteration.
	-unsigned int *center_x:		This is the center of the x-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int *center_y:		This is the center of the y-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int real_dim_x:		This is the x-dimension of the bigger array to be used as a reference.
	-unsigned int real_dim_x:		This is the x-dimension of the bigger array to be used as a reference.
	-unsigned int iteration_num:	This is the total number of iterations that need to be performed.
	-unsigned int iteration_curr:	This is the current iteration, or recursive call.
	-float shrinkFactor:			This is the amount that the square grid will be diminished in future iterations.
Outputs:	
	-float *dev_frame_long:  	Mentioned in inputs.
	-unsigned int *max:			""
	-unsigned int *center_x:	""
	-unsigned int *center_y:	""
*/
__global__ void g_Contracting_Max(	float *dev_frame_long,
									unsigned int *max, 									
									unsigned int *center_x,
									unsigned int *center_y,
									unsigned int real_dim_x,
									unsigned int real_dim_y,
									unsigned int iteration_num, 
									unsigned int iteration_curr,
									float shrinkFactor);

/*
Title: g_Contracting_Max_Broad
Description: 	This is a special recursive implementation of the Contracting Grid algorithm.  It takes in different parameters and outputs the 
				broad center's coordinates on the dev_frame_long variable.
Inputs:		
	-float *dev_frame_long:			This is the bigger array that the dev_frame data will be copied onto.  Therefore, this is an output.
	-unsigned int *max:				This is the max value of each iteration of the Contracting Grid algorithm.  The greatest value will be returned after the last iteration.
	-unsigned int *center_x:		This is the center of the x-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int *center_y:		This is the center of the y-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int real_dim_x:		This is the x-dimension of the bigger array to be used as a reference.
	-unsigned int real_dim_y:		This is the y-dimension of the bigger array to be used as a reference.
	-unsigned int iteration_num:	This is the total number of iterations that need to be performed.
	-unsigned int iteration_curr:	This is the current iteration, or recursive call.
	-float shrinkFactor:			This is the amount that the square grid will be diminished in future iterations.
Outputs:	
	-float *dev_frame_long:  	Mentioned in inputs.
	-unsigned int *max:			""
	-unsigned int *center_x:	""
	-unsigned int *center_y:	""
*/									
__global__ void g_Contracting_Max_Broad_Thrust(	float *dev_frame_long,
												unsigned int *max, 									
												unsigned int *center_x,
												unsigned int *center_y,
												unsigned int real_dim_x,
												unsigned int real_dim_y,
												unsigned int iteration_num, 
												unsigned int iteration_curr,
												float shrinkFactor);									
									
/*
Title: g_Contracting_Max_Broad
Description: 	This is a special recursive implementation of the Contracting Grid algorithm.  It takes in different parameters and outputs the 
				broad center's coordinates on the dev_frame_long variable.
Inputs:		
	-float *dev_frame_long:			This is the bigger array that the dev_frame data will be copied onto.  Therefore, this is an output.
	-unsigned int *max:				This is the max value of each iteration of the Contracting Grid algorithm.  The greatest value will be returned after the last iteration.
	-unsigned int *center_x:		This is the center of the x-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int *center_y:		This is the center of the y-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int real_dim_x:		This is the x-dimension of the bigger array to be used as a reference.
	-unsigned int real_dim_y:		This is the y-dimension of the bigger array to be used as a reference.
	-unsigned int iteration_num:	This is the total number of iterations that need to be performed.
	-unsigned int iteration_curr:	This is the current iteration, or recursive call.
	-float shrinkFactor:			This is the amount that the square grid will be diminished in future iterations.
Outputs:	
	-float *dev_frame_long:  	Mentioned in inputs.
	-unsigned int *max:			""
	-unsigned int *center_x:	""
	-unsigned int *center_y:	""
*/									
__global__ void g_Contracting_Max_Broad(	float *dev_frame_long,
											unsigned int *max, 									
											unsigned int *center_x,
											unsigned int *center_y,
											unsigned int real_dim_x,
											unsigned int real_dim_y,
											unsigned int iteration_num, 
											unsigned int iteration_curr,
											float shrinkFactor);

/*
Title: g_Contracting_Max_Broad
Description: 	This is a special recursive implementation of the Contracting Grid algorithm.  It takes in different parameters and outputs the 
				broad center's coordinates on the dev_frame_long variable.
Inputs:		
	-float *dev_frame_long:			This is the bigger array that the dev_frame data will be copied onto.  Therefore, this is an output.
	-unsigned int *max:				This is the max value of each iteration of the Contracting Grid algorithm.  The greatest value will be returned after the last iteration.
	-unsigned int *center_x:		This is the center of the x-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int *center_y:		This is the center of the y-dimension of the square grid that is being used to calculate the max value.  This is in terms of the bigger array.
	-unsigned int real_dim_x:		This is the x-dimension of the bigger array to be used as a reference.
	-unsigned int real_dim_y:		This is the y-dimension of the bigger array to be used as a reference.
	-unsigned int iteration_num:	This is the total number of iterations that need to be performed.
	-unsigned int iteration_curr:	This is the current iteration, or recursive call.
	-float shrinkFactor:			This is the amount that the square grid will be diminished in future iterations.
Outputs:	
	-float *dev_frame_long:  	Mentioned in inputs.
	-unsigned int *max:			""
	-unsigned int *center_x:	""
	-unsigned int *center_y:	""
*/									
__global__ void g_Contracting_Max_Broad_ArrayFire(	float *dev_frame_long,
													unsigned int *max, 									
													unsigned int *center_x,
													unsigned int *center_y,
													unsigned int real_dim_x,
													unsigned int real_dim_y,
													unsigned int iteration_num, 
													unsigned int iteration_curr,
													float shrinkFactor);											
									
/*
Title: g_Zero_Array
Description:  Takes the array and initializes it to 0.
Inputs:		
	-float *dev_frame: 			This is the frame data that has been copied onto the GPU.
	-float *dev_frame_long:		This is the bigger array that the dev_frame data will be copied onto.  Therefore, this is an output.
	-unsigned int real_dim_x:	This is the x-dimension of the bigger array to be used as a reference.
	-unsigned int real_dim_y:	This is the y-dimension of the bigger array to be used as a reference.
Outputs:	
	-float *dev_frame_long:  	Mentioned in inputs.
*/
__global__ void g_Zero_Array(		float *dev_frame_long);

/*
Title: g_GetValues
Description:  This is the beginning of the algorithm that returns the broadCenter data points.  It takes in the frame, which dimension (x or y) to calculate, and a 1D float* array where the sum of each row in the frame is stored. The 1D array is then returned for the algorithm to use in its next phase.
Inputs:		
	-float *dev_frame: 	This is the frame data that has been copied onto the GPU.
	-int dim:			This tells the algorithm whether to calculate the x or y direction.
	-float *dev_x:		This is the output array that is the sum of each value in the rows/columns.
Outputs:	
	-float *dev_x:  	Mentioned in inputs.
*/									
__global__ void g_GetValues(		float *dev_frame,	// This is the full 2D array
									int dim,			// 0 is x, 1 is y
									float *dev_x);		// This is the output array

/*
Title: g_GetCenters
Description:  This is the second half of the getBroadCenter algorithm.  It takes the sums that were calculated in g_GetValues, dev_x, sums up 21 values at a time from dev_x, and puts them in another array.  From that array, the highest sum is taken and outputted as dev_broad_x.
Inputs:		
	-float *dev_x:				This is the array that is the sum of each value in the rows/columns.
	-unsigned int *dev_broad_x: This is the output value that gives the center of the broad peak. 
Outputs:	
	-unsigned int *dev_broad_x: Mentioned in inputs.
*/									
__global__ void g_GetCenters(		float *dev_x,
									unsigned int *dev_broad_x);

/*
Title: g_RemoveSharpCenter
Description:  Takes the coordinates for the sharp peak's center, and replaces its value, as well as the value of 5 pixels all around it, with pixels surrounding it.
Inputs:		
	-float *dev_frame:				The frame being read from and modified.
	-unsigned int sharp_x: 			The x coordinate of the sharp peak.
	-unsigned int sharp_y: 			The y coordinate of the sharp peak.
Outputs:	
	-float *dev_frame: Mentioned in inputs.
*/										
__global__ void g_RemoveSharpCenter(	float *dev_frame, 
										unsigned int dim_x,
										unsigned int dim_y,
										unsigned int sharp_x,
										unsigned int sharp_y);									

/*
Title: g_SumBroad
Description:  Takes in a 2D array.  Given the pixel points passed in, it sums the values around that pixel and sums it up into a value that is returned.
Inputs:		
	-float *dev_frame_long:	The 2D frame stored in a float pointer.
	-float *pixel:			Pixel array that is in shared memory on the getBroadCenter_Rec kernal.
	-unsigned int pix_x:	X-dimension index of the pixel to calculate.
	-unsigned int pix_y:	Y-dimension index of the pixel to calculate.	
	-unsigned int idx_x:	X index of the float* pixel to store the sum into.
	-unsigned int idx_y:	Y index of the float* pixel to store the sum into.
Outputs:	
	-float *pixel: 			Mentioned in inputs.
*/
__global__ void g_PutInArray(	float *dev_frame_long,
								float *dev_array,
								unsigned int pix_x,
								unsigned int pix_y);									
									
__global__ void g_PutInArray_1(	float *dev_frame_long,
								float *dev_array,
								unsigned int pix_x,
								unsigned int pix_y);									
									
#endif