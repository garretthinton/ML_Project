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
// Garrett Hinton		Jun 12 2017
//
//
// Description of File
// This is a worker class that takes a .cine file and turns it into an array of frames. 
// It also has the ability to turn a .cine frame into a .eps file (Latec).  This was taken 
// Luca's code.

#ifndef CINE_H
#define CINE_H

//#include "Contracting_Grid.h"

class Cine
{
	private:
		
		// These are the variables for the x, y, and z dimensions.
		// The x and y are dimensions for each frame, while the z variable shows how 
		// 	 many frames there are.
		unsigned int dim_x, dim_y, dim_z;
		
		// The filename of the .cine file being used.
		const char *filename;
		
		// This stores the raw data of the .cine file
		float *data_all;
		float *data_frame;
	
	public:
	
		/*
		Title: Default constructor
		Description:  This constructor needs no parameters modified
		Inputs:		None
		Outputs:	None
		*/
		Cine();
	
		/*
		Title: Overloaded Constructor
		Description:  This constructor makes a Cine Class that is functional and ready to use
		Inputs:		
			-char *filename : This is the name of the file to be read.
		Outputs:	None
		*/
		Cine(const char *filename);
	
		/*
		Title: write_eps
		Description:  Writes out one frame of data from the 'data' variable to the file named,
			'filename'.  This output file is in .eps format 
		Inputs:  
			-float *data:  Though this is an input parameter, this is actually the output
			-const unsigned int dim_x:  The x dimension of the frame
			-const unsigned int dim_y:  The y dimension of the frame
			-const char *filename:  The name of the .cine file to be written.
		Outputs:
			-None, as the output is in a file.
		*/
		void write_eps(	 const char *filename_out  );
		
		/*
		Title: write_dat_3d
		Description:  Writes out the data from the 'data' variable to the file named,
			'filename'.
		Inputs:  
			-float *data:  Though this is an input parameter, this is actually the output
			-const unsigned int dim_x:  The x dimension of the frame
			-const unsigned int dim_y:  The y dimension of the frame
			-const unsigned int dim_z:  The z dimension, which is how many total 
				frames there are.
			-const char *filename:  The name of the .cine file to be written.
		Outputs:
			-None, as the output is in a file.
		*/
		void write_dat_3d(	const char *filename_out  ); 
	
		/*
		Title: read_cine_3d
		Description:  Outputs all frames from the .cine file onto the variable 'data'.  This
			is taken from the .cine file, 'filename', and is frame number 'index'.
		Inputs:  
			-float *data:  Though this is an input parameter, this is actually the output
			-const unsigned int dim_x:  The x dimension of the frame
			-const unsigned int dim_y:  The y dimension of the frame
			-const unsigned int dim_z:  The z dimension, which is how many total 
				frames there are.
			-const char *filename:  The name of the .cine file to be read.
		Outputs:
			-None:  But the true output is the float *data variable.
		*/
		void read_cine_3d(  );
	
		/*
		Title: read_cine_info
		Description:  Inputs the filename, and outputs the x, y, and z dimensions of the 
			.cine file.
		Inputs:  
			-unsigned int dim_x:  The x dimension of the frame
			-unsigned int dim_y:  The y dimension of the frame
			-unsigned int dim_z:  The z dimension, which is how many total frames there are.
			-const char *filename:  The name of the .cine file to be read.
		Outputs:
			-None:  But the true outputs are the dim_* variables.
		*/
		void read_cine_info(	);
		
		/*
		Title: read_cine_frame
		Description:  Outputs a frame from the .cine file onto the variable 'data'.  This frame  
			is taken from the .cine file, 'filename', and is frame number 'index'.
		***Need to be sure to check that index is less that dim_z
		Inputs:  
			-float *data:  Though this is an input parameter, this is actually the output
			-const unsigned int dim_x:  The x dimension of the frame
			-const unsigned int dim_y:  The y dimension of the frame
			-const unsigned int index:  The index of the frame, or in other words, the 	
				frame number.
			-const char *filename:  The name of the .cine file to be read.
		Outputs:
			-float *:  It will output the data frame.  But the true output is the float *data variable.
		*/
		float* read_cine_frame(	const unsigned int index  );
		
		/*
		Title: print_dim
		Description:  This is a helper function for the user to be sure the correct dimensions 
			are obtained
		Inputs:		None
		Outputs:  	None:  It will be a console output
		*/
		void print_dim();
		
		
		// Getters and Setters for the private variables
		void Dim_x(unsigned int num);
		unsigned int Dim_x();
		
		void Dim_y(unsigned int num);
		unsigned int Dim_y();
		
		void Dim_z(unsigned int num);
		unsigned int Dim_z();
		
		void Filename(char *name);
		const char* Filename();
		
		// Getter for the data frame
		float *Data_Frame();
		
		// Getter for the data frame
		float *Data_All();
		
		// This is needed to get the frame and dimension data to the Contracting_Grid class.
		//friend void Contracting_Grid::getFrameData(Cine c);
};



#endif