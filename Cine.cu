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
	

// #define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string>
#include "Cine.h"

// Macro used to terminate the program with an error message
#define ABORT(__message) do {											\
  perror(__message);												\
  fprintf(stderr, "File: \"%s\", function: \"%s\", line: %d\n", __FILE__, __FUNCTION__, __LINE__);		\
  exit(EXIT_FAILURE);												\
} while(0)


// Macro used to terminate the program with an error message
#define ERROR(__message) do {											\
  fprintf(stderr, "%s\n", __message);										\
  fprintf(stderr, "File: \"%s\", function: \"%s\", line: %d\n", __FILE__, __FUNCTION__, __LINE__);		\
  exit(EXIT_FAILURE);												\
} while(0)


#define MAX(__a, __b)					(((__a) > (__b)) ? (__a) : (__b))
#define MIN(__a, __b)					(((__a) < (__b)) ? (__a) : (__b))
#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))
#define MAP_3D(__dimx, __dimy, __dimz, __x, __y, __z)	(((__z) * (__dimy) + (__y)) * (__dimx) + (__x))

	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct time64_t {
  uint32_t fractions;
  uint32_t seconds;
} __attribute__((__packed__));


struct cine_header_t {
  uint16_t type;
  uint16_t header_size;
  uint16_t compression;
  uint16_t version;
  int32_t first_image;
  uint32_t image_count;
  int32_t first_saved;
  uint32_t saved_count;
  uint32_t offset_header;
  uint32_t offset_setup;
  uint32_t offset_offsets;
  struct time64_t start_time;
} __attribute__((__packed__));


struct bitmap_header_t {
  uint32_t header_size;
  int32_t frame_width;
  int32_t frame_height;
  uint16_t num_planes;
  uint16_t bit_count;
  uint32_t compression;
  uint32_t frame_size_bytes;
  int32_t pixel_per_meter_x;
  int32_t pixel_per_meter_y;
  uint32_t num_color_indexed;
  uint32_t num_color_required;
} __attribute__((__packed__));


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Default Constructor referred to in Cine.h
	Cine::Cine() {
	
		dim_x = 0;
		dim_y = 0;
		dim_z = 0;
		filename = "";
	}
	
	// Overloaded constructor of Cine class
	Cine::Cine( const char *name ){
		
		dim_x = 0;
		dim_y = 0;
		dim_z = 0;
		
		// Make deep copy of filename
		filename = new char[strlen(name)];		
		std::string temp = std::string(name);
		filename = temp.c_str();
		
		// this will populate the dim_x, dim_y, and dim_z fields
		read_cine_info();		
		data_all = new float[dim_x * dim_y * dim_z];
		data_frame = new float[dim_x * dim_y];
		
	}
	
	// details referred to in Cine.h
	void Cine::write_dat_3d(	const char *filename_out	) {
							
		char err_msg[256];
		int fd;

		printf("Size: %d %d %d\n", dim_x, dim_y, dim_z);
		
		const unsigned int x = dim_x;
		const unsigned int y = dim_y;
		const unsigned int z = dim_z;
		
		printf("Size: %d %d %d\n", dim_x, dim_y, dim_z);
		
		// Create the file
		fd = creat(filename_out, 0644);
		if(fd == -1) {
		snprintf(err_msg, 256, "Cannot open file \"%s\"", filename_out);
		ERROR(err_msg);
		}
		
		printf("Size: %d %d %d\n", dim_x, dim_y, dim_z);
		// Write the size along the X dimension
		write(fd, & x, sizeof(unsigned int));
		// Write the size along the Y dimension
		write(fd, & y, sizeof(unsigned int));
		// Write the size along the Z dimension
		write(fd, & z, sizeof(unsigned int));
		// Write the data
		write(fd, data_all, dim_x * dim_y * dim_z * sizeof(float));
		// Close the file
		close(fd);
		return;
}
	
	void Cine::write_eps(	const char *filename_out	) {
								
		register unsigned int i, j, tmp;
		char err_msg[256];
		float min, max;
		time_t now;
		FILE *fid;

		now = time(NULL);
		// Find the minimum and maximum of the data so
		// the data will use the whole [0, 255] range
		min = max = data_frame[0];
		for(i = 0; i < (dim_x * dim_y); ++i) {
		if(min > data_frame[i]) {
		  min = data_frame[i];
		}
		if(max < data_frame[i]) {
		  max = data_frame[i];
		}
		}
		// Open the file for writing
		fid = fopen(filename_out, "w");
		if(fid == NULL) {
		snprintf(err_msg, 256, "Cannot open file \"%s\"", filename_out);
		ABORT(err_msg);
		}
		// Write the header
		fprintf(fid, "%%!PS-Adobe-3.0 EPSF-3.0\n");
		fprintf(fid, "%%%%Title: %s\n", filename_out);
		fprintf(fid, "%%%%Creator: Luca Caucci, caucci@email.arizona.edu\n");
		fprintf(fid, "%%%%CreationDate: %s", ctime(& now));
		fprintf(fid, "%%%%BoundingBox: 0 0 %d %d\n", dim_x, dim_y);
		fprintf(fid, "%%%%LanguageLevel: 1\n");
		fprintf(fid, "%%%%Pages: 0\n");
		fprintf(fid, "%%%%DocumentData: Clean7Bit\n");
		fprintf(fid, "%d %d scale\n", dim_x, dim_y);
		fprintf(fid, "%d %d 8 [%d 0 0 -%d 0 %d]\n", dim_x, dim_y, dim_x, dim_y, dim_y);
		fprintf(fid, "{currentfile %d string readhexstring pop} bind\n", dim_x);
		fprintf(fid, "image\n");
		// Check if the image is constant
		if(min != max) {
		for(i = 0; i < dim_y; ++i) {
		  for(j = 0; j < dim_x; ++j) {
			// Scale the datum and convert it to an integer number
			tmp = (int) floor(256.00 * (((float) (data_frame[MAP_2D(dim_x, dim_y, j, i)] - min)) / (max - min)));
			// Write the datum to the file as a hexadecimal number
			fprintf(fid, "%02x", MIN(tmp, 255));
		  }
		  fprintf(fid, "\n");
		}
		} else {
		// If the data are constant, the image is filled with a constant value
		for(i = 0; i < dim_y; ++i) {
		  for(j = 0; j < dim_x; ++j) {
			fprintf(fid, "%02x", 128);
		  }
		  fprintf(fid, "\n");
		}
		}
		fprintf(fid, "%%%%EOF\n");
		fclose(fid);
		return;
}
	
	// read_cine_3d referred to in Cine.h
	void Cine::read_cine_3d(	) {
		
		struct bitmap_header_t bitmap_header;
		struct cine_header_t cine_header;
		uint32_t frame_size_bytes;
		uint32_t annot_size;
		char err_msg[256];
		int64_t *offsets;
		uint16_t *frame;
		int i, j, k;
		int fd;
		
		fd = open(filename, O_RDONLY);
		if(fd == -1) {
		snprintf(err_msg, 256, "Cannot open file \"%s\"", filename);
		ERROR(err_msg);
		}
		
		read(fd, & cine_header, sizeof(cine_header));
		if(cine_header.type != ('C' | ('I' << 8))) {
		ERROR("Invalid CINE file!");
		}
		if(cine_header.header_size != sizeof(cine_header)) {
		ERROR("Unexpected CINE header length!");
		}
		if(cine_header.compression != 0) {
		ERROR("Unsupported CINE compression method!");
		}
		if(cine_header.version != 1) {
		ERROR("Unsupported CINE version!");
		}
		if(dim_z != cine_header.saved_count) {
		ERROR("Unexpected data size in the Z dimension!");
		}
		
		read(fd, & bitmap_header, sizeof(bitmap_header));
		if(bitmap_header.header_size != sizeof(bitmap_header)) {
		ERROR("Unexpected bitmap header length!");
		}
		if(dim_x != bitmap_header.frame_width) {
		ERROR("Unexpected data size in the X dimension!");
		}
		if(dim_y != bitmap_header.frame_height) {
		ERROR("Unexpected data size in the Y dimension!");
		}
		if(bitmap_header.num_planes != 1) {
		ERROR("Unsupported number of planes!");
		}
		if(bitmap_header.bit_count != 16) {
		ERROR("Unsupported number of bits per pixel!");
		}
		if(bitmap_header.frame_size_bytes != (dim_x * dim_y * sizeof(*frame))) {
		ERROR("Unexpected byte size of image frames!");
		}
		
		lseek(fd, cine_header.offset_offsets, SEEK_SET);
		offsets = (int64_t *)calloc(dim_z, sizeof(*offsets));
		read(fd, offsets, dim_z * sizeof(*offsets));
		frame = (uint16_t *)calloc(dim_x * dim_y, sizeof(*frame));
		
		//data_all = new float[dim_x * dim_y * dim_z];
		
		for(k = 0; k < dim_z; ++k) {
			lseek64(fd, offsets[k], SEEK_SET);
			read(fd, & annot_size, sizeof(annot_size));
			if(annot_size > 8) {
			  lseek(fd, annot_size - 8, SEEK_CUR);
			}
			
			read(fd, & frame_size_bytes, sizeof(frame_size_bytes));
			if(frame_size_bytes != (dim_x * dim_y * sizeof(*frame))) {
			  ERROR("Unexpected byte size of image frames!");
			}
			read(fd, frame, frame_size_bytes);
			
			for(i = 0; i < dim_x; ++i) {
			  for(j = 0; j < dim_y; ++j) {
				data_all[MAP_3D(dim_x, dim_y, dim_z, i, j, k)] = (float) frame[MAP_2D(dim_x, dim_y, i, dim_y - j - 1)];
			  }
			}
		}
		
		free(frame);
		free(offsets);
		close(fd);
		return;
	}
	
	// read_cine_info referred to in Cine.h 
	void Cine::read_cine_info(	) {
		
		struct bitmap_header_t bitmap_header;
		struct cine_header_t cine_header;
		char err_msg[256];
		int fd;
		  
		fd = open(filename, O_RDONLY);
		if(fd == -1) {
			snprintf(err_msg, 256, "Cannot open file \"%s\"", filename);
			ERROR(err_msg);
		}
		read(fd, & cine_header, sizeof(cine_header));
		if(cine_header.type != ('C' | ('I' << 8))) {
			ERROR("Invalid CINE file!");
		}
		if(cine_header.header_size != sizeof(cine_header)) {
			ERROR("Unexpected CINE header length!");
		}
		if(cine_header.compression != 0) {
			ERROR("Unsupported CINE compression method!");
		}
		if(cine_header.version != 1) {
			ERROR("Unsupported CINE version!");
		}
		read(fd, & bitmap_header, sizeof(bitmap_header));
		if(bitmap_header.header_size != sizeof(bitmap_header)) {
			ERROR("Unexpected bitmap header length!");
		}
		if(bitmap_header.num_planes != 1) {
			ERROR("Unsupported number of planes!");
		}
		if(bitmap_header.bit_count != 16) {
			ERROR("Unsupported number of bits per pixel!");
		}
		dim_x = bitmap_header.frame_width;
		dim_y = bitmap_header.frame_height;
		dim_z = cine_header.saved_count;
		if(bitmap_header.frame_size_bytes != (dim_x * dim_y * sizeof(uint16_t))) {
			ERROR("Unexpected byte size of image frames!");
		}
		close(fd);
		return;
	}
	
	// read_cine_frame referred to in Cine.h
	void Cine::read_cine_frame(	const unsigned int index) {
		
	  struct bitmap_header_t bitmap_header;
	  struct cine_header_t cine_header;
	  uint32_t frame_size_bytes;
	  uint32_t annot_size;
	  char err_msg[256];
	  uint16_t *frame;
	  int64_t offset;
	  int i, j;
	  int fd;
	  
		
		
	  fd = open(filename, O_RDONLY);
	  if(fd == -1) {
		snprintf(err_msg, 256, "Cannot open file \"%s\"", filename);
		ERROR(err_msg);
	  }
	  
	  printf("Size2: %d %d %d\n", dim_x, dim_y, dim_z);
	  
	  read(fd, & cine_header, sizeof(cine_header));
	  if(cine_header.type != ('C' | ('I' << 8))) {
		ERROR("Invalid CINE file!");
	  }
	  if(cine_header.header_size != sizeof(cine_header)) {
		ERROR("Unexpected CINE header length!");
	  }
	  if(cine_header.compression != 0) {
		ERROR("Unsupported CINE compression method!");
	  }
	  if(cine_header.version != 1) {
		ERROR("Unsupported CINE version!");
	  }
	  if(index >= cine_header.saved_count) {
		ERROR("Invalid frame index!");
	  }
	  
	  printf("Size3: %d %d %d\n", dim_x, dim_y, dim_z);
	  
	  read(fd, & bitmap_header, sizeof(bitmap_header));
	  if(bitmap_header.header_size != sizeof(bitmap_header)) {
		ERROR("Unexpected bitmap header length!");
	  }
	  if(dim_x != bitmap_header.frame_width) {
		ERROR("Unexpected data size in the X dimension!");
	  }
	  if(dim_y != bitmap_header.frame_height) {
		ERROR("Unexpected data size in the Y dimension!");
	  }
	  if(bitmap_header.num_planes != 1) {
		ERROR("Unsupported number of planes!");
	  }
	  if(bitmap_header.bit_count != 16) {
		ERROR("Unsupported number of bits per pixel!");
	  }
	  if(bitmap_header.frame_size_bytes != (dim_x * dim_y * sizeof(*frame))) {
		ERROR("Unexpected byte size of image frames!");
	  }
	  printf("Size4: %d %d %d\n", dim_x, dim_y, dim_z);
	  lseek(fd, cine_header.offset_offsets + sizeof(offset) * index, SEEK_SET);
	  printf("Size5: %d %d %d\n", dim_x, dim_y, dim_z);
	  read(fd, & offset, sizeof(offset));
	  
	  frame = (uint16_t *)calloc(dim_x * dim_y, sizeof(*frame));
	  
	  lseek64(fd, offset, SEEK_SET);
	  printf("Size6: %d %d %d\n", dim_x, dim_y, dim_z);
	  read(fd, & annot_size, sizeof(annot_size));
	  if(annot_size > 8) {
		lseek(fd, annot_size - 8, SEEK_CUR);
	  }
	  printf("Size7: %d %d %d\n", dim_x, dim_y, dim_z);
	  read(fd, & frame_size_bytes, sizeof(frame_size_bytes));
	  if(frame_size_bytes != (dim_x * dim_y * sizeof(*frame))) {
		ERROR("Unexpected byte size of image frames!");
	  }
	  printf("Size8: %d %d %d\n", dim_x, dim_y, dim_z);
	  read(fd, frame, frame_size_bytes);
	  
	  for(i = 0; i < dim_x; ++i) {
		for(j = 0; j < dim_y; ++j) {
		  data_frame[MAP_2D(dim_x, dim_y, i, j)] = (float) frame[MAP_2D(dim_x, dim_y, i, dim_y - j - 1)];
		}
	  }
	  printf("Size9: %d %d %d\n", dim_x, dim_y, dim_z);
	  free(frame);
	  close(fd);
	  return;
	}
	
	void Cine::print_dim(){
		printf("Size: %d %d %d\n", dim_x, dim_y, dim_z);
	}

	// These are inline getters and setters reffered to in Cine.h.
	// The setters for the dim_x, dim_y, and dim_z variables may not be necessary or desirable.
	void Cine::Dim_x(unsigned int num){dim_x = num;}
		
	unsigned int Cine::Dim_x(){return dim_x;}
		
	void Cine::Dim_y(unsigned int num){dim_y = num;}
		
	unsigned int Cine::Dim_y(){return dim_y;}
		
	void Cine::Dim_z(unsigned int num){dim_z = num;}
		
	unsigned int Cine::Dim_z(){return dim_z;}
		
	void Cine::Filename(char *name){
		filename = name; 
		read_cine_info();
	}
		
	const char * Cine::Filename(){return filename;}
