#include <iostream>
#include "Contracting_Grid.h"
#include "Cine.h"
#include "Device.h"

#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x))
#define MAP_3D(__dimx, __dimy, __dimz, __x, __y, __z)	(((__z) * (__dimy) + (__y)) * (__dimx) + (__x))

using namespace std;

int main(){
	
	Cine cine("../../C/test.cine");
	
	// Shows that the class is working
	//cine.print_dim();
	
	cine.read_cine_3d();
	
	//cine.write_dat_3d("testdata.dat");
	
	unsigned int dim_x = cine.Dim_x();
	unsigned int dim_y = cine.Dim_y();
	unsigned int dim_z = cine.Dim_z();
	
	//float *all_frames = new float[cine.Dim_x() * cine.Dim_y() * cine.Dim_z()]; 
	//all_frames = cine.Data_All();
	
	Contracting_Grid cg(cine);
	
	//for(unsigned int i = 0; i < dim_z ; i++){
		cg.Frame(cine.read_cine_frame(0));
		cg.find_4D_Center();
		cg.find_4D_Center_Cuda();
		//cg.findSharpCenter_Rec();
		
		//cg.findBroadCenter_brute();
		//cg.findDirection();
	//}
	
	cine.print_dim();
	
	return 0;
}
