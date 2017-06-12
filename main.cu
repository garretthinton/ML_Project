#include <iostream>
#include "Cine.h"

using namespace std;

int main(){
	
	// Cine cine("test.cine");
	Cine cine("fourth.cine");

	// Shows that the class is working
	cine.print_dim();
	
	cine.read_cine_3d();
	
	cine.print_dim();
	
	cine.write_dat_3d("cine_fourth.dat");
	
	
	//for(int i = 0; i<cine.dim_z;i++){
	//	cine.read_cine_frame();
	//}
	
	
	return 0;
}
