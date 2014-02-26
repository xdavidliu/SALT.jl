#include <slepc.h>


void printname(char** names, int N){
	SlepcInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL); 

	int i;
	for(i =0; i<N; i++)
	PetscPrintf(PETSC_COMM_WORLD, "Name %i = %s\n", i, names[i]);
	
}


/*

names = ["David", "Dustin", "Noam"];
ccall((:printname, "/Users/daveliu/Documents/saltc/foonames2"), Void, (Ptr{Ptr{Uint8}}, Cint), names, length(names));

*/