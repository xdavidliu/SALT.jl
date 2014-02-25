#include <slepc.h>


class thing{

	public:
	Vec x;
	thing(){};
	thing(int N){ 
		PetscPrintf(PETSC_COMM_WORLD, "x before create = %i\n", x);
		VecCreateSeq(PETSC_COMM_SELF, N, &x);
		PetscPrintf(PETSC_COMM_WORLD, "x after create = %i\n", x);
	}
	
	~thing(){
		PetscPrintf(PETSC_COMM_WORLD, "x before destroy = %i\n", x);	
		VecDestroy(&x);
		PetscPrintf(PETSC_COMM_WORLD, "x after destroy = %i\n", x);			
	}
	
};


int main(int argc, char** argv){

	SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

	thing t;
	PetscPrintf(PETSC_COMM_WORLD, "x before everything = %i\n", t.x);
	t = thing(2);
	PetscPrintf(PETSC_COMM_WORLD, "x after everything = %i\n", t.x);
	
//	int N;
//	VecGetSize(t.x, &N);

	SlepcFinalize();
}
