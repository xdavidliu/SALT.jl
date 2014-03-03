#include <slepc.h>

int main(int argc, char** argv){

	SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

	SlepcFinalize();
}
