all:

include ${SLEPC_DIR}/conf/slepc_common
#PCC=${CXX}
#PCC_LINKER=${CXX}


BASIC_OFILES = Geometry.o Tools.o Mode.o Moperator.o Pml.o Passive.o
NEWTON_OFILES = Newton.o Jacobian.o Creeper.o

NLOPT = /usr/include/nlopt.hpp
CHEBPACK = /Users/daveliu/chebpack
#CFLAGS = -I${NLOPT} -I${CHEBPACK} -L${CHEBPACK} -lnlopt -lchebpack -lfftw3
CFLAGS = -lnlopt
CLEANFILES = PassiveOut SaltOut log* #.o files removed by default!
COMMAND = ${CLINKER} $^ -o $@Out ${SLEPC_LIB}

		
SaltShared: Salt.o ${BASIC_OFILES} ${NEWTON_OFILES}
	${CLINKER} $^ -shared -o saltlib.dylib ${SLEPC_LIB}

main: main.o
	${CLINKER} $^ -o SaltOut saltlib.dylib ${SLEPC_LIB}

Test: Test.o
	${CLINKER} $^ -o TestOut saltlib.dylib ${SLEPC_LIB}


