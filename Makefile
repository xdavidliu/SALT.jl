all:

include ${SLEPC_DIR}/conf/slepc_common
PCC=${CXX}
PCC_LINKER=${CXX}


BASIC_OFILES = Geometry.o Tools.o Mode.o Moperator.o Pml.o
NEWTON_OFILES = Newton.o Jacobian.o

NLOPT = /usr/include/nlopt.hpp
CHEBPACK = /Users/daveliu/chebpack
#CFLAGS = -I${NLOPT} -I${CHEBPACK} -L${CHEBPACK} -lnlopt -lchebpack -lfftw3
CFLAGS = -lnlopt
CLEANFILES = PassiveOut SaltOut log* #.o files removed by default!
COMMAND = ${CLINKER} $^ -o $@Out ${SLEPC_LIB}


Passive: Passive.o ${BASIC_OFILES}
	${COMMAND}

Salt: Salt.o ${BASIC_OFILES} ${NEWTON_OFILES}
	${COMMAND}

Opt: Opt.o ${BASIC_OFILES} ${NEWTON_OFILES}
	${COMMAND}


Test: Test.o
	${COMMAND}


