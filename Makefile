all:

include ${SLEPC_DIR}/conf/slepc_common


BASIC_OFILES = Geometry.o Tools.o Mode.o Moperator.o Pml.o Passive.o
NEWTON_OFILES = Newton.o Jacobian.o Creeper.o

UNAME = $(shell uname)
ifeq (${UNAME},Darwin)
	EXTENSION = .dylib
else
	EXTENSION = .so
endif

SALTLIB = ${PWD}/saltlib${EXTENSION}
SALTOUT = SaltOut
TESTOUT = TestOut
CLEANFILES = *.dylib *.so *Out

		
Salt: Salt.o ${BASIC_OFILES} ${NEWTON_OFILES}
	${CLINKER} $^ -shared -o ${SALTLIB} ${SLEPC_LIB}
	rm *.o;

Main: main.o
	${CLINKER} $^ -o ${SALTOUT} ${SALTLIB} ${SLEPC_LIB}
	echo ${UNAME}

Test: Test.o
	${CLINKER} $^ -o ${TESTOUT} ${SALTLIB} ${SLEPC_LIB}

cleanout:
	rm ${SALTLIB} ${SALTOUT} ${TESTOUT}