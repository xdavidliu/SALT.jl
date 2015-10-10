all: Salt Shared
	rm *.o

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
CLEANFILES = *.dylib *.so *Out

	
Shared: ${BASIC_OFILES} ${NEWTON_OFILES}
	${CLINKER} $^ -shared -o ${SALTLIB} ${SLEPC_LIB}

Salt: main.o ${BASIC_OFILES} ${NEWTON_OFILES}
	${CLINKER} $^ -o ${SALTOUT} ${SLEPC_LIB}

clean_data:
	rm *_file.m
