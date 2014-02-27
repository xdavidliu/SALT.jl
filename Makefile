all:

include ${SLEPC_DIR}/conf/slepc_common


BASIC_OFILES = Geometry.o Tools.o Mode.o Moperator.o Pml.o Passive.o
NEWTON_OFILES = Newton.o Jacobian.o Creeper.o

SALTLIB = saltlib.dylib
SALTOUT = SaltOut
TESTOUT = TestOut
#CLEANFILES = 

		
Salt: Salt.o ${BASIC_OFILES} ${NEWTON_OFILES}
	${CLINKER} $^ -shared -o ${SALTLIB} ${SLEPC_LIB}

Main: main.o
	${CLINKER} $^ -o ${SALTOUT} ${SALTLIB} ${SLEPC_LIB}

Test: Test.o
	${CLINKER} $^ -o ${TESTOUT} ${SALTLIB} ${SLEPC_LIB}

cleanout:
	rm ${SALTLIB} ${SALTOUT} ${TESTOUT}