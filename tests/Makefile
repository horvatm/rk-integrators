CXX = g++

DEPS = $(wildcard include/*.h)
INCLUDE = -I. -I../include   
LDFLAGS = -L. 
CXXFLAGS = -O3 -Wall -fopenmp -std=c++17

.PHONY: all clean

all: example_C.e  example_R.e test_check_order.e  test_integrator.e  test_method.e

%.e : %.cpp ${DEPS}
	${CXX} ${CXXFLAGS} ${INCLUDE} $< -o $@ ${LDFLAGS}

clean:
	rm -rf *.e  
