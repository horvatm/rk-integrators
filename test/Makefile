CC=g++
CXXFLAGS=-O3 -Wall

all: test_check_order  test_integrator  test_method

%:%.cpp *.h
	$(CC) $(CXXFLAGS) $*.cpp -o $*
