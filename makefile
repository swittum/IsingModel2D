GCC = g++
CFLAGS = -std=c++98 -Wall -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib
LDLIBS = -lhdf5

all: ising

ising: ising.cpp
	$(GCC) $(CFLAGS) $(LDFLAGS) -o ising ising.cpp $(LDLIBS)

run: ising
	./ising
	python plotting.py
	open ./ani.mp4

clean:
	rm -rf gas
