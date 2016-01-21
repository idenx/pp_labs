sources := main.cpp
headers := matrix.hpp eigenvalue_searcher.hpp matrix.cpp timer.hpp main.hpp
objects := $(patsubst %.cpp, %.o, $(sources))

CFLAGS := -Wall -Werror -ggdb3 -O3 -fopenmp
LDFLAGS := $(CFLAGS)

%.o: %.cpp $(headers) Makefile
	$(CXX) $(CFLAGS) -c $<

all: $(objects)
	$(CXX) $(LDFLAGS) $^ -o ./main

clean:
	@rm *.o
	@rm main
