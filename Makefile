CXX=g++-7
CXFLAGS=-W -Wall -O2 -std=c++11
LDFLAGS=-fopenmp
EXEC=simu

SRC=$(wildcard *.cpp)
HEADER=$(wildcard *.hpp)
OBJ=$(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXFLAGS) -o $@ -c $< $(LDFLAGS)

.PHONY: clean cleanall

clean:
	rm -rf $(OBJ)

cleanall: clean
	rm -rf $(EXEC)
