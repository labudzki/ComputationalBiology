CC=g++
OPTS= -Wall --pedantic -std=c++17
SRC=SpatialSSA_setup_1.cpp SpatialSSA_io.cpp SpatialSSA_model.cpp
LIBS=-lm

all: $(SRC:.cpp=.o)
	$(CC) $(OPTS) -o SpatialSSA_setup_1 $(SRC:.cpp=.o) $(LIBS)

%.o: %.cpp
	$(CC) $(OPTS) -c $< -o $@

clean:
	rm -f *.o
	rm -f SpatialSSA
