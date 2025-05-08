CC=g++
OPTS= -Wall --pedantic -std=c++17
LIBS= -lsfml-graphics -lsfml-window -lsfml-system
HDR_PATH = -I...\Libraries\SFML-2.5.1-64\include
LIB_PATH = -L...\Libraries\SFML-2.5.1-64\lib
SRC=Density.cpp

all: $(SRC:.cpp=.o)
	$(CC) $(OPTS) -o Density $(SRC:.cpp=.o) $(LIBS) $(LIB_PATH)

%.o: %.cpp
	$(CC) $(OPTS) -c $< -o $@ $(HDR_PATH)

clean:
	rm -f *.o
	rm -f Density.exe
