CC=g++
OPTS= -Wall --pedantic -std=c++17
LIBS= -lmingw32 -lsfml-graphics -lsfml-window -lsfml-system -lsfml-main -mwindows
HDR_PATH = -IC:\SFML-2.6.1\include
LIB_PATH = -LC:\SFML-2.6.1\lib
SRC=SpatialSSA.cpp SpatialSSA_io.cpp SpatialSSA_model.cpp SpatialSSA_rendering.cpp

all: $(SRC:.cpp=.o)
	$(CC) $(OPTS) $(LIB_PATH) -o SpatialSSA $(SRC:.cpp=.o) $(LIBS)

%.o: %.cpp
	$(CC) $(OPTS) $(HDR_PATH) -c $< -o $@ 

# leaky program stays alive after shutdown | kill before fresh build
kill:
	TASKKILL /f /im SpatialSSA.exe 

# rm is a GNU command which is part of make environment
clear:
	rm *.o SpatialSSA.exe

# 
all_manual:
	g++ -IC:\SFML-2.6.1\include -c ${SRC}

link_manual:
	g++ -LC:\SFML-2.6.1\lib -o SpatialSSA.exe SpatialSSA.o SpatialSSA_io.o SpatialSSA_model.o SpatialSSA_rendering.o -lmingw32 -lsfml-graphics -lsfml-window -lsfml-system -lsfml-main -mwindows
