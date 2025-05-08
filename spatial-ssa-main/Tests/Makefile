CC=g++
OPTS= -Wall --pedantic -std=c++17
SRC=random_numbers.cpp

start: $(SRC:.cpp=.o)
	$(CC) $(OPTS) -o test $(SRC:.cpp=.o)

%.o: %.cpp
	$(CC) $(OPTS) -c $< -o $@

clean:
	rm -f *.o
	rm -f test.exe
	rm -f test

