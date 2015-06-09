# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -O3 -Wall -std=c++0x 

# the build target executable:
TARGET = test createData test_dblp

#LIBS = -larmadillo
RM = rm -f

test: test.o Algebra.o Sketch.o
	$(CC) $(CFLAGS) -o test test.o Algebra.o Sketch.o
test_dblp: test_dblp.o Algebra.o Sketch.o
	$(CC) $(CFLAGS) -o test_dblp test_dblp.o Algebra.o Sketch.o	

test_dblp.o: test_dblp.cpp Algebra.h Sketch.h utils.h
	$(CC) $(CFLAGS)  -c test_dblp.cpp *.h

createData: createData.o Algebra.o
	$(CC) $(CFLAGS) -o createData createData.o Algebra.o

createData.o: createData.cpp Algebra.h
	$(CC) $(CFLAGS)  -c createData.cpp *.h

test.o: test.cpp Algebra.h Sketch.h utils.h
	$(CC) $(CFLAGS)  -c test.cpp *.h

Algebra.o: Algebra.h Algebra.cpp
	$(CC) $(CFLAGS)  -c Algebra.cpp Algebra.h


Sketch.o: Sketch.h Algebra.h utils.h Sketch.cpp
	$(CC) $(CFLAGS) -c *.h Sketch.cpp

clean:
	$(RM) $(TARGET) *~ *.o *.gch
