# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -std=c++0x

# the build target executable:
TARGET = test

#LIBS = -larmadillo
RM = rm -f


all: test.o Algebra.o Sketch.o
	$(CC) $(CFLAGS) -o test *.o

$(TARGET).o: $(TARGET).cpp Algebra.h Sketch.h
	$(CC) $(CFLAGS) -c $(TARGET).cpp *.h

Algebra.o: Algebra.h Algebra.cpp
	$(CC) $(CFLAGS) -c Algebra.cpp Algebra.h


Sketch.o: Sketch.h Algebra.h Sketch.cpp
	$(CC) $(CFLAGS) -c Sketch.h Algebra.h Sketch.cpp

clean:
	$(RM) $(TARGET) *~ *.o *.gch
