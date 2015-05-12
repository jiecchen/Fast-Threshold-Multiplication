# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# the build target executable:
TARGET = test

#LIBS = -larmadillo
RM = rm -f


all: test.o Algebra.o
	$(CC) $(CFLAGS) -o test *.o

$(TARGET).o: $(TARGET).cpp Algebra.h
	$(CC) $(CFLAGS) -c $(TARGET).cpp Algebra.h

Algebra.o: Algebra.h Algebra.cpp
	$(CC) $(CFLAGS) -c Algebra.cpp Algebra.h

clean:
	$(RM) $(TARGET) *~ *.o *.gch
