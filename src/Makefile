CC=g++
#CFLAGS=-g -Wall
INCLUDES=-I /home/nils/Development/seqan-src/include/
LIBS=-lgmp -lgmpxx

TARGET=ML-2015-82599-posterior

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(TARGET).cpp $(LIBS)

clean:
	$(RM) $(TARGET)
