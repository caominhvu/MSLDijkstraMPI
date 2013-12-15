UNAME := $(shell uname -s)

ifeq ($(UNAME), Darwin)
	CXX = openmpicxx
endif
ifeq ($(UNAME), Linux)
	CXX = mpic++
endif


CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		src/PrimAlgo/PrimAlgo.o \
			src/DijkstraAlgo/DijkstraAlgo.o \
			src/MPIMain.o

LIBS =

TARGET =	MPIMain

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
