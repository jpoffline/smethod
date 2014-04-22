OBJECTS = FD_FFT_3D.o
INC = -Wall -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas 
LIBS_FFT = -lfftw3 -lm
LIBS_BOOST = -L/usr/local/boost -lboost_filesystem -lboost_system -lboost_timer
CC = g++

all: main

%.o:     %.cpp
	$(CC) -c $(INC) $< -o $@

main: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) $(LIBS_FFT) $(LIBS_BOOST) -o $@

clean:
	-rm -f $(OBJECTS) main

