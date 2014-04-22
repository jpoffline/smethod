OBJECTS = main.o initialconditions.o kgsolve.o setup.o write.o 
INC = -Wall -I/usr/local/include
LIBS = -L/usr/local/lib -lm  -L/usr/local/boost 
GSL_LIBS = -lgsl -lgslcblas
BOOST_LIBS = -lboost_filesystem -lboost_system -lboost_timer
CC = g++

all: main

%.o:    src/%.cpp
	$(CC) -c $(INC) $< -o $@

main: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) $(GSL_LIBS) $(BOOST_LIBS) -o $@

clean:
	-rm -f $(OBJECTS) main

