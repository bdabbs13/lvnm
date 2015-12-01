OBJECTS = mmsbm.o lsm.o
INCLUDES = -I/usr/share/R/include
FLAGS = -DNDEBUG -fpic -O2 -g
LIBS = -L/usr/lib/R/lib -lR
all: mmsbm.so

mmsbm.so: $(OBJECTS)
	g++ -shared -o mmsbm.so $(OBJECTS) $(LIBS)

mmsbm.o: mmsbm.cpp mmsbm.h
	g++ $(INCLUDES) $(FLAGS) -c mmsbm.cpp -o mmsbm.o

lsm.o: lsm.cpp lsm.h
	g++ $(INCLUDES) $(FLAGS) -c lsm.cpp -o lsm.o

clean:
	-rm $(OBJECTS) mmsbm.so