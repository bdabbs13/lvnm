CXX=g++
CXXFLAGS= -DNDEBUG -fpic -O2 -g
OBJECTS = lvnm.o sbm.o mmsbm.o lsm.o helper.o
LVNMINC = ../include
INCLUDES = -I/Library/Frameworks/R.framework/Resources/include -I$(LVNMINC)
LIBS = -L/Library/Frameworks/R.framework/Resources/lib -lR
VPATH = $(LVNMINC)

all: lvnm.so

lvnm.so: $(OBJECTS)
	$(CXX) -shared -o lvnm.so $(OBJECTS) $(LIBS)

lvnm.o: lvnm.cpp lvnm.h 
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c lvnm.cpp -o lvnm.o

sbm.o: sbm.cpp sbm.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c sbm.cpp -o sbm.o

mmsbm.o: mmsbm.cpp mmsbm.h helper.o
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c mmsbm.cpp -o mmsbm.o

lsm.o: lsm.cpp lsm.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c lsm.cpp -o lsm.o

helper.o: helper.cpp helper.h
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c helper.cpp -o helper.o



clean:
	-rm $(OBJECTS) lvnm.so
	-rm *~
	-rm *#