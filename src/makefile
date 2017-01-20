CXX=g++
CXXFLAGS= -DNDEBUG -fpic -O2 -g
OBJECTS = lvnm.o sbm.o wsbm.o dynamic-sbm.o mmsbm.o lsm.o helper.o
LVNMINC = ../include
INCLUDES = -I/Library/Frameworks/R.framework/Resources/include -I$(LVNMINC)
LIBS = -L/Library/Frameworks/R.framework/Resources/lib -lR
VPATH = $(LVNMINC)

all: lvnm.so

lvnm.so: $(OBJECTS)
	$(CXX) -shared -o lvnm.so $(OBJECTS) $(LIBS)

lvnm.o: lvnm.cpp lvnm.h helper.o
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c lvnm.cpp -o lvnm.o

sbm.o: sbm.cpp sbm.h helper.o
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c sbm.cpp -o sbm.o

wsbm.o: wsbm.cpp wsbm.h helper.o
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c wsbm.cpp -o wsbm.o

dynamic-sbm.o: dynamic-sbm.cpp dynamic-sbm.h helper.o wsbm.o
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c dynamic-sbm.cpp -o dynamic-sbm.o

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