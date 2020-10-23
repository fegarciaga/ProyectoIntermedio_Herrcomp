CXX=g++
CXXFLAGS=-std=c++17
SANITFLAGS =-fsanitize=leak -fsanitize=undefined -fsanitize=address
OPFLAGS= -O2
GPROFFLAGS= -pg
VALGRINDFLAGS=--tool=memcheck --track-origins=yes --leak-check=full


SOURCES = main.cpp entropy.cpp
OBJ = $(SOURCES:.cpp=.o)
DEPS = entropy.h


all: main.x 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(SANITFLAGS)  $(OPFLAGS)  $^ -c $<
main.x: main.o entropy.o
	$(CXX) $(CXXFLAGS) $(SANITFLAGS) $(OPFLAGS) $^ -o $@

gprof: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OPFLAGS) $(GPROFFLAGS) $^ -c $<
	$(CXX) $(CXXFLAGS) $(OPFLAGS) $(GPROFFLAGS) $(OBJ)
	./a.out 1e6 400 1 20 10
	gprof ./a.out > gprof-report.txt

cachegrind: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OPFLAGS) $^ -c $<
	$(CXX) $(CXXFLAGS) $(OPFLAGS) $(OBJ)
	valgrind --tool=cachegrind ./a.out 1e6 400 1 20 10
	cg_annotate $$(ls -t cachegrind.out.* | head -n 1) > cachegrind-report.txt

memcheck: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(OPFLAGS) $^ -c $<
	$(CXX) $(CXXFLAGS) $(OPFLAGS) $(OBJ)
	valgrind $(VALGRINDFLAGS) ./a.out 1e6 400 1 20 10
clean:
	rm -f *.x *.o a.out *.out.* *.out
