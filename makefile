all: main.x

%.o: %.cpp
	g++ -c $<
main.x: main.o entropy.o
	g++ $^ -o $@

clean:
	rm -f *.x *.o a.out
