all: SuSiEx

SuSiEx: main.cpp data.o model.o
	mkdir -p ../bin
	g++ main.cpp data.o model.o -O3 -o ../bin/SuSiEx -fopenmp

data.o: data.cpp data.hpp
	g++ -c data.cpp -fopenmp

model.o: model.cpp model.hpp
	g++ -c model.cpp -fopenmp

clean:
	rm -rf *.o ../bin/SuSiEx
