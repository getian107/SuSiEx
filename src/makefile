CXXFLAGS=-O3 -Wall -Wextra
LDFLAGS=-O3 -flto

all: SuSiEx

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CXXFLAGS+=-Xpreprocessor -fopenmp
	LDLIBS+=-lomp
else
	CXXFLAGS=-fopenmp
endif

SuSiEx: main.cpp data.o model.o
	mkdir -p ../bin
	$(CXX) $(LDFLAGS) main.cpp data.o model.o -o ../bin/SuSiEx $(LDLIBS)

data.o: data.cpp data.hpp
	$(CXX) $(CXXFLAGS) -c data.cpp

model.o: model.cpp model.hpp
	$(CXX) $(CXXFLAGS) -c model.cpp

.PHONY: clean
clean:
	rm -rf *.o ../bin/SuSiEx
