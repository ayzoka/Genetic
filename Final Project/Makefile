QKDOpt: Algebra.o Complex.o quantum.o attack.o main.o
	g++ -O3 -o QKDOpt -framework accelerate Algebra.o Complex.o  quantum.o attack.o main.o

main.o: main.cpp GA.h
	g++ -O3 -c main.cpp

Algebra.o: lib-quantum/Algebra.h lib-quantum/Algebra.cpp
	g++ -O3 -c lib-quantum/Algebra.cpp -I /usr/include

Complex.o: lib-quantum/Complex.hpp lib-quantum/Complex.cpp
	g++ -O3 -c lib-quantum/Complex.cpp

quantum.o: lib-quantum/quantum.hpp lib-quantum/quantum.cpp
	g++ -O3 -c lib-quantum/quantum.cpp

attack.o: attack.h attack.cpp
	g++ -O3 -c attack.cpp

