all: mtwistlib compile_main
	echo "Everything is done!"
compile_main:
	g++ -o main main.cpp library/mtwist-1.5/mtwist.o library/mtwist-1.5/randistrs.o -lm -Ofast -Wall

mtwistlib: mtwist randistrs
	echo "Library Compiled" 

mtwist:
	gcc -c -O3 library/mtwist-1.5/mtwist.c

randistrs:
	gcc -c -O3 library/mtwist-1.5/randistrs.c

