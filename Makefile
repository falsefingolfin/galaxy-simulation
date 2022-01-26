CXX=gcc
CXXFLAGS?=-o -lm


aarseth: aarseth.c
	gcc -O3 -o aarseth aarseth.c -lm

clean:
	rm -f aarseth *.o