main: main.c functions.c functions.h
	gcc main.c functions.c -lgsl -lgslcblas -lm -W -Wall -Wextra -O3 -o main
