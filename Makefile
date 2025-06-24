CC := clang
CFLAGS := -g -Wall -fsanitize=address -mcpu=apple-m1

all : chaos_comparison

clean :
	rm -rf chaos_comparison chaos_comparison.dSYM

chaos_comparison : chaos_comparison.c phys_math.c phys_math.h butcher_tableau.h
	$(CC) $(CFLAGS) -o chaos_comparison chaos_comparison.c phys_math.c -lm -Iexternal/linmath.h

phys_math.o : phys_math.c phys_math.h butcher_tableau.h
	$(CC) $(CFLAGS) -c phys_math.c

.PHONY: all clean
