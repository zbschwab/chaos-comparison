CC := clang
CPPFLAGS := -I/opt/homebrew/opt/openblas/include
CFLAGS := -g -Wall -fsanitize=address -mcpu=apple-m1
LDFLAGS := -L/opt/homebrew/opt/openblas/lib -lopenblas -lm 

all : chaos_comparison

clean :
	rm -rf chaos_comparison chaos_comparison.dSYM

chaos_comparison : chaos_comparison.c phys_math.c phys_math.h butcher_tableau.h
	$(CC) $(CFLAGS) $(CPPFLAGS) chaos_comparison.c phys_math.c $(LDFLAGS) -o chaos_comparison

phys_math.o : phys_math.c phys_math.h butcher_tableau.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c phys_math.c 

.PHONY: all clean
