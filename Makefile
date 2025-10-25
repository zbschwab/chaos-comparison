CC := clang
CPPFLAGS := -I/opt/homebrew/opt/openblas/include
CFLAGS := -g -Wall -Wno-unused-variable -fsanitize=address -mcpu=apple-m1
LDFLAGS := -L/opt/homebrew/opt/openblas/lib -lm #lopenblas 

all : chaos_comparison

clean :
	rm -rf main main.dSYM

main : main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) main.c $(LDFLAGS) -o main

.PHONY: all clean
