# AR1KLT -- exact fast AR(1) KLT: fixed modules N = 2..8 and the
# general-N recursive factorization (strict C89)
CC      ?= cc
CFLAGS  ?= -std=c89 -pedantic -Wall -Wextra -O2
CPPFLAGS = -Iinclude
LDLIBS   = -lm

all: test_ar1klt test_ar1klt_gen

ar1klt.o: src/ar1klt.c include/ar1klt.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/ar1klt.c -o ar1klt.o

ar1klt_gen.o: src/ar1klt_gen.c include/ar1klt.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/ar1klt_gen.c -o ar1klt_gen.o

test_ar1klt: test/test_ar1klt.c ar1klt.o include/ar1klt.h
	$(CC) $(CFLAGS) $(CPPFLAGS) test/test_ar1klt.c ar1klt.o -o test_ar1klt $(LDLIBS)

test_ar1klt_gen: test/test_ar1klt_gen.c ar1klt.o ar1klt_gen.o include/ar1klt.h
	$(CC) $(CFLAGS) $(CPPFLAGS) test/test_ar1klt_gen.c ar1klt.o ar1klt_gen.o -o test_ar1klt_gen $(LDLIBS)

test: test_ar1klt test_ar1klt_gen
	./test_ar1klt
	./test_ar1klt_gen

docs:
	doxygen Doxyfile

clean:
	rm -f ar1klt.o ar1klt_gen.o test_ar1klt test_ar1klt_gen
	rm -rf docs/html

.PHONY: all test docs clean
