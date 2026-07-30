# AR1KLT -- exact fast AR(1) KLT modules, N = 2..8
CC      ?= cc
CFLAGS  ?= -std=c89 -pedantic -Wall -Wextra -O2
CPPFLAGS = -Iinclude
LDLIBS   = -lm

all: test_ar1klt

ar1klt.o: src/ar1klt.c include/ar1klt.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/ar1klt.c -o ar1klt.o

test_ar1klt: test/test_ar1klt.c ar1klt.o include/ar1klt.h
	$(CC) $(CFLAGS) $(CPPFLAGS) test/test_ar1klt.c ar1klt.o -o test_ar1klt $(LDLIBS)

test: test_ar1klt
	./test_ar1klt

docs:
	doxygen Doxyfile

clean:
	rm -f ar1klt.o test_ar1klt
	rm -rf docs/html

.PHONY: all test docs clean
