# qchem-vqe build
CC      ?= gcc
GIT_REV := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
LIBCINT ?= third_party/libcint

CFLAGS  := -std=c17 -O2 -g -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare \
           -Iinclude -I$(LIBCINT)/include -I/usr/include/x86_64-linux-gnu \
           -DQC_GIT_REV=\"$(GIT_REV)\" -DQC_DATA_DIR=\"$(CURDIR)/data\"
LDLIBS  := -L$(LIBCINT)/lib -Wl,-rpath,$(CURDIR)/$(LIBCINT)/lib -lcint -lquadmath \
           -llapacke -lopenblas -lnlopt -lcjson -lcurl -lm

SRC_LIB := $(shell find src -name '*.c' ! -name main.c)
OBJ_LIB := $(SRC_LIB:src/%.c=build/obj/%.o)
OBJ_MAIN := build/obj/main.o
TESTS   := $(wildcard tests/test_*.c)
TEST_BIN := $(TESTS:tests/%.c=build/tests/%)

.PHONY: all test clean fetch-libcint
all: build/qchem-vqe

fetch-libcint:
	sh scripts/fetch_libcint.sh

build/obj/%.o: src/%.c include/qchem.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

build/libqchem.a: $(OBJ_LIB)
	ar rcs $@ $^

build/qchem-vqe: $(OBJ_MAIN) build/libqchem.a
	$(CC) $(CFLAGS) -o $@ $(OBJ_MAIN) build/libqchem.a $(LDLIBS)

build/tests/%: tests/%.c build/libqchem.a
	@mkdir -p build/tests
	$(CC) $(CFLAGS) -o $@ $< build/libqchem.a $(LDLIBS)

test: $(TEST_BIN)
	@fail=0; for t in $(TEST_BIN); do echo "== $$t"; ./$$t || fail=1; done; \
	if [ $$fail -ne 0 ]; then echo "SOME TESTS FAILED"; exit 1; else echo "ALL TESTS PASSED"; fi

clean:
	rm -rf build
