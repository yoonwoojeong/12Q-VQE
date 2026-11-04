# qchem-vqe build & deployment
CC      ?= gcc
PREFIX  ?= /usr/local
GIT_REV := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
LIBCINT ?= third_party/libcint
SANITIZE ?= 0

CFLAGS  := -std=c17 -O2 -g -fPIC -fopenmp -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare \
           -Iinclude -I$(LIBCINT)/include -I/usr/include/x86_64-linux-gnu \
           -DQC_GIT_REV=\"$(GIT_REV)\" -DQC_DATA_DIR=\"$(CURDIR)/data\"
RPATH_FLAGS := -Wl,-rpath,'$$ORIGIN/../lib' -Wl,-rpath,'$$ORIGIN/lib' -Wl,-rpath,$(CURDIR)/$(LIBCINT)/lib
LDLIBS  := -fopenmp -L$(LIBCINT)/lib $(RPATH_FLAGS) -lcint -lquadmath \
           -llapacke -lopenblas -lnlopt -lcjson -lcurl -lm

ifeq ($(SANITIZE),1)
CFLAGS  += -fsanitize=address,undefined -fno-omit-frame-pointer
LDLIBS  += -fsanitize=address,undefined
endif

SRC_LIB := $(shell find src -name '*.c' ! -name main.c)
OBJ_LIB := $(SRC_LIB:src/%.c=build/obj/%.o)
OBJ_MAIN := build/obj/main.o
TESTS   := $(wildcard tests/test_*.c)
TEST_BIN := $(TESTS:tests/%.c=build/tests/%)

.PHONY: all test sanitize install dist clean fetch-libcint
all: build/qchem-vqe build/libqchem.a build/libqchem.so

fetch-libcint:
	tr -d '\r' < scripts/fetch_libcint.sh | sh

build/obj/%.o: src/%.c include/qchem.h
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

build/libqchem.a: $(OBJ_LIB)
	ar rcs $@ $^

build/libqchem.so: $(OBJ_LIB)
	$(CC) -shared $(CFLAGS) -o $@ $(OBJ_LIB) $(LDLIBS)

build/qchem-vqe: $(OBJ_MAIN) build/libqchem.a
	$(CC) $(CFLAGS) -o $@ $(OBJ_MAIN) build/libqchem.a $(LDLIBS)

build/tests/%: tests/%.c build/libqchem.a
	@mkdir -p build/tests
	$(CC) $(CFLAGS) -o $@ $< build/libqchem.a $(LDLIBS)

test: $(TEST_BIN)
	@fail=0; for t in $(TEST_BIN); do echo "== $$t"; ./$$t || fail=1; done; \
	if [ $$fail -ne 0 ]; then echo "SOME TESTS FAILED"; exit 1; else echo "ALL TESTS PASSED"; fi

sanitize:
	$(MAKE) clean
	$(MAKE) test SANITIZE=1
	$(MAKE) clean
	$(MAKE) all

install: all
	install -d $(DESTDIR)$(PREFIX)/bin
	install -d $(DESTDIR)$(PREFIX)/lib/pkgconfig
	install -d $(DESTDIR)$(PREFIX)/include
	install -d $(DESTDIR)$(PREFIX)/share/qchem-vqe/data/basis
	install -m 0755 build/qchem-vqe $(DESTDIR)$(PREFIX)/bin/qchem-vqe
	install -m 0644 build/libqchem.a $(DESTDIR)$(PREFIX)/lib/libqchem.a
	install -m 0755 build/libqchem.so $(DESTDIR)$(PREFIX)/lib/libqchem.so
	cp -a $(LIBCINT)/lib/libcint.so* $(DESTDIR)$(PREFIX)/lib/
	install -m 0644 include/qchem.h $(DESTDIR)$(PREFIX)/include/qchem.h
	install -m 0644 data/basis/*.nw $(DESTDIR)$(PREFIX)/share/qchem-vqe/data/basis/
	@printf "prefix=%s\nexec_prefix=\$${prefix}\nlibdir=\$${exec_prefix}/lib\nincludedir=\$${prefix}/include\n\nName: qchem\nDescription: Native C17 Quantum Chemistry VQE, ADAPT-VQE, and SQD Engine\nVersion: 0.2.0\nLibs: -L\$${libdir} -Wl,-rpath,\$${libdir} -lqchem -lcint -llapacke -lopenblas -lnlopt -lcjson -lcurl -lm -fopenmp\nCflags: -I\$${includedir} -fopenmp\n" "$(PREFIX)" > $(DESTDIR)$(PREFIX)/lib/pkgconfig/qchem.pc

dist: all
	@rm -rf build/dist && mkdir -p build/dist/qchem-vqe/bin build/dist/qchem-vqe/lib build/dist/qchem-vqe/include build/dist/qchem-vqe/data/basis
	cp build/qchem-vqe build/dist/qchem-vqe/bin/
	cp build/libqchem.a build/libqchem.so build/dist/qchem-vqe/lib/
	cp -a $(LIBCINT)/lib/libcint.so* build/dist/qchem-vqe/lib/
	cp include/qchem.h build/dist/qchem-vqe/include/
	cp data/basis/*.nw build/dist/qchem-vqe/data/basis/
	cp README.md LICENSE build/dist/qchem-vqe/ 2>/dev/null || true
	tar -C build/dist -czf build/qchem-vqe-linux-x86_64.tar.gz qchem-vqe
	@echo "Created portable release tarball: build/qchem-vqe-linux-x86_64.tar.gz"

clean:
	rm -rf build/obj build/tests build/qchem-vqe build/libqchem.a build/libqchem.so build/dist build/qchem-vqe-linux-x86_64.tar.gz
