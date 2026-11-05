# Copyright (C) 2026 Antoine Genitrini -- SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE)
# Makefile -- random sampler of Pólya structures
CC      ?= cc
CFLAGS  ?= -O3 -march=native -std=c99 -Wall -Wextra
LDLIBS  = -lm

all: polyagen

polyagen: polyagen.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: polyagen
	./polyagen --selftest

bench: polyagen
	./bench.sh

clean:
	rm -f polyagen

.PHONY: all test bench clean
