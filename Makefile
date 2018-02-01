CC=gcc
CFLAGS=-O3 -Wall
# don't have time to chase all those up
# CFLAGS=-O3
DBGCFLAGS=-g -Wall
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!

EXES=gffsimp astarp

gffsimp: gffsimp.c
	${CC} ${DBGCFLAGS} -o $@ $^

# derived from gffsimp, reading in array star novoaign style text files
astarp: astarp.c
	${CC} ${DBGCFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
