CFLAGS=-g
LDFLAGS=-lm

all: dummy_sim dummy_sense

dummy_sim: dummy_sim.o common.o toml.o
	mpicc $^ -o $@ ${LDFLAGS}

dummy_sense: dummy_sense.o common.o toml.o
	mpicc $^ -o $@ ${LDFLAGS}

%.o: %.c
	mpicc ${CFLAGS} -c $<

clean:
	rm -f *.o dummy_sim dummy_sense 
