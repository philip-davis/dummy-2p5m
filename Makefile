CFLAGS=-g
LDFLAGS=-lm

all: dummy_sim dummy_sense

dummy_sim: dummy_sim.o common.o toml.o
	gcc $^ -o $@ ${LDFLAGS}

dummy_sense: dummy_sense.o common.o toml.o
	gcc $^ -o $@ ${LDFLAGS}

%.o: %.c
	gcc ${CFLAGS} -c $<

clean:
	rm -f *.o dummy_sim dummy_sense 
