CC=gcc
CFLAGS=-Wall -O2 -DDEBUG
INCLUDE=-I/opt/local/include -I.
LIBS=-L/opt/local/lib -lpthread -lhdf5

.PHONY: all clean
all: dpo5054 wavedump analyze_pe
dpo5054: main.c hdf5io.o fifo.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
analyze_pe: analysis/analyze_pe.c hdf5io.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
analyze_int: analysis/analyze_int.c hdf5io.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
wavedump: analysis/wavedump.c hdf5io.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
hdf5io.o: hdf5io.c hdf5io.h
	$(CC) $(CFLAGS) -DH5_NO_DEPRECATED_SYMBOLS $(INCLUDE) -c $<
hdf5io: hdf5io.c hdf5io.h
	$(CC) $(CFLAGS) -DH5_NO_DEPRECATED_SYMBOLS $(INCLUDE) -DHDF5IO_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
fifo.o: fifo.c fifo.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $<
fifo_test: fifo.c fifo.h
	$(CC) $(CFLAGS) $(INCLUDE) -DFIFO_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
clean:
	rm -f *.o
