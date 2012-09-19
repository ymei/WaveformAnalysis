CC=gcc
CFLAGS=-Wall -O2 -DDEBUG
INCLUDE=-I/opt/local/include -I.
LIBS=-L/opt/local/lib -lpthread -lhdf5

GSLLIBS = $(shell gsl-config --libs)

LIBS   += $(GSLLIBS)
LIBS   += -lfftw3

.PHONY: all clean
all: analyzeWaveform wavedump analyze_pe
analyzeWaveform: main.c hdf5RawWaveformIo.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
analyze_pe: analysis/analyze_pe.c hdf5io.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
analyze_int: analysis/analyze_int.c hdf5io.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
wavedump: analysis/wavedump.c hdf5io.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
hdf5rawWaveformIo.o: hdf5rawWaveformIo.c hdf5rawWaveformIo.h
	$(CC) $(CFLAGS) -DH5_NO_DEPRECATED_SYMBOLS $(INCLUDE) -c $<
hdf5rawWaveformIo: hdf5rawWaveformIo.c hdf5rawWaveformIo.h
	$(CC) $(CFLAGS) -DH5_NO_DEPRECATED_SYMBOLS $(INCLUDE) -DHDF5RAWWAVEFORMIO_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
clean:
	rm -f *.o

