OSTYPE = $(shell uname)
ARCH   = $(shell uname -m)
##################################### Defaults ################################
CC             := gcc
INCLUDE        := -I. -I/usr/local/include
CFLAGS         := -Wall -Wno-overlength-strings -std=c99 -pedantic -O3
CFLAGS_32      := -m32
SHLIB_CFLAGS   := -fPIC -shared
SHLIB_EXT      := .so
LIBS           := -L/opt/local/lib -lm
LDFLAGS        :=
############################# Library add-ons #################################
TINYSCHEME_FEATURES := -DUSE_DL=1 -DUSE_MATH=1 -DUSE_ASCII_NAMES=0
INCLUDE += -I/opt/local/include -I./tinyscheme/trunk
LIBS    += -L/opt/local/lib -lpthread -lhdf5
CFLAGS  += -DH5_NO_DEPRECATED_SYMBOLS
GSLLIBS  = $(shell gsl-config --libs)
LIBS    += $(GSLLIBS)
LIBS    += -lfftw3 -lfftw3_threads
LIBS    += -lreadline ./tinyscheme/trunk/libtinyscheme.a
GLLIBS   =
############################# OS & ARCH specifics #############################
ifneq ($(OSTYPE), Linux)
  ifeq ($(OSTYPE), Darwin)
    CC      = clang
    CFLAGS       += -Wno-gnu-zero-variadic-macro-arguments
    GLLIBS       += -framework GLUT -framework OpenGL -framework Cocoa
    SHLIB_CFLAGS := -dynamiclib
    SHLIB_EXT    := .dylib
    TINYSCHEME_FEATURES += -DUSE_STRLWR=1 -D__APPLE__=1 -DOSX=1
    ifeq ($(shell sysctl -n hw.optional.x86_64), 1)
      ARCH       := x86_64
    endif
  else ifeq ($(OSTYPE), FreeBSD)
    CC      = clang
    GLLIBS += -lGL -lGLU -lglut
  else ifeq ($(OSTYPE), SunOS)
      CFLAGS := -c -Wall -std=c99 -pedantic
  else
      # Let's assume this is win32
      SHLIB_EXT           := .dll
      TINYSCHEME_FEATURES += -DUSE_STRLWR=0
  endif
else
  CFLAGS              += -D_GNU_SOURCE
  TINYSCHEME_FEATURES += -DSUN_DL=1
  GLLIBS              += -lGL -lGLU -lglut
endif

ifneq ($(ARCH), x86_64)
  CFLAGS_32 += -m32
endif

# Are all G5s ppc970s?
ifeq ($(ARCH), ppc970)
  CFLAGS += -m64
endif
############################ Define targets ###################################
EXE_TARGETS = analyzeWaveform waveview spectraFilter
DEBUG_EXE_TARGETS = hdf5rawWaveformIo utils filters runScriptNGetConfig
# SHLIB_TARGETS = XXX$(SHLIB_EXT)

ifeq ($(ARCH), x86_64) # compile a 32bit version on 64bit platforms
  # SHLIB_TARGETS += XXX_m32$(SHLIB_EXT)
endif

.PHONY: exe_targets shlib_targets debug_exe_targets dox clean
exe_targets: $(EXE_TARGETS)
shlib_targets: $(SHLIB_TARGETS)
debug_exe_targets: $(DEBUG_EXE_TARGETS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

spectraFilter: spectraFilter.o filters.o mreadarray.o utils.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
analyzeWaveform: main.o utils.o filters.o peakFinder.o runScriptNGetConfig.o hdf5rawWaveformIo.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
main.o: main.c hdf5rawWaveformIo.h common.h
utils.o: utils.c utils.h common.h
utils: utils.c utils.h common.h
	$(CC) $(CFLAGS) $(INCLUDE) -DUTILS_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
filters.o: filters.c filters.h common.h
filters: filters.c utils.o filters.h common.h
	$(CC) $(CFLAGS) $(INCLUDE) -DFILTERS_DEBUG_ENABLEMAIN $< utils.o $(LIBS) $(LDFLAGS) -o $@
peakFinder.o: peakFinder.c peakFinder.h filters.h common.h
scmfpreg.o: scmfpreg.c scmfpreg.h common.h
	$(CC) $(CFLAGS) $(INCLUDE) $(TINYSCHEME_FEATURES) -c $< -o $@
scmfpreg: scmfpreg.c scmfpreg.h common.h
	$(CC) $(CFLAGS) $(INCLUDE) $(TINYSCHEME_FEATURES) -DSCMFPREG_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
runScriptNGetConfig.o: runScriptNGetConfig.c runScriptNGetConfig.h common.h
	$(CC) $(CFLAGS) $(INCLUDE) $(TINYSCHEME_FEATURES) -c $< -o $@
runScriptNGetConfig: runScriptNGetConfig.c runScriptNGetConfig.h common.h
	$(CC) $(CFLAGS) $(INCLUDE) $(TINYSCHEME_FEATURES) -DRUNSCRIPTNGETCONFIG_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
hdf5rawWaveformIo.o: hdf5rawWaveformIo.c hdf5rawWaveformIo.h common.h
hdf5rawWaveformIo: hdf5rawWaveformIo.c hdf5rawWaveformIo.h
	$(CC) $(CFLAGS) $(INCLUDE) -DHDF5RAWWAVEFORMIO_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
demux: demux.c hdf5rawWaveformIo.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
demuxdump: demuxdump.c hdf5rawWaveformIo.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
waveview: waveview.c hdf5rawWaveformIo.o
	$(CC) $(CFLAGS) $(INCLUDE) -Wno-deprecated-declarations $^ $(LIBS) $(GLLIBS) $(LDFLAGS) -o $@
# libmreadarray$(SHLIB_EXT): mreadarray.o
# 	$(CC) $(SHLIB_CFLAGS) $(CFLAGS) $(LIBS) -o $@ $<
# mreadarray.o: mreadarray.c
# 	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<
# mreadarray: mreadarray.c
# 	$(CC) $(CFLAGS) -DENABLEMAIN $(INCLUDE) $(LIBS) -o $@ $<
# libmreadarray_m32$(SHLIB_EXT): mreadarray.c
# 	$(CC) -m32 $(SHLIB_CFLAGS) $(CFLAGS) $(CFLAGS_32) -o $@ $<

dox:
	doxygen
clean:
	rm -f *.o *.so *.dylib *.dll *.bundle
	rm -f $(SHLIB_TARGETS) $(EXE_TARGETS) $(DEBUG_EXE_TARGETS)
