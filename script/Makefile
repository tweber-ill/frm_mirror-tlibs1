CC = gcc
#CC = clang

DEBUG = 0
USE_FFTW = 0

ifeq (${DEBUG}, 0)
	OPT_FLAGS = -O3 -fexpensive-optimizations -march=native -DNDEBUG -Wall
	STRIP = strip
else
	OPT_FLAGS = -ggdb -DDEBUG -rdynamic -Wall
	STRIP = echo debug version: no stripping of
endif

DEFINES = -DNO_QT -DNO_COMP
INCS =
FLAGS = -std=c++14 -fopenmp ${OPT_FLAGS} ${DEFINES} ${INCS}
#FLAGS += -Wall -Wpedantic

LIBS = -L/usr/local/lib

# Boost
LIBS = -lboost_iostreams

# FFTW
ifeq (${USE_FFTW}, 1)
	DEFINES += -DUSE_FFTW
	LIBS += -lfftw3
	FFTW_OBJ = obj/fftw.o
else
	FFTW_OBJ =
endif

# Minuit
INCS += -I/usr/include/root
LIBS += -L/usr/lib64/root -lMinuit2

# std
LIBS += -lstdc++ -lm

# Lapack/e
INCS += -I/usr/include/lapacke
LIBS += -llapack -llapacke

YY = bison
#YY = /usr/local/bin/bison

.PHONY: all clean

all: hermelin


hermelin: obj/script.o obj/lexer.o obj/main.o obj/yylexer.o obj/symbol.o obj/node.o \
		obj/node_eval.o obj/calls.o obj/handles.o obj/gnuplot.o obj/file.o \
		obj/string.o obj/loadtxt.o obj/rand.o obj/fourier.o obj/script_helper.o \
		obj/calls_basic.o obj/info.o \
		obj/calls_plot.o obj/calls_math.o obj/calls_file.o obj/calls_thread.o \
		obj/calls_fit.o obj/globals.o obj/flags.o obj/minuit.o \
		obj/interpolation.o obj/spec_char.o obj/linalg2.o \
		obj/log.o obj/node_opt.o ${FFTW_OBJ}
	${CC} ${FLAGS} -o $@ $+ ${LIBS}
	${STRIP} $@

obj/script.o: lang/script.yy
	${YY} --defines=tokens.h --output-file=script.cc lang/script.yy
	${CC} ${FLAGS} -x c++ -c -o obj/script.o script.cc
obj/lexer.o: lang/lexer.cpp lang/lexer.h
	${CC} ${FLAGS} -c -o $@ $<
obj/yylexer.o: lang/yylexer.cpp lang/yylexer.h obj/script.o
	${CC} ${FLAGS} -c -o $@ $<
obj/main.o: lang/main.cpp
	${CC} ${FLAGS} -c -o $@ $<
obj/info.o: lang/info.cpp lang/info.h
	${CC} ${FLAGS} -c -o $@ $<
obj/node.o: lang/node.cpp lang/node.h
	${CC} ${FLAGS} -c -o $@ $<
obj/node_eval.o: lang/node_eval.cpp lang/node.h
	${CC} ${FLAGS} -c -o $@ $<
obj/symbol.o: lang/symbol.cpp lang/symbol.h
	${CC} ${FLAGS} -c -o $@ $<
obj/node_opt.o: lang/node_opt.cpp lang/node.h
	${CC} ${FLAGS} -c -o $@ $<
obj/calls.o: lang/calls.cpp lang/calls.h
	${CC} ${FLAGS} -c -o $@ $<
obj/calls_basic.o: runtime/calls_basic.cpp runtime/calls_basic.h
	${CC} ${FLAGS} -c -o $@ $<
obj/calls_plot.o: runtime/calls_plot.cpp runtime/calls_plot.h
	${CC} ${FLAGS} -c -o $@ $<
obj/calls_math.o: runtime/calls_math.cpp runtime/calls_math.h
	${CC} ${FLAGS} -c -o $@ $<
obj/calls_file.o: runtime/calls_file.cpp runtime/calls_file.h
	${CC} ${FLAGS} -c -o $@ $<
obj/calls_thread.o: runtime/calls_thread.cpp runtime/calls_thread.h
	${CC} ${FLAGS} -c -o $@ $<
obj/calls_fit.o: runtime/calls_fit.cpp runtime/calls_fit.h
	${CC} ${FLAGS} -c -o $@ $<
	${CC} ${FLAGS} -c -o obj/calls_fit.o runtime/calls_fit.cpp
obj/handles.o: lang/handles.cpp lang/handles.h
	${CC} ${FLAGS} -c -o $@ $<
obj/gnuplot.o: ../gfx/gnuplot.cpp ../gfx/gnuplot.h
	${CC} ${FLAGS} -c -o $@ $<
obj/file.o: ../file/file.cpp ../file/file.h
	${CC} ${FLAGS} -c -o $@ $<
obj/string.o: ../string/string.cpp ../string/string.h
	${CC} ${FLAGS} -c -o $@ $<
obj/spec_char.o: ../string/spec_char.cpp ../string/spec_char.h
	${CC} ${FLAGS} -c -o $@ $<
obj/loadtxt.o: ../file/loadtxt.cpp ../file/loadtxt.h
	${CC} ${FLAGS} -c -o $@ $<
obj/rand.o: ../math/rand.cpp ../math/rand.h
	${CC} ${FLAGS} -c -o $@ $<
obj/fourier.o: ../math/fourier.cpp ../math/fourier.h ../math/dft.h
	${CC} ${FLAGS} -c -o $@ $<
obj/fftw.o: ../math/fftw.cpp ../math/fftw.h ../math/dft.h
	${CC} ${FLAGS} -c -o $@ $<
obj/script_helper.o: lang/script_helper.cpp lang/script_helper.h
	${CC} ${FLAGS} -c -o $@ $<
obj/globals.o: lang/globals.cpp lang/globals.h
	${CC} ${FLAGS} -c -o $@ $<
obj/flags.o: ../helper/flags.cpp ../helper/flags.h
	${CC} ${FLAGS} -c -o $@ $<
obj/log.o: ../helper/log.cpp ../helper/log.h
	${CC} ${FLAGS} -c -o $@ $<
obj/minuit.o: ../math/minuit.cpp ../math/minuit.h
	${CC} ${FLAGS} -c -o $@ $<
obj/interpolation.o: ../math/interpolation.cpp ../math/interpolation.h
	${CC} ${FLAGS} -c -o $@ $<
obj/linalg2.o: ../math/linalg2.cpp ../math/linalg2.h
	${CC} ${FLAGS} -c -o $@ $<


clean:
	rm -f obj/*.o
	rm -f hermelin
	rm -f script.cc
	rm -f script.tab.cc
	rm -f tokens.h
	rm -f *.hh
