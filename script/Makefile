CC = gcc
#CC = clang
OPT_FLAGS = -O2 -march=native
#OPT_FLAGS = -ggdb -DDEBUG
DEFINES = -DNO_QT -DNO_COMP
INCS = 
FLAGS = -std=c++11 -fopenmp ${OPT_FLAGS} ${DEFINES} ${INCS}
#FLAGS += -Wall -Wpedantic

LIBS = -L/usr/local/lib

# Boost
LIBS = -lboost_iostreams -lboost_regex

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

hermelin: obj/script.o obj/lexer.o obj/main.o obj/yylexer.o obj/symbol.o obj/node.o \
		obj/node_eval.o obj/calls.o obj/handles.o obj/gnuplot.o obj/file.o \
		obj/string.o obj/loadtxt.o obj/rand.o obj/fourier.o obj/script_helper.o \
		obj/calls_basic.o \
		obj/calls_plot.o obj/calls_math.o obj/calls_file.o obj/calls_thread.o \
		obj/calls_fit.o obj/globals.o obj/flags.o obj/fitter.o obj/chi2.o \
		obj/interpolation.o obj/spec_char.o obj/linalg.o obj/linalg2.o \
		obj/log.o
	${CC} ${FLAGS} -o hermelin obj/script.o obj/lexer.o obj/main.o obj/yylexer.o \
		obj/symbol.o obj/node.o obj/node_eval.o obj/calls.o obj/handles.o \
		obj/gnuplot.o obj/file.o obj/string.o obj/loadtxt.o obj/rand.o \
		obj/fourier.o obj/script_helper.o obj/calls_plot.o obj/calls_math.o \
		obj/calls_basic.o \
		obj/calls_file.o obj/calls_thread.o obj/calls_fit.o obj/globals.o \
		obj/flags.o obj/fitter.o obj/chi2.o obj/interpolation.o \
		obj/spec_char.o obj/linalg.o obj/linalg2.o \
		obj/log.o \
			${LIBS}
	strip hermelin

obj/script.o: script.yy
	${YY} --defines=tokens.h script.yy
	${CC} ${FLAGS} -x c++ -c -o obj/script.o script.tab.cc

obj/lexer.o: lexer.h lexer.cpp
	${CC} ${FLAGS} -c -o obj/lexer.o lexer.cpp

obj/yylexer.o: yylexer.h yylexer.cpp obj/script.o
	${CC} ${FLAGS} -c -o obj/yylexer.o yylexer.cpp

obj/main.o: main.cpp
	${CC} ${FLAGS} -c -o obj/main.o main.cpp

obj/node.o: node.cpp node.h
	${CC} ${FLAGS} -c -o obj/node.o node.cpp

obj/node_eval.o: node_eval.cpp node.h
	${CC} ${FLAGS} -c -o obj/node_eval.o node_eval.cpp

obj/symbol.o: symbol.cpp symbol.h
	${CC} ${FLAGS} -c -o obj/symbol.o symbol.cpp


obj/calls.o: calls.cpp calls.h
	${CC} ${FLAGS} -c -o obj/calls.o calls.cpp

obj/calls_basic.o: runtime/calls_basic.cpp runtime/calls_basic.h
	${CC} ${FLAGS} -c -o obj/calls_basic.o runtime/calls_basic.cpp

obj/calls_plot.o: runtime/calls_plot.cpp runtime/calls_plot.h
	${CC} ${FLAGS} -c -o obj/calls_plot.o runtime/calls_plot.cpp

obj/calls_math.o: runtime/calls_math.cpp runtime/calls_math.h
	${CC} ${FLAGS} -c -o obj/calls_math.o runtime/calls_math.cpp

obj/calls_file.o: runtime/calls_file.cpp runtime/calls_file.h
	${CC} ${FLAGS} -c -o obj/calls_file.o runtime/calls_file.cpp

obj/calls_thread.o: runtime/calls_thread.cpp runtime/calls_thread.h
	${CC} ${FLAGS} -c -o obj/calls_thread.o runtime/calls_thread.cpp

obj/calls_fit.o: runtime/calls_fit.cpp runtime/calls_fit.h
	${CC} ${FLAGS} -c -o obj/calls_fit.o runtime/calls_fit.cpp


obj/handles.o: handles.cpp handles.h
	${CC} ${FLAGS} -c -o obj/handles.o handles.cpp

obj/gnuplot.o: helper/gnuplot.cpp helper/gnuplot.h
	${CC} ${FLAGS} -c -o obj/gnuplot.o helper/gnuplot.cpp

obj/file.o: helper/file.cpp helper/file.h
	${CC} ${FLAGS} -c -o obj/file.o helper/file.cpp

obj/string.o: helper/string.cpp helper/string.h
	${CC} ${FLAGS} -c -o obj/string.o helper/string.cpp

obj/spec_char.o: helper/spec_char.cpp helper/spec_char.h
	${CC} ${FLAGS} -c -o obj/spec_char.o helper/spec_char.cpp

obj/loadtxt.o: loader/loadtxt.cpp loader/loadtxt.h
	${CC} ${FLAGS} -c -o obj/loadtxt.o loader/loadtxt.cpp

obj/rand.o: helper/rand.cpp helper/rand.h
	${CC} ${FLAGS} -c -o obj/rand.o helper/rand.cpp

obj/fourier.o: helper/fourier.cpp helper/fourier.h
	${CC} ${FLAGS} -c -o obj/fourier.o helper/fourier.cpp

obj/script_helper.o: script_helper.cpp script_helper.h
	${CC} ${FLAGS} -c -o obj/script_helper.o script_helper.cpp

obj/globals.o: globals.cpp globals.h
	${CC} ${FLAGS} -c -o obj/globals.o globals.cpp

obj/flags.o: helper/flags.cpp helper/flags.h
	${CC} ${FLAGS} -c -o obj/flags.o helper/flags.cpp

obj/log.o: helper/log.cpp helper/log.h
	${CC} ${FLAGS} -c -o obj/log.o helper/log.cpp


obj/fitter.o: fitter/fitter.cpp fitter/fitter.h
	${CC} ${FLAGS} -c -o obj/fitter.o fitter/fitter.cpp

obj/chi2.o: fitter/chi2.cpp fitter/chi2.h
	${CC} ${FLAGS} -c -o obj/chi2.o fitter/chi2.cpp

obj/interpolation.o: fitter/models/interpolation.cpp fitter/models/interpolation.h
	${CC} ${FLAGS} -c -o obj/interpolation.o fitter/models/interpolation.cpp


obj/linalg.o: helper/linalg.cpp helper/linalg.h
	${CC} ${FLAGS} -c -o obj/linalg.o helper/linalg.cpp

obj/linalg2.o: helper/linalg2.cpp helper/linalg2.h
	${CC} ${FLAGS} -c -o obj/linalg2.o helper/linalg2.cpp


clean:
	rm -f obj/*.o
	rm -f hermelin
	rm -f script.tab.cc
	rm -f tokens.h
	rm -f *.hh
