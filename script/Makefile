CC = gcc
#OPT_FLAGS = -O2 -march=native
OPT_FLAGS = -ggdb -DDEBUG
DEFINES = -DNO_QT -DNO_COMP
FLAGS = -std=c++11 -fopenmp ${OPT_FLAGS} ${DEFINES}
LIBS = -lstdc++ -lm -lboost_iostreams-mt

hermelin: script.o lexer.o main.o yylexer.o symbol.o node.o node_eval.o calls.o handles.o gnuplot.o file.o string.o loadtxt.o rand.o fourier.o script_helper.o
	${CC} ${FLAGS} -o hermelin script.o lexer.o main.o yylexer.o symbol.o node.o node_eval.o calls.o handles.o gnuplot.o file.o string.o loadtxt.o rand.o fourier.o script_helper.o ${LIBS}

script.o: script.yy
	bison --defines=tokens.h script.yy
	${CC} ${FLAGS} -x c++ -c -o script.o script.tab.cc

lexer.o: lexer.h lexer.cpp
	${CC} ${FLAGS} -c -o lexer.o lexer.cpp

yylexer.o: yylexer.h yylexer.cpp
	${CC} ${FLAGS} -c -o yylexer.o yylexer.cpp

main.o: main.cpp
	${CC} ${FLAGS} -c -o main.o main.cpp

node.o: node.cpp node.h
	${CC} ${FLAGS} -c -o node.o node.cpp

node_eval.o: node_eval.cpp node.h
	${CC} ${FLAGS} -c -o node_eval.o node_eval.cpp

symbol.o: symbol.cpp symbol.h
	${CC} ${FLAGS} -c -o symbol.o symbol.cpp

calls.o: calls.cpp calls.h
	${CC} ${FLAGS} -c -o calls.o calls.cpp

handles.o: handles.cpp handles.h
	${CC} ${FLAGS} -c -o handles.o handles.cpp

gnuplot.o: helper/gnuplot.cpp helper/gnuplot.h
	${CC} ${FLAGS} -c -o gnuplot.o helper/gnuplot.cpp

file.o: helper/file.cpp helper/file.h
	${CC} ${FLAGS} -c -o file.o helper/file.cpp

string.o: helper/string.cpp helper/string.h
	${CC} ${FLAGS} -c -o string.o helper/string.cpp

loadtxt.o: loader/loadtxt.cpp loader/loadtxt.h
	${CC} ${FLAGS} -c -o loadtxt.o loader/loadtxt.cpp

rand.o: helper/rand.cpp helper/rand.h
	${CC} ${FLAGS} -c -o rand.o helper/rand.cpp

fourier.o: helper/fourier.cpp helper/fourier.h
	${CC} ${FLAGS} -c -o fourier.o helper/fourier.cpp

script_helper.o: script_helper.cpp script_helper.h
	${CC} ${FLAGS} -c -o script_helper.o script_helper.cpp

clean:
	rm -f *.o
	rm -f hermelin
	rm -f script.tab.cc
	rm -f tokens.h
	rm -f *.hh
