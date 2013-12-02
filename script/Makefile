#OPT_FLAGS = -O2 -march=native
OPT_FLAGS = -ggdb
FLAGS = -std=c++11 -fopenmp ${OPT_FLAGS}

hermelin: script.o lexer.o main.o yylexer.o symbol.o node.o node_eval.o calls.o handles.o
	gcc ${FLAGS} -o hermelin script.o lexer.o main.o yylexer.o symbol.o node.o node_eval.o calls.o handles.o -lstdc++ -lm

script.o: script.yy
	bison --defines=tokens.h script.yy
	gcc ${FLAGS} -x c++ -c -o script.o script.tab.cc

lexer.o: lexer.h lexer.cpp
	gcc ${FLAGS} -c -o lexer.o lexer.cpp

yylexer.o: yylexer.h yylexer.cpp
	gcc ${FLAGS} -c -o yylexer.o yylexer.cpp

main.o: main.cpp
	gcc ${FLAGS} -c -o main.o main.cpp

node.o: node.cpp node.h
	gcc ${FLAGS} -c -o node.o node.cpp

node_eval.o: node_eval.cpp node.h
	gcc ${FLAGS} -c -o node_eval.o node_eval.cpp

symbol.o: symbol.cpp symbol.h
	gcc ${FLAGS} -c -o symbol.o symbol.cpp

calls.o: calls.cpp calls.h
	gcc ${FLAGS} -c -o calls.o calls.cpp

handles.o: handles.cpp handles.h
	gcc ${FLAGS} -c -o handles.o handles.cpp

clean:
	rm -f *.o
	rm -f hermelin
	rm -f script.tab.cc
	rm -f tokens.h
