OPT_FLAGS = -O2 -march=native
FLAGS = -std=c++11 ${OPT_FLAGS}

hermelin: script.o lexer.o main.o yylexer.o symbol.o node.o calls.o
	gcc ${FLAGS} -lstdc++ -lm -o hermelin script.o lexer.o main.o yylexer.o symbol.o node.o calls.o

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

symbol.o: symbol.cpp symbol.h
	gcc ${FLAGS} -c -o symbol.o symbol.cpp

calls.o: calls.cpp calls.h
	gcc ${FLAGS} -c -o calls.o calls.cpp

clean:
	rm -f *.o
	rm -f hermelin
	rm -f script.tab.cc
	rm -f tokens.h
