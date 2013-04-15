FLAGS = -std=c++11

script: script.o lexer.o main.o yylexer.o symbol.o node.o
	gcc ${FLAGS} -lstdc++ -o script script.o lexer.o main.o yylexer.o symbol.o node.o

script.o: script.yy
	bison --defines=tokens.h script.yy
	gcc ${FLAGS} -x c++ -c -o script.o script.tab.cc

lexer.o: lexer.h lexer.cpp
	gcc ${FLAGS} -c -o lexer.o lexer.cpp

yylexer.o: yylexer.h yylexer.cpp tokens.h
	gcc ${FLAGS} -c -o yylexer.o yylexer.cpp

main.o: main.cpp
	gcc ${FLAGS} -c -o main.o main.cpp

node.o: node.cpp node.h
	gcc ${FLAGS} -c -o node.o node.cpp

symbol.o: symbol.cpp symbol.h
	gcc ${FLAGS} -c -o symbol.o symbol.cpp

clean:
	rm -f *.o
	rm -f script
	rm -f script.tab.cc
	rm -f tokens.h
