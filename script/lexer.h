/*
 * Lexer
 * @author tweber
 */

#ifndef __MIEZE_LEXER__
#define __MIEZE_LEXER__

#include <istream>
#include <string>
#include <vector>

enum TokenType
{
	LEX_TOKEN_INVALID,
	LEX_TOKEN_END,
	
	LEX_TOKEN_DOUBLE,
	LEX_TOKEN_INT,
	LEX_TOKEN_STRING,
	LEX_TOKEN_IDENT,
	LEX_TOKEN_CHAROP,

	LEX_TOKEN_IF,
	LEX_TOKEN_ELSE,
	LEX_TOKEN_FOR,
	LEX_TOKEN_WHILE,

	LEX_TOKEN_RETURN,
	LEX_TOKEN_BREAK,
	LEX_TOKEN_CONTINUE,

	LEX_TOKEN_LOG_AND,
	LEX_TOKEN_LOG_OR,
	LEX_TOKEN_LOG_NOT,
	LEX_TOKEN_LOG_EQ,
	LEX_TOKEN_LOG_NEQ,
	LEX_TOKEN_LOG_LESS,
	LEX_TOKEN_LOG_GREATER,
	LEX_TOKEN_LOG_LEQ,
	LEX_TOKEN_LOG_GEQ,
};

struct Token
{
	TokenType type;
	
	char cOp;
	int iVal;
	double dVal;
	std::string strVal;
	
	Token()
	{
		type = LEX_TOKEN_INVALID;
		cOp = 0;
		dVal = 0.;
	}
};

class Lexer
{
protected:
	std::string m_strWhitespace, m_strSep;
	
	unsigned int m_iNumToks;
	unsigned int m_iLexPos;
	std::vector<Token> m_vecToks;
	Token m_tokEnd;
	
	void FixTokens();
	
public:
	Lexer();
	Lexer(const std::istream& istr);
	Lexer(const std::string& str);
	virtual ~Lexer();
	
	void load(const std::string& strInput);
	void print();
	const Token& lex();
	
	unsigned int GetNumTokens() const { return m_vecToks.size(); }
	const Token& GetToken(unsigned int i) const { return m_vecToks[i]; }

	static std::string RemoveComments(const std::string& strInput);
	static std::vector<std::string> GetStringTable(const std::string& strInput);
	static void ReplaceExcapes(std::string& str);
};

#endif
