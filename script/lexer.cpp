/*
 * Lexer
 * @author tweber
 * @date 2013
 */

#include "lexer.h"
#include <iostream>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <ctype.h>
#include "helper/string.h"

Lexer::Lexer() : m_bOk(1), 
		m_strWhitespace(" \t\r"), m_strSep("=+-*/\%^{}[]();,\":\n"),
		m_iLexPos(0), m_iNumToks(0)
{
	m_tokEnd.type = LEX_TOKEN_END;
}

Lexer::Lexer(const std::string& strInput, const char* pcFile)
		: Lexer()
{
	if(pcFile) m_strFile = pcFile;
	load(strInput);
}

Lexer::~Lexer()
{}

std::string Lexer::RemoveComments(const std::string& strInput)
{
	std::string strRet;
	std::istringstream istr(strInput);

	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		std::size_t iPos = strLine.find('#');

		strLine = strLine.substr(0, iPos);

		strRet += strLine;
		strRet += "\n";
	}

	//std::cout << strRet << std::endl;
	return strRet;
}

static unsigned char ctoi(unsigned char c)
{
	return c-'0';
}

void Lexer::ReplaceEscapes(std::string& str)
{
	find_all_and_replace(str, "\\n", "\n");
	find_all_and_replace(str, "\\t", "\t");
	find_all_and_replace(str, "\\v", "\v");
	find_all_and_replace(str, "\\f", "\f");
	find_all_and_replace(str, "\\b", "\b");
	find_all_and_replace(str, "\\a", "\a");
	find_all_and_replace(str, "\\r", "\r");

// TODO: handle "" in string
	find_all_and_replace(str, "\\\"", "\"");
	find_all_and_replace(str, "\\\'", "\'");
	find_all_and_replace(str, "\\\\", "\\");


	// octal numbers
	if(str.length()>=4) for(int i=0; i<str.length()-3; ++i)
	{
		unsigned char c0 = str[i+1];
		unsigned char c1 = str[i+2];
		unsigned char c2 = str[i+3];

		if(str[i]=='\\' && isdigit(c0) && isdigit(c1) && isdigit(c2))
		{
			//std::cout << "Found: " << str.substr(i, 4) << std::endl;
			char c[2];
			c[0] = ctoi(c0)*8*8 + ctoi(c1)*8 + ctoi(c2);
			c[1] = 0;
			str.replace(i, 4, c);
		}
	}

	// hex numbers
	if(str.length()>=4) for(int i=0; i<str.length()-3; ++i)
	{
		unsigned char c0 = str[i+2];
		unsigned char c1 = str[i+3];

		if(str[i]=='\\' && tolower(str[i+1])=='x' && isdigit(c0) && isdigit(c1))
		{
			//std::cout << "Found: " << str.substr(i, 4) << std::endl;
			char c[2];
			c[0] = ctoi(c0)*16 + ctoi(c1);
			c[1] = 0;
			str.replace(i, 4, c);
		}
	}
}

std::vector<std::string> Lexer::GetStringTable(const std::string& strInput)
{
	std::vector<std::string> vecStr;

	bool bInString = 0;
	std::string str;
	for(char c : strInput)
	{
		if(c == '\"')
		{
			bInString = !bInString;
			if(bInString)
				str = "";
			else
			{
				ReplaceEscapes(str);
				vecStr.push_back(str);
			}
			continue;
		}

		if(bInString)
			str += c;
	}

	//for(const std::string& str : vecStr)
	//	std::cout << "String: " << str << std::endl;
	return vecStr;
}

void Lexer::load(const std::string& _strInput)
{
	std::string strInput = RemoveComments(_strInput);

	// Lexer cannot yet handle \" directly -> replace it
	find_all_and_replace(strInput, "\\\"", "\'");

	std::vector<std::string> vecStr = GetStringTable(strInput);

	typedef boost::char_separator<char> t_sep;
	typedef boost::tokenizer<t_sep> t_tok;
	
	t_sep sep(m_strWhitespace.c_str(), m_strSep.c_str());
	t_tok tok(strInput, sep);

	bool bInString = 0;
	unsigned int iStringIdx = 0;
	unsigned int iCurLine = 1;
	for(const std::string& str : tok)
	{
		if(str.length() == 0) continue;
		if(str == "\n")
		{
			++iCurLine;
			continue;
		}

		if(str=="\"")
		{
			bInString = !bInString;

			if(!bInString)
			{
				Token tokStr;
				tokStr.type = LEX_TOKEN_STRING;
				if(iStringIdx >= vecStr.size())
				{
					m_bOk = 0;
					std::cerr << "Error: String index exceeds string table size."
								<< std::endl;

					continue;
				}
				tokStr.strVal = vecStr[iStringIdx];
				tokStr.iLine = iCurLine;
				m_vecToks.push_back(tokStr);

				++iStringIdx;
			}
			continue;
		}

		if(bInString)
			continue;

		//std::cout << "token: " << str << std::endl;


		Token tok;
		tok.iLine = iCurLine;

		if(str.length()==1 && m_strSep.find(str)!=std::string::npos)
		{
			tok.type = LEX_TOKEN_CHAROP;
			tok.cOp = str[0];
		}
		else if(isalpha(str[0]) || str[0]=='_')
		{
			tok.type = LEX_TOKEN_IDENT;
			tok.strVal = str;

			if(str == "if")				tok.type = LEX_TOKEN_IF;
			else if(str == "else")		tok.type = LEX_TOKEN_ELSE;
			else if(str == "for")		tok.type = LEX_TOKEN_FOR;
			else if(str == "while")		tok.type = LEX_TOKEN_WHILE;
			else if(str == "return")	tok.type = LEX_TOKEN_RETURN;
			else if(str == "break")		tok.type = LEX_TOKEN_BREAK;
			else if(str == "continue")	tok.type = LEX_TOKEN_CONTINUE;
			else if(str == "and")		tok.type = LEX_TOKEN_LOG_AND;
			else if(str == "or")		tok.type = LEX_TOKEN_LOG_OR;
			else if(str == "not")		tok.type = LEX_TOKEN_LOG_NOT;
			else if(str == "eq")		tok.type = LEX_TOKEN_LOG_EQ;
			else if(str == "neq")		tok.type = LEX_TOKEN_LOG_NEQ;
			else if(str == "less")		tok.type = LEX_TOKEN_LOG_LESS;
			else if(str == "greater")	tok.type = LEX_TOKEN_LOG_GREATER;
			else if(str == "leq")		tok.type = LEX_TOKEN_LOG_LEQ;
			else if(str == "geq")		tok.type = LEX_TOKEN_LOG_GEQ;
			else if(str == "global")	tok.type = LEX_TOKEN_GLOBAL;
		}
		else if(isdigit(str[0]))
		{
			tok.type = LEX_TOKEN_DOUBLE;
			tok.strVal = str;
		}
		else
		{
			m_bOk = 0;
			std::string strFile;
			if(m_strFile != "")
				strFile = std::string(" in \"") + m_strFile + "\"";

			std::cerr << "Error (line " << iCurLine << strFile << "): "
					<< "Unknown token: \"" << str << "\"." << std::endl;
			continue;
		}

		m_vecToks.push_back(tok);
	}
	
	FixTokens();
	//print();
}

void Lexer::FixTokens()
{
	for(unsigned int i=0; i<m_vecToks.size(); ++i)
	{
		Token& tok = m_vecToks[i];
		std::string strFullDouble = tok.strVal;
		
		if(tok.type==LEX_TOKEN_DOUBLE)
		{
			if(tok.strVal[tok.strVal.length()-1]=='e' || tok.strVal[tok.strVal.length()-1]=='E')		// 12.3e
			{
				if(m_vecToks[i+1].type==LEX_TOKEN_CHAROP) 						// 12.3e-4
				{
					strFullDouble += m_vecToks[i+1].cOp + m_vecToks[i+2].strVal;
					m_vecToks.erase(m_vecToks.begin()+i+2);
					m_vecToks.erase(m_vecToks.begin()+i+1);
				}
				else if(m_vecToks[i+1].type==LEX_TOKEN_DOUBLE)					// 12.3e3
				{
					strFullDouble += m_vecToks[i+1].strVal;
					m_vecToks.erase(m_vecToks.begin()+i+1);
				}
			}

			// value is actually an int
			if(strFullDouble.find_first_of(".eE") == std::string::npos)
				tok.type = LEX_TOKEN_INT;

			std::istringstream istr(strFullDouble);

			if(tok.type == LEX_TOKEN_DOUBLE)
				istr >> tok.dVal;
			else if(tok.type == LEX_TOKEN_INT)
				istr >> tok.iVal;

			//std::cout << "val " << tok.type << ": " << strFullDouble << std::endl;

			tok.strVal = "";
		}
	}
}

void Lexer::print()
{
	while(1)
	{
		const Token& tok = lex();
		if(tok.type == LEX_TOKEN_END) 
			break;
		
		std::cout << "type: " << tok.type 
					<< ", cOp: " << tok.cOp 
					<< ", dVal: " << tok.dVal 
					<< ", strVal: " << tok.strVal
					<< std::endl;
	}
}

const Token& Lexer::lex()
{	
	if(m_iLexPos >= m_vecToks.size())
		return m_tokEnd;

	return m_vecToks[m_iLexPos++];
}


/*
int main()
{
	Lexer lex("d=5.e-12.3; \n a=2 * d;");
	lex.print();
	
	return 0;
}
*/

