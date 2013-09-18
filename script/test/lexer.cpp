/*
 * Lexer
 * @author tweber
 */

#include "lexer.h"
#include <iostream>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <ctype.h>


Lexer::Lexer() : m_strWhitespace(" \t\n\r"), m_strSep("=+-*/^{}();,\""),
				m_iLexPos(0), m_iNumToks(0)
{
	m_tokEnd.type = LEX_TOKEN_END;
}

Lexer::Lexer(const std::istream& istr) : Lexer()
{
	// TODO
}

Lexer::Lexer(const std::string& strInput) : Lexer()
{
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

static void find_and_replace(std::string& str1, const std::string& str_old,
                                                const std::string& str_new)
{
	std::size_t pos = str1.find(str_old);
	if(pos==std::string::npos)
			return;
	str1.replace(pos, str_old.length(), str_new);
}

void Lexer::ReplaceExcapes(std::string& str)
{
	find_and_replace(str, "\\n", "\n");
	find_and_replace(str, "\\t", "\t");
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
				ReplaceExcapes(str);
				vecStr.push_back(str);
			}
			continue;
		}

		if(bInString)
			str += c;
	}

	for(const std::string& str : vecStr)
		std::cout << "String: " << str << std::endl;
	return vecStr;
}

void Lexer::load(const std::string& _strInput)
{
	std::string strInput = RemoveComments(_strInput);
	std::vector<std::string> vecStr = GetStringTable(strInput);

	typedef boost::char_separator<char> t_sep;
	typedef boost::tokenizer<t_sep> t_tok;
	
	t_sep sep(m_strWhitespace.c_str(), m_strSep.c_str());
	t_tok tok(strInput, sep);

	bool bInString = 0;
	unsigned int iStringIdx = 0;
	for(const std::string& str : tok)
	{
		if(str.length() == 0) continue;

		if(str=="\"")
		{
			bInString = !bInString;

			if(!bInString)
			{
				Token tokStr;
				tokStr.type = LEX_TOKEN_STRING;
				if(iStringIdx >= vecStr.size())
				{
					std::cerr << "Error: String index exceeds string table size."
								<< std::endl;
					continue;
				}
				tokStr.strVal = vecStr[iStringIdx];
				m_vecToks.push_back(tokStr);

				++iStringIdx;
			}
			continue;
		}

		if(bInString)
			continue;

		//std::cout << "token: " << str << std::endl;


		Token tok;

		if(str.length()==1 && m_strSep.find(str)!=std::string::npos)
		{
			tok.type = LEX_TOKEN_CHAROP;
			tok.cOp = str[0];
		}
		else if(isalpha(str[0]) || str[0]=='_')
		{
			tok.type = LEX_TOKEN_IDENT;
			tok.strVal = str;
		}
		else if(isdigit(str[0]))
		{
			tok.type = LEX_TOKEN_DOUBLE;
			tok.strVal = str;
		}
		else
		{
			std::cerr << "Unknown token: " << str << std::endl;
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

