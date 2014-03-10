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
#include "helper/spec_char.h"

Lexer::Lexer() : m_bOk(1), 
		m_strWhitespace(T_STR" \t\r"), m_strSep(T_STR"=+-*/\%^{}[]();,\":\n"),
		m_iLexPos(0), m_iNumToks(0)
{
	m_tokEnd.type = LEX_TOKEN_END;
}

Lexer::Lexer(const t_string& strInput, const t_char* pcFile)
		: Lexer()
{
	if(pcFile) m_strFile = pcFile;
	load(strInput);
}

Lexer::~Lexer()
{}

t_string Lexer::RemoveComments(const t_string& strInput)
{
	t_string strRet;
	t_istringstream istr(strInput);

	while(!istr.eof())
	{
		t_string strLine;
		std::getline(istr, strLine);
		std::size_t iPos = strLine.find('#');

		strLine = strLine.substr(0, iPos);

		strRet += strLine;
		strRet += T_STR"\n";
	}

	//std::cout << strRet << std::endl;
	return strRet;
}

static unsigned char ctoi(unsigned char c)
{
	return c-'0';
}

void Lexer::ReplaceEscapes(t_string& str)
{
	const t_mapSpecChars& mapSpec = get_spec_chars();

	static bool s_bEscapesInited = 0;
	static std::map<t_string, t_string> s_mapstrSpecial;
	static std::map<t_string, t_string> s_mapstrEscapes
	{
		{ t_string(T_STR"\\n"), t_string(T_STR"\n") },
		{ t_string(T_STR"\\t"), t_string(T_STR"\t") },
		{ t_string(T_STR"\\v"), t_string(T_STR"\v") },
		{ t_string(T_STR"\\f"), t_string(T_STR"\f") },
		{ t_string(T_STR"\\b"), t_string(T_STR"\b") },
		{ t_string(T_STR"\\a"), t_string(T_STR"\a") },
		{ t_string(T_STR"\\r"), t_string(T_STR"\r") },

		// TODO: handle "" in string
		{ t_string(T_STR"\\\""), t_string(T_STR"\"") },
		{ t_string(T_STR"\\\'"), t_string(T_STR"\'") },
		{ t_string(T_STR"\\\\"), t_string(T_STR"\\") }
	};

	if(!s_bEscapesInited)
	for(const t_mapSpecChars::value_type& pair : mapSpec)
	{
		const t_string strKey = T_STR("\\") + pair.first;
		const t_string& strVal = T_STR(pair.second.strUTF8);

		s_mapstrSpecial.insert(std::pair<t_string, t_string>(strKey, strVal));
		s_bEscapesInited = 1;
	}

	for(const auto& pair : s_mapstrSpecial)
		find_all_and_replace(str, pair.first, pair.second);

	for(const auto& pair : s_mapstrEscapes)
		find_all_and_replace(str, pair.first, pair.second);




	// octal numbers
	if(str.length()>=4) for(int i=0; i<str.length()-3; ++i)
	{
		t_uchar c0 = str[i+1];
		t_uchar c1 = str[i+2];
		t_uchar c2 = str[i+3];

		if(str[i]=='\\' && isdigit(c0) && isdigit(c1) && isdigit(c2))
		{
			//std::cout << "Found: " << str.substr(i, 4) << std::endl;
			t_char c[2];
			c[0] = ctoi(c0)*8*8 + ctoi(c1)*8 + ctoi(c2);
			c[1] = 0;
			str.replace(i, 4, c);
		}
	}

	// hex numbers
	if(str.length()>=4) for(int i=0; i<str.length()-3; ++i)
	{
		t_uchar c0 = str[i+2];
		t_uchar c1 = str[i+3];

		if(str[i]=='\\' && tolower(str[i+1])=='x' && isdigit(c0) && isdigit(c1))
		{
			//std::cout << "Found: " << str.substr(i, 4) << std::endl;
			t_char c[2];
			c[0] = ctoi(c0)*16 + ctoi(c1);
			c[1] = 0;
			str.replace(i, 4, c);
		}
	}
}

std::vector<t_string> Lexer::GetStringTable(const t_string& strInput)
{
	std::vector<t_string> vecStr;

	bool bInString = 0;
	t_string str;
	for(t_char c : strInput)
	{
		if(c == '\"')
		{
			bInString = !bInString;
			if(bInString)
				str = T_STR"";
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

void Lexer::load(const t_string& _strInput)
{
	t_string strInput = RemoveComments(_strInput);
	//G_COUT << strInput << std::endl;

	// Lexer cannot yet handle \" directly -> replace it
	find_all_and_replace(strInput, t_string(T_STR"\\\""), t_string(T_STR"\'"));

	std::vector<t_string> vecStr = GetStringTable(strInput);

	typedef boost::char_separator<t_char> t_sep;
	typedef boost::tokenizer<t_sep, t_string::const_iterator, t_string> t_tok;
	
	t_sep sep(m_strWhitespace.c_str(), m_strSep.c_str());
	t_tok tok(strInput, sep);

	bool bInString = 0;
	unsigned int iStringIdx = 0;
	unsigned int iCurLine = 1;
	for(const t_string& str : tok)
	{
		if(str.length() == 0) continue;
		if(str == T_STR"\n")
		{
			++iCurLine;
			continue;
		}

		if(str == T_STR"\"")
		{
			bInString = !bInString;

			if(!bInString)
			{
				Token tokStr;
				tokStr.type = LEX_TOKEN_STRING;
				if(iStringIdx >= vecStr.size())
				{
					m_bOk = 0;
					G_CERR << "Error: String index exceeds string table size."
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

		if(str.length()==1 && m_strSep.find(str)!=t_string::npos)
		{
			tok.type = LEX_TOKEN_CHAROP;
			tok.cOp = str[0];
		}
		else if(isalpha(str[0]) || str[0]=='_')
		{
			tok.type = LEX_TOKEN_IDENT;
			tok.strVal = str;

			if(str == T_STR"if")			tok.type = LEX_TOKEN_IF;
			else if(str == T_STR"else")		tok.type = LEX_TOKEN_ELSE;
			else if(str == T_STR"for")		tok.type = LEX_TOKEN_FOR;
			else if(str == T_STR"while")	tok.type = LEX_TOKEN_WHILE;
			else if(str == T_STR"return")	tok.type = LEX_TOKEN_RETURN;
			else if(str == T_STR"break")	tok.type = LEX_TOKEN_BREAK;
			else if(str == T_STR"continue")	tok.type = LEX_TOKEN_CONTINUE;
			else if(str == T_STR"and")		tok.type = LEX_TOKEN_LOG_AND;
			else if(str == T_STR"or")		tok.type = LEX_TOKEN_LOG_OR;
			else if(str == T_STR"not")		tok.type = LEX_TOKEN_LOG_NOT;
			else if(str == T_STR"eq")		tok.type = LEX_TOKEN_LOG_EQ;
			else if(str == T_STR"neq")		tok.type = LEX_TOKEN_LOG_NEQ;
			else if(str == T_STR"less")		tok.type = LEX_TOKEN_LOG_LESS;
			else if(str == T_STR"greater")	tok.type = LEX_TOKEN_LOG_GREATER;
			else if(str == T_STR"leq")		tok.type = LEX_TOKEN_LOG_LEQ;
			else if(str == T_STR"geq")		tok.type = LEX_TOKEN_LOG_GEQ;
			else if(str == T_STR"globaT_STR")	tok.type = LEX_TOKEN_GLOBAL;
		}
		else if(isdigit(str[0]))
		{
			tok.type = LEX_TOKEN_DOUBLE;
			tok.strVal = str;
		}
		else
		{
			m_bOk = 0;
			t_string strFile;
			if(m_strFile != T_STR"")
				strFile = t_string(T_STR" in \"") + m_strFile + T_STR"\"";

			G_CERR << "Error (line " << iCurLine << strFile << "): "
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
		t_string strFullDouble = tok.strVal;
		
		if(tok.type==LEX_TOKEN_DOUBLE)
		{
			if(tok.strVal[tok.strVal.length()-1]=='e' ||
					tok.strVal[tok.strVal.length()-1]=='E')		// 12.3e
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
			if(strFullDouble.find_first_of(T_STR".eE") == t_string::npos)
				tok.type = LEX_TOKEN_INT;

			t_istringstream istr(strFullDouble);

			if(tok.type == LEX_TOKEN_DOUBLE)
				istr >> tok.dVal;
			else if(tok.type == LEX_TOKEN_INT)
				istr >> tok.iVal;

			//std::cout << "val " << tok.type << ": " << strFullDouble << std::endl;

			tok.strVal = T_STR"";
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
		
		G_COUT << "type: " << tok.type
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
// gcc -o tst0 lexer.cpp -lstdc++ -std=c++11
int main()
{
	Lexer lex(T_STR"d=5.6e-12.3; \n a=2 * d;");
	lex.print();

	return 0;
}
*/
