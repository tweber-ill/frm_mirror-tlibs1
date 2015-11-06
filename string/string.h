/*
 * string helper
 * @author tweber
 * @date 06-mar-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIB_STRINGS__
#define __TLIB_STRINGS__

#include <string>
#include <cstring>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <sstream>
#include <locale>
#include <map>

#include "../helper/exception.h"
#include "../helper/misc.h"


namespace tl {

// -----------------------------------------------------------------------------

template<class t_str=std::string> const t_str& get_dir_seps();
template<class t_str=std::string> const t_str& get_trim_chars();

template<> inline const std::string& get_dir_seps()
{
	static const std::string strSeps("\\/");
	return strSeps;
}
template<> inline const std::wstring& get_dir_seps()
{
	static const std::wstring strSeps(L"\\/");
	return strSeps;
}

template<> inline const std::string& get_trim_chars()
{
	static const std::string strC(" \t\r");
	return strC;
}
template<> inline const std::wstring& get_trim_chars()
{
	static const std::wstring strC(L" \t\r");
	return strC;
}


// -----------------------------------------------------------------------------


static inline std::wstring str_to_wstr(const std::string& str)
{
	return std::wstring(str.begin(), str.end());
}

static inline std::string wstr_to_str(const std::wstring& str)
{
	return std::string(str.begin(), str.end());
}


// -----------------------------------------------------------------------------

template<class t_str=std::string>
t_str get_file_noext(const t_str& str)
{
	std::size_t iPos = str.find_last_of('.');

	if(iPos == t_str::npos)
		return str;
	return str.substr(0, iPos);
}

template<class t_str=std::string>
t_str get_fileext(const t_str& str)
{
	std::size_t iPos = str.find_last_of('.');

	if(iPos == t_str::npos)
		return t_str();
	return str.substr(iPos+1);
}

// e.g. returns "tof" for "123.tof.bz2"
template<class t_str=std::string>
t_str get_fileext2(const t_str& str)
{
	std::size_t iPos = str.find_last_of('.');
	if(iPos == t_str::npos || iPos == 0)
		return t_str();

	t_str strFile = str.substr(0, iPos);
	return get_fileext(strFile);
}

template<class t_str=std::string>
t_str get_dir(const t_str& str)
{
	std::size_t iPos = str.find_last_of(get_dir_seps<t_str>());

	if(iPos == t_str::npos)
		return t_str();
	return str.substr(0, iPos);
}

template<class t_str=std::string>
t_str get_file(const t_str& str)
{
	std::size_t iPos = str.find_last_of(get_dir_seps<t_str>());

	if(iPos == t_str::npos)
		return t_str();
	return str.substr(iPos+1);
}


// -----------------------------------------------------------------------------


template<class t_str=std::string>
t_str insert_before(const t_str& str,
					const t_str& strChar,
					const t_str& strInsert)
{
	std::size_t pos = str.find(strChar);
	if(pos==t_str::npos)
			return str;

	t_str strRet = str;
	strRet.insert(pos, strInsert);

	return strRet;
}

template<class t_str=std::string>
t_str str_to_lower(const t_str& str)
{
	typedef typename std::string::value_type t_char;

	t_str strOut;
	strOut.reserve(str.length());

	for(t_char ch : str)
		strOut.push_back(std::tolower(ch));

	return strOut;
}

template<class t_str=std::string>
t_str str_to_upper(const t_str& str)
{
	typedef typename std::string::value_type t_char;

	t_str strOut;
	strOut.reserve(str.length());

	for(t_char ch : str)
		strOut.push_back(std::toupper(ch));

	return strOut;
}

template<class t_str=std::string>
bool str_is_equal(const t_str& str0, const t_str& str1, bool bCase=0)
{
	if(str0.size() != str1.size())
		return false;

	if(bCase) return str0==str1;

	for(unsigned int i=0; i<str0.size(); ++i)
	{
		if(std::tolower(str0[i]) != std::tolower(str1[i]))
			return false;
	}
	return true;
}

template<class t_str=std::string>
bool str_is_equal_to_either(const t_str& str0,
	const std::initializer_list<t_str>& lststr1, bool bCase=0)
{
	for(const t_str& str1 : lststr1)
		if(str_is_equal<t_str>(str0, str1, bCase))
			return true;
	return false;
}

template<class t_str=std::string>
t_str str_tolower(const t_str& str)
{
	t_str strLower;
	strLower.reserve(str.length());

	for(typename t_str::value_type ch : str)
		strLower.push_back(std::tolower(ch));

	return strLower;
}

template<class t_str=std::string>
bool str_contains(const t_str& str, const t_str& strSub, bool bCase=0)
{
	if(bCase)
		return str.find(strSub) != t_str::npos;

	t_str strLower = str_tolower(str);
	t_str strSubLower = str_tolower(strSub);

	return strLower.find(strSubLower) != t_str::npos;
}


template<class t_str=std::string>
void trim(t_str& str)
{
	std::size_t posLast = str.find_last_not_of(get_trim_chars<t_str>());
	if(posLast == std::string::npos)
			posLast = str.length();
	else
		++posLast;

	str.erase(str.begin()+posLast, str.end());


	std::size_t posFirst = str.find_first_not_of(get_trim_chars<t_str>());
	if(posFirst == t_str::npos)
		posFirst = str.length();

	str.erase(str.begin(), str.begin()+posFirst);
}

template<class t_str=std::string>
t_str trimmed(const t_str& str)
{
	t_str strret = str;
	trim(strret);
	return strret;
}

template<class t_str=std::string>
t_str remove_char(const t_str& str, typename t_str::value_type ch)
{
	t_str strRet;

	for(typename t_str::value_type c : str)
		if(c != ch)
			strRet.push_back(c);

	return strRet;
}

/**
 * Removes substring between strStart and strEnd
 * @return Number of removed substrings
 */
template<class t_str = std::string>
unsigned int string_rm(t_str& str, const t_str& strStart, const t_str& strEnd)
{
	unsigned int iNumFound = 0;

	while(1)
	{
		std::size_t iStart = str.find(strStart);
		std::size_t iEnd = str.find(strEnd);

		if(iStart == t_str::npos || iEnd == t_str::npos)
			break;
		if(iStart >= iEnd)
			break;

		str.erase(iStart, iEnd-iStart+strEnd.length());
		++iNumFound;
	}

	return iNumFound;
}

template<class t_str=std::string>
bool find_and_replace(t_str& str1, const t_str& str_old,
						const t_str& str_new)
{
	std::size_t pos = str1.find(str_old);
	if(pos==t_str::npos)
			return false;

	str1.replace(pos, str_old.length(), str_new);
	return true;
}

template<class t_str=std::string>
void find_all_and_replace(t_str& str1, const t_str& str_old,
						const t_str& str_new)
{
	std::size_t pos=0;
	while(pos < str1.length())
	{
		pos = str1.find(str_old, pos);
		if(pos==t_str::npos)
			break;
		str1.replace(pos, str_old.length(), str_new);
		pos += str_new.length();
	}
}


template<class t_str=std::string>
bool begins_with(const t_str& str, const t_str& strBeg)
{
	if(str.length() < strBeg.length())
		return false;

	for(unsigned int i=0; i<strBeg.length(); ++i)
		if(str[i] != strBeg[i])
			return false;

	return true;
}


template<class t_str=std::string>
std::pair<t_str, t_str>
split_first(const t_str& str, const t_str& strSep, bool bTrim=0)
{
	t_str str1, str2;

	std::size_t pos = str.find_first_of(strSep);
	if(pos != t_str::npos)
	{
		str1 = str.substr(0, pos);
		if(pos+1 < str.length())
			str2 = str.substr(pos+1, t_str::npos);
	}

	if(bTrim)
	{
		trim(str1);
		trim(str2);
	}

	return std::pair<t_str, t_str>(str1, str2);
}


template<typename T, class t_str=std::string, bool bTIsStr=0>
struct _str_to_var_impl
{
	inline T operator()(const t_str&) const { throw Err("No implementation for str_to_var"); }
};

template<typename T, class t_str>
struct _str_to_var_impl<T, t_str, 1>
{
	inline const T& operator()(const t_str& str) const
	{
		return str;
	}
};

template<typename T, class t_str>
struct _str_to_var_impl<T, t_str, 0>
{
	inline T operator()(const t_str& str) const
	{
		typedef typename t_str::value_type t_char;
		std::basic_istringstream<t_char> istr(str);

		T t;
		istr >> t;
		return t;
	}
};


template<class T, class t_str=std::string, class t_cont=std::vector<T>>
void get_tokens(const t_str& str, const t_str& strDelim, t_cont& vecRet)
{
	typedef typename t_str::value_type t_char;

	boost::char_separator<t_char> delim(strDelim.c_str());
	boost::tokenizer<boost::char_separator<t_char> > tok(str, delim);

	typename boost::tokenizer<boost::char_separator<t_char> >::iterator iter;
	for(iter=tok.begin(); iter!=tok.end(); ++iter)
	{
		T t = _str_to_var_impl<T, t_str, !std::is_fundamental<T>::value>()(*iter);
		vecRet.push_back(std::move(t));
	}
}

template<typename T, class t_str=std::string>
T str_to_var(const t_str& str)
{
	return _str_to_var_impl<T, t_str, !std::is_fundamental<T>::value>()(str);
}

template<typename T, class t_str=std::string>
t_str var_to_str(const T& t, std::streamsize iPrec=10, int iGroup=-1)
{
	typedef typename t_str::value_type t_char;

	std::basic_ostringstream<t_char> ostr;
	ostr.precision(iPrec);

	class Sep : public std::numpunct<t_char>
	{
	public:
		Sep() : std::numpunct<t_char>(1) {}
		~Sep() { /*std::cout << "~Sep();" << std::endl;*/ }
	protected:
		virtual t_char do_thousands_sep() const override { return ' ';}
		virtual std::string do_grouping() const override { return "\3"; }
	};
	Sep *pSep = nullptr;

	if(iGroup > 0)
	{
		pSep = new Sep();
		ostr.imbue(std::locale(ostr.getloc(), pSep));
	}


	// prevents printing "-0"
	T t0 = t;
	if(t0==T(-0)) t0 = T(0);

	ostr << t0;
	t_str str = ostr.str();

	if(pSep)
	{
		ostr.imbue(std::locale());
		delete pSep;
	}
	return str;
}


template<typename t_char=char>
bool skip_after_line(std::basic_istream<t_char>& istr,
		const std::basic_string<t_char>& strLineBegin,
		bool bTrim=true, bool bCase=0)
{
	while(!istr.eof())
	{
		std::basic_string<t_char> strLine;
		std::getline(istr, strLine);
		if(bTrim)
			trim(strLine);

		if(strLine.size() < strLineBegin.size())
			continue;

		std::basic_string<t_char> strSub = strLine.substr(0, strLineBegin.size());

		if(str_is_equal<std::basic_string<t_char>>(strSub, strLineBegin, bCase))
			return true;
	}
	return false;
}

template<typename t_char=char>
void skip_after_char(std::basic_istream<t_char>& istr, t_char ch, bool bCase=0)
{
	if(!bCase) ch = std::tolower(ch);

	while(!istr.eof())
	{
		t_char c;
		istr.get(c);

		if(!bCase) c = std::tolower(c);

		if(c == ch)
			break;
	}
}


}
#endif
