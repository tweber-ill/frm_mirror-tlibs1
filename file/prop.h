/*
 * property tree wrapper
 * @author tweber
 * @date dec-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __PROP_FILES_H__
#define __PROP_FILES_H__

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <type_traits>
#include <map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "../string/string.h"
#include "../helper/traits.h"
#include "../helper/misc.h"


namespace tl {

namespace prop = ::boost::property_tree;

enum class PropType
{
	XML,
	JSON,
	INFO,
	INI,
};

template<class _t_str = std::string>
class Prop
{
public:
	using t_str = _t_str;
	using t_ch = typename t_str::value_type;

protected:
	prop::basic_ptree<t_str, t_str> m_prop;
	t_ch m_chSep = '/';

public:
	Prop() = default;
	virtual ~Prop() = default;

	void SetSeparator(t_ch ch) { m_chSep = ch; }


	bool Load(const t_ch* pcFile)
	{
		t_str strFile = pcFile;
		t_str strExt = str_to_lower<t_str>(get_fileext<t_str>(strFile));

		bool bOk = 0;
		if(strExt == "xml")
			bOk = Load(pcFile, PropType::XML);
		else if(strExt == "json")
			bOk = Load(pcFile, PropType::JSON);
		else if(strExt == "info")
			bOk = Load(pcFile, PropType::INFO);
		else if(strExt == "ini")
			bOk = Load(pcFile, PropType::INI);

		return bOk;
	}

	bool Load(const t_ch* pcFile, PropType ty)
	{
		std::basic_ifstream<t_ch> ifstr(pcFile);
		if(!ifstr) return false;

		return Load(ifstr, ty);
	}

	bool Load(std::basic_istream<t_ch>& istr, PropType ty)
	{
		try
		{
			switch(ty)
			{
				case PropType::XML: 
					prop::read_xml(istr, m_prop);
					break;
				case PropType::JSON: 
					prop::read_json(istr, m_prop);
					break;
				case PropType::INFO: 
					prop::read_info(istr, m_prop);
					break;
				case PropType::INI: 
					prop::read_ini(istr, m_prop);
					break;
				default:
					return false;
			}
		}
		catch(const prop::file_parser_error& err)
		{
			return false;
		}

		return true;
	}


	bool Save(const t_ch* pcFile) const
	{
		t_str strFile = pcFile;
		t_str strExt = str_to_lower<t_str>(get_fileext<t_str>(strFile));

		bool bOk = 0;
		if(strExt == "xml")
			bOk = Save(pcFile, PropType::XML);
		else if(strExt == "json")
			bOk = Save(pcFile, PropType::JSON);
		else if(strExt == "info")
			bOk = Save(pcFile, PropType::INFO);
		else if(strExt == "ini")
			bOk = Save(pcFile, PropType::INI);

		return bOk;
	}

	bool Save(const t_ch* pcFile, PropType ty) const
	{
		std::basic_ofstream<t_ch> ofstr(pcFile);
		if(!ofstr) return false;

		Save(ofstr, ty);
		return true;
	}

	bool Save(std::basic_ofstream<t_ch>& ofstr, PropType ty) const
	{
		#if BOOST_VERSION >= 105700
			using t_writer = t_str;
		#else
			using t_writer = t_ch;
		#endif

		try
		{
			switch(ty)
			{
				case PropType::XML: 
					prop::write_xml(ofstr, m_prop, prop::xml_writer_settings<t_writer>('\t',1));
					break;
				case PropType::JSON: 
					prop::write_json(ofstr, m_prop);
					break;
				case PropType::INFO: 
					prop::write_info(ofstr, m_prop);
					break;
				case PropType::INI: 
					prop::write_ini(ofstr, m_prop);
					break;
				default:
					return false;
			}
		}
		catch(const prop::file_parser_error& err)
		{
			return false;
		}

		return true;
	}


	template<typename T>
	T Query(const t_str& _strAddr, const T* pDef=nullptr, bool *pbOk=nullptr) const
	{
		t_str strAddr = _strAddr;
		trim(strAddr);

		if(strAddr.length() == 0)
		{
			if(pbOk) *pbOk = 0;
			return T();
		}

		if(strAddr[0] == m_chSep)
			strAddr = strAddr.substr(1);

		T tOut;
		try
		{
			prop::string_path<t_str, prop::id_translator<t_str>> path(strAddr, m_chSep);
			//tOut = m_prop.template get<T>(path);
			tOut = tl::str_to_var<T, t_str>(m_prop.template get<t_str>(path));
		}
		catch(const prop::ptree_bad_path& ex)
		{
			if(pbOk) *pbOk = 0;
			if(pDef) return *pDef;
			return T();
		}

		if(std::is_same<t_str, T>::value)
			trim(*reinterpret_cast<t_str*>(&tOut));

		if(pbOk) *pbOk = 1;
		return tOut;
	}

	template<typename T>
	T Query(const t_str& _strAddr, const T def, bool *pbOk=nullptr) const
	{
		return Query<T>(_strAddr, &def, pbOk);
	}


	bool Exists(const t_str& strAddr) const
	{
		bool bOk = 0;
		t_str strQuery = Query<t_str>(strAddr, nullptr, &bOk);
		if(strQuery.length() == 0)
			bOk = 0;

		return bOk;
	}


	/*template<class T = t_str>
	void Add(T&& tKey, T&& tVal)
	{
		prop::string_path<t_str, prop::id_translator<t_str>>
			path(std::forward<t_str>(tKey), m_chSep);
		m_prop.add(std::move(path), std::forward<t_str>(tVal));
	}*/

	template<class T>
	void Add(const t_str& strKey, T&& tVal)
	{
		prop::string_path<t_str, prop::id_translator<t_str>>
			path(strKey, m_chSep);

		//std::cout << "type: " << get_typename<remove_constref_t<T>>() << std::endl;
		if(std::is_convertible<T, t_str>::value)
		{
			//std::cout << "string: " << tVal << std::endl;
			m_prop.add(std::move(path), std::forward<T>(tVal));
		}
		else
		{
			t_str strVal = var_to_str<T, t_str>(tVal);
			m_prop.add(std::move(path), std::move(strVal));
		}
	}

	template<class t_map = std::map<t_str, t_str>>
	void Add(const t_map& map)
	{
		using t_key = typename t_map::key_type;
		using t_val = typename t_map::mapped_type;
		using t_pair = typename t_map::value_type;

		for(const t_pair& pair : map)
		{
			t_str strKey = var_to_str<t_key, t_str>(pair.first);
			t_str strVal = var_to_str<t_val, t_str>(pair.second);

			Add<t_str>(std::move(strKey), std::move(strVal));
		}
	}
};

}
#endif
