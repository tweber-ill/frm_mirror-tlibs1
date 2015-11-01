/*
 * xml helper
 * @author tweber
 * @date 23-apr-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_XML__
#define __TLIBS_XML__

#include <string>
#include <sstream>
#include <map>
//#include <mutex>
#include <boost/property_tree/ptree.hpp>

namespace tl {

class Xml
{
protected:
	boost::property_tree::ptree m_xml;
	bool m_bOK = 0;
	//mutable std::mutex m_mutex;

public:
	Xml() = default;
	virtual ~Xml() = default;

	bool Load(const char* pcFile);

	std::string QueryString(const char* pcAdr, const char* pcDef, bool *pbOk=0);
	bool Exists(const char* pcAdr);

	template<typename T>
	T Query(const char* pcAdr, T tDef, bool *pbOk=0)
	{
		if(!m_bOK) return tDef;

		std::ostringstream ostrDef;
		ostrDef << tDef;

		bool bOk = 0;
		std::string str = QueryString(pcAdr, ostrDef.str().c_str(), &bOk);
		if(pbOk) *pbOk = bOk;

		if(bOk)
		{
			T tRes(0);

			std::istringstream istr(str);
			istr >> tRes;

			return tRes;
		}
		else
			return tDef;
	}


	static bool SaveMap(const char* pcFile, const std::map<std::string, std::string>& mapXml);
};

}

#endif
