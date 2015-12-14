/*
 * xml helper
 * @author tweber
 * @date 23-apr-2013
 * @license GPLv2 or GPLv3
 */

#include "xml.h"
#include "../string/string.h"
#include "../helper/misc.h"
#include <fstream>
#include <list>
#include <vector>
#include <boost/property_tree/xml_parser.hpp>


namespace tl {

// -------------------------------------------------------------------------------
// loading

namespace prop = boost::property_tree;


bool Xml::Load(const char* pcFile)
{
	try
	{
		prop::read_xml(pcFile, m_xml);
	}
	catch(const prop::xml_parser_error&)
	{
		m_bOK = 0;
		return false;
	}

	m_bOK = 1;
	return true;
}

std::string Xml::QueryString(const char* pcAdr, const char* pcDef, bool *pbOk) const
{
	if(!m_bOK) return pcDef;

	std::string strAdr(pcAdr);
	trim(strAdr);
	if(strAdr.length()==0)
	{
		if(pbOk) *pbOk = 0;
		return "";
	}

	if(strAdr[0] == '/')
		strAdr = strAdr.substr(1);

	prop::string_path<std::string, prop::id_translator<std::string>> path(strAdr, '/');
	std::string strOut;
	try
	{
		//std::lock_guard<std::mutex> _lck(m_mutex);
		strOut = m_xml.get<std::string>(path);
	}
	catch(const prop::ptree_bad_path&)
	{
		if(pbOk) *pbOk = 0;
		if(pcDef) return pcDef;
		return "";
	}

	trim(strOut);
	if(pbOk) *pbOk = 1;
	return strOut;
}

bool Xml::Exists(const char* pcAdr) const
{
	bool bOk = 0;
	std::string strQuery = QueryString(pcAdr, "", &bOk);
	if(strQuery.length() == 0)
		bOk = 0;

	return bOk;
}



// -------------------------------------------------------------------------------
// saving

struct XmlNode
{
	std::string strName, strVal;
	std::vector<XmlNode*> vecChildren;

	XmlNode() {}

	~XmlNode()
	{
		for(XmlNode* pNode : vecChildren)
			delete pNode;
	}

	XmlNode* GetChild(const std::string& strKey)
	{
		for(XmlNode* pNode : vecChildren)
		{
			if(pNode->strName == strKey)
				return pNode;
		}

		XmlNode *pNode = new XmlNode;
		vecChildren.push_back(pNode);
		pNode->strName = strKey;
		return pNode;
	}

	void SortIntoTree(const std::list<std::string>& lstKey, const std::string& strVal)
	{
		if(lstKey.empty())
			return;

		XmlNode* pNode = GetChild(*lstKey.begin());
		if(lstKey.size()==1)
		{
			pNode->strVal = strVal;
		}
		else
		{
			std::list<std::string> lstKeyNext = lstKey;
			lstKeyNext.pop_front();

			pNode->SortIntoTree(lstKeyNext, strVal);
		}
	}

	void Print(std::ostream& ostr, unsigned int iLevel=1) const
	{
		const XmlNode& node = *this;

		ostr << "<" << node.strName << ">";
		if(node.vecChildren.size()==0)
			ostr << " " << node.strVal << " ";
		else
		{
			ostr << "\n";
			for(const XmlNode* pNode : node.vecChildren)
			{
				for(unsigned int iLev=0; iLev<iLevel; ++iLev)
					ostr << "\t";
				pNode->Print(ostr, iLevel+1);
				ostr << "\n";
			}
		}

		if(node.vecChildren.size()!=0)
			for(unsigned int iLev=0; iLev<iLevel-1; ++iLev)
				ostr << "\t";
		ostr << "</" << node.strName << ">";
	}
};

std::ostream& operator<<(std::ostream& ostr, const XmlNode& node)
{
	if(node.vecChildren.size() == 1)
		node.vecChildren[0]->Print(ostr);
 	return ostr;
}

bool Xml::SaveMap(const char* pcFile, const std::map<std::string, std::string>& mapXml)
{
	std::ofstream ofstr(pcFile);
	if(!ofstr.is_open())
		return false;

	XmlNode node;
	node.strName = "root";
	for(const std::pair<std::string, std::string>& pairXml : mapXml)
	{
		const std::string& strKey = pairXml.first;
		const std::string& strVal = pairXml.second;
		//std::cout << "key: " << strKey << ", val: " << strVal << std::endl;

		std::vector<std::string> vecKey;
		get_tokens<std::string>(strKey, std::string("/"), vecKey);

		std::list<std::string> lstKey = vector_to_list(vecKey);
		node.SortIntoTree(lstKey, strVal);
	}

	ofstr << node;
	return true;
}

}
