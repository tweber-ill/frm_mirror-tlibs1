/*
 * kd trees
 * @author tweber
 * @date jun-2015
 * @copyright GPLv2 or GPLv3
 */
 
// TODO

#ifndef __TLIBS_KD_H__
#define __TLIBS_KD_H__

#include <vector>
#include <list>

template<class T=double>
class Kd
{
protected:
	unsigned int m_iDim = 3;

public:
	void Unload()
	{
		
	}

	void Load(const std::list<std::vector<T>>& lstPoints)
	{
		Unload();
	}

	Kd() = default;
	Kd(const std::list<std::vector<T>>& lstPoints) { Load(lstPoints); }
	virtual ~Kd() { Unload(); }
};

#endif
