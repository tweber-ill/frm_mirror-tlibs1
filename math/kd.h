/*
 * kd trees
 * @author tweber
 * @date jun-2015
 * @license GPLv2 or GPLv3
 */
 
#ifndef __TLIBS_KD_H__
#define __TLIBS_KD_H__

#include <vector>
#include <list>
#include <algorithm>
#include <iostream>


template<class T=double>
struct KdNode
{
	std::vector<T> vecMid;
	KdNode<T> *pLeft = nullptr;
	KdNode<T> *pRight = nullptr;
	
	void print(std::ostream& ostr, unsigned int iLevel=0) const
	{
		for(unsigned int i=0; i<iLevel; ++i) ostr << "\t";
		ostr << "mid: ";
		for(T t : vecMid)
			ostr << t << ", ";

		if(pLeft) 
		{
			ostr << "\n";
			for(unsigned int i=0; i<iLevel; ++i) ostr << "\t";
			ostr << "left:\n";
			pLeft->print(ostr, iLevel+1);
		}
		if(pRight)
		{
			ostr << "\n";
			for(unsigned int i=0; i<iLevel; ++i) ostr << "\t";
			ostr << "right:\n";
			pRight->print(ostr, iLevel+1);
		}
		
		ostr << "\n";
	}
};


template<class T=double>
class Kd
{
private:
	static KdNode<T>* make_kd(std::list<std::vector<T>>& lstPoints, unsigned int iLevel=0)
	{
		const unsigned int iSize = lstPoints.size();
		if(iSize == 0) return nullptr;
		
		const unsigned int iDim = lstPoints.begin()->size();
		const unsigned int iAxis = iLevel % iDim;
		
		KdNode<T> *pNode = new KdNode<T>;
		
		if(iSize == 1)
		{
			pNode->vecMid = *lstPoints.begin();
			return pNode;
		}
		
		//std::cout << "level: " << iLevel << ", size: " << iSize << ", mid: " << iSize/2 << std::endl;
		
		lstPoints.sort(
			[iAxis](const std::vector<T>& vec0, const std::vector<T>& vec1) -> bool
			{
				return vec0[iAxis] <= vec1[iAxis];
			});

		typename std::list<std::vector<T>>::iterator iterMid = std::next(lstPoints.begin(), iSize/2);
		pNode->vecMid = *iterMid;
		std::list<std::vector<T>> lstLeft(lstPoints.begin(), iterMid);
		std::list<std::vector<T>> lstRight(std::next(iterMid), lstPoints.end());

		pNode->pLeft = make_kd(lstLeft, iLevel+1);
		pNode->pRight = make_kd(lstRight, iLevel+1);

		return pNode;
	}

	static void clear_kd(KdNode<T> *pNode)
	{
		if(!pNode) return;
		
		if(pNode->pLeft) clear_kd(pNode->pLeft);
		if(pNode->pRight) clear_kd(pNode->pRight);
		delete pNode;
	}

protected:
	KdNode<T> *m_pNode = nullptr;
	
public:
	void Unload()
	{
		clear_kd(m_pNode);
		m_pNode = nullptr;
	}

	// alters lstPoints!
	void Load(std::list<std::vector<T>>& lstPoints)
	{
		Unload();
		m_pNode = make_kd(lstPoints);
	}

	Kd() = default;
	Kd(std::list<std::vector<T>>& lstPoints) { Load(lstPoints); }
	virtual ~Kd() { Unload(); }
	
	const KdNode<T>* GetRootNode() const { return m_pNode; }
};

#endif
