/*
 * Simple Script
 * @author tweber
 */

#include "node.h"

template<typename T> T plus_op(T a, T b) { return a+b; }
template<typename T> T minus_op(T a, T b) { return a-b; }
template<typename T> T mult_op(T a, T b) { return a*b; }
template<typename T> T div_op(T a, T b) { return a/b; }

std::map<NodeType, double (*)(double,double)> g_mapBinOps_d = 
{
	std::pair<NodeType, double (*)(double,double)>(NODE_PLUS, plus_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_MINUS, minus_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_MULT, mult_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_DIV, div_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_POW, pow)
};

Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op)
{
	if(!pSymLeft || !pSymRight) return 0;
	if(pSymLeft->GetType() != pSymRight->GetType())
		; 	//TODO

	if(pSymLeft->GetType() == SYMBOL_DOUBLE)
	{
		SymbolDouble *pResult = new SymbolDouble();
		pResult->strName = "<op_result>";
		pResult->dVal = g_mapBinOps_d[op](((SymbolDouble*)pSymLeft)->dVal, 
						((SymbolDouble*)pSymRight)->dVal);

		return pResult;
	}
	
	return 0;
}

void safe_delete(Symbol *pSym, const SymbolTable* pSymTab)
{
	if(!pSym) return;
	
	bool bIsInTable = pSymTab->IsPtrInMap(pSym);
	if(!bIsInTable)
		delete pSym;
}
