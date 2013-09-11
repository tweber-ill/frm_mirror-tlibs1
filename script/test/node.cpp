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


//--------------------------------------------------------------------------------


Symbol* NodeIdent::eval(SymbolTable *pSym) const
{
	return pSym->GetSymbol(m_strIdent);
}

Symbol* NodeCall::eval(SymbolTable *pSym) const
{
	if(m_pIdent->m_type != NODE_IDENT)
		return 0;
	//if(m_pArgs->m_type != NODE_ARGS)
	//	return 0;

	NodeIdent* pIdent = (NodeIdent*) m_pIdent;
	NodeBinaryOp* pArgs = (NodeBinaryOp*) m_pArgs;

	std::vector<Node*> vecArgs;
	if(pArgs)
		vecArgs = pArgs->flatten(NODE_ARGS);

	std::string strFkt = pIdent->m_strIdent;
	//std::cout << "call to " << strFkt << " with " << vecArgs.size() << " arguments." << std::endl;
	for(Node* pNode : vecArgs)
	{
		Symbol *pSymbol = pNode->eval(pSym);
		//std::cout << "argument: " << pSymbol->print() << std::endl;
	}

	return 0;
}

Symbol* NodeDouble::eval(SymbolTable *pSym) const
{
	SymbolDouble *pSymbol = new SymbolDouble;
	pSymbol->dVal = m_dVal;
	pSymbol->strName = "<const>";
	return pSymbol;
}

Symbol* NodeUnaryOp::eval(SymbolTable *pSym) const
{
	switch(m_type)
	{
		case NODE_UMINUS:
		{
			Symbol *pSymbolEval = m_pChild->eval(pSym);
			Symbol *pSymbol = pSymbolEval->clone();
			safe_delete(pSymbolEval, pSym);

			if(pSymbol->GetType() == SYMBOL_DOUBLE)
				((SymbolDouble*)pSymbol)->dVal = -((SymbolDouble*)pSymbol)->dVal;

			return pSymbol;
		}
	}
}

Symbol* NodeBinaryOp::eval(SymbolTable *pSym) const
{
	switch(m_type)
	{
		case NODE_STMTS:
		case NODE_ARGS:
			if(m_pLeft)
			{
				Symbol *pSymbol = m_pLeft->eval(pSym);
				safe_delete(pSymbol, pSym);
			}
			if(m_pRight)
			{
				Symbol *pSymbol = m_pRight->eval(pSym);
				safe_delete(pSymbol, pSym);
			}
			return 0;
		case NODE_ASSIGN:
			Symbol *pSymbol = m_pRight->eval(pSym);
			const std::string& strIdent = ((NodeIdent*)m_pLeft)->m_strIdent;

			pSym->InsertSymbol(strIdent, pSymbol);
			return pSymbol;
	};

	Symbol *pSymbolLeft = m_pLeft->eval(pSym);
	Symbol *pSymbolRight = m_pRight->eval(pSym);
	Symbol *pSymbol = Op(pSymbolLeft, pSymbolRight, m_type);
	safe_delete(pSymbolLeft, pSym);
	safe_delete(pSymbolRight, pSym);

	return pSymbol;
}

std::vector<Node*> NodeBinaryOp::flatten(NodeType ntype)
{
	//std::cout << m_type << ", " << ntype << std::endl;

	std::vector<Node*> vecNodes;

	NodeBinaryOp *pLeft = (NodeBinaryOp*) m_pLeft;
	NodeBinaryOp *pRight = (NodeBinaryOp*) m_pRight;

	if(m_type == ntype)
	{
		if(m_pLeft)
		{
			if(m_pLeft->m_type == ntype)
			{
				std::vector<Node*> vecLeft = pLeft->flatten(ntype);
				vecNodes.insert(vecNodes.begin(), vecLeft.begin(), vecLeft.end());
			}
			else
			{
				vecNodes.push_back(m_pLeft);
			}
		}

		if(m_pRight)
		{
			if(m_pRight->m_type == ntype)
			{
				std::vector<Node*> vecRight = pRight->flatten(ntype);
				vecNodes.insert(vecNodes.end(), vecRight.begin(), vecRight.end());
			}
			else
			{
				vecNodes.push_back(m_pRight);
			}
		}
	}
	else
	{
		vecNodes.push_back(this);
	}

	return vecNodes;
}
