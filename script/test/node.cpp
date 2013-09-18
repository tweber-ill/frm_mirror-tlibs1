/*
 * Simple Script
 * @author tweber
 */

#include "node.h"

template<typename T> T plus_op(T a, T b) { return a+b; }
template<typename T> T minus_op(T a, T b) { return a-b; }
template<typename T> T mult_op(T a, T b) { return a*b; }
template<typename T> T div_op(T a, T b) { return a/b; }

static int pow_int(int a, int b)
{
	return int(pow(a,b));
}

std::map<NodeType, double (*)(double,double)> g_mapBinOps_d = 
{
	std::pair<NodeType, double (*)(double,double)>(NODE_PLUS, plus_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_MINUS, minus_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_MULT, mult_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_DIV, div_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_POW, pow)
};

std::map<NodeType, int (*)(int,int)> g_mapBinOps_i =
{
	std::pair<NodeType, int (*)(int,int)>(NODE_PLUS, plus_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_MINUS, minus_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_MULT, mult_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_DIV, div_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_POW, pow_int)
};

Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op)
{
	if(!pSymLeft || !pSymRight) return 0;

	const Symbol *pLeft = pSymLeft;
	const Symbol *pRight = pSymRight;
	bool bCleanLeft = 0;
	bool bCleanRight = 0;

	if(pSymLeft->GetType() != pSymRight->GetType())
	{
		if(pLeft->GetType()==SYMBOL_INT && pRight->GetType()==SYMBOL_DOUBLE)
		{
			pLeft = pLeft->ToType(SYMBOL_DOUBLE);
			bool bCleanLeft = 1;
		}

		if(pLeft->GetType()==SYMBOL_DOUBLE && pRight->GetType()==SYMBOL_INT)
		{
			pRight = pRight->ToType(SYMBOL_DOUBLE);
			bool bCleanRight = 1;
		}
	}

	if(pLeft->GetType() == SYMBOL_DOUBLE)
	{
		SymbolDouble *pResult = new SymbolDouble();
		pResult->strName = "<op_result>";
		pResult->dVal = g_mapBinOps_d[op](((SymbolDouble*)pLeft)->dVal,
						((SymbolDouble*)pRight)->dVal);

		return pResult;
	}
	else if(pLeft->GetType() == SYMBOL_INT)
	{
		SymbolInt *pResult = new SymbolInt();
		pResult->strName = "<op_result>";
		pResult->iVal = g_mapBinOps_i[op](((SymbolInt*)pLeft)->iVal,
						((SymbolInt*)pRight)->iVal);

		return pResult;
	}
	
	if(bCleanLeft) delete pLeft;
	if(bCleanRight) delete pRight;

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


Symbol* NodeIdent::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	return pSym->GetSymbol(m_strIdent);
}

Symbol* NodeCall::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
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


	NodeFunction *pFkt = 0;
	for(NodeFunction *pFktIter : vecFuncs)
	{
		if(pFktIter && pFktIter->GetName()==strFkt)
			pFkt = pFktIter;
	}
	if(!pFkt)
	{
		std::cerr << "Error: Trying to call unknown function \" << strFkt << \"."
					<< std::endl;
		return 0;
	}


	for(Node* pNode : vecArgs)
	{
		Symbol *pSymbol = pNode->eval(pSym, vecFuncs);
		//std::cout << "argument: " << pSymbol->print() << std::endl;
	}


	return pFkt->eval(pSym, vecFuncs);
}

Symbol* NodeDouble::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	SymbolDouble *pSymbol = new SymbolDouble;
	pSymbol->dVal = m_dVal;
	pSymbol->strName = "<const>";
	return pSymbol;
}

Symbol* NodeInt::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	SymbolInt *pSymbol = new SymbolInt;
	pSymbol->iVal = m_iVal;
	pSymbol->strName = "<const>";
	return pSymbol;
}

Symbol* NodeString::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	SymbolString *pSymbol = new SymbolString;
	pSymbol->strVal = m_strVal;
	pSymbol->strName = "<const>";
	return pSymbol;
}

Symbol* NodeUnaryOp::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	switch(m_type)
	{
		case NODE_UMINUS:
		{
			Symbol *pSymbolEval = m_pChild->eval(pSym, vecFuncs);
			Symbol *pSymbol = pSymbolEval->clone();
			safe_delete(pSymbolEval, pSym);

			if(pSymbol->GetType() == SYMBOL_DOUBLE)
				((SymbolDouble*)pSymbol)->dVal = -((SymbolDouble*)pSymbol)->dVal;
			else if(pSymbol->GetType() == SYMBOL_DOUBLE)
				((SymbolInt*)pSymbol)->iVal = -((SymbolInt*)pSymbol)->iVal;

			return pSymbol;
		}

		case NODE_STMTS:
		{
			if(m_pChild)
			{
				Symbol *pSymbol = m_pChild->eval(pSym, vecFuncs);
				safe_delete(pSymbol, pSym);
			}
			return 0;
		}
	}
}

Symbol* NodeBinaryOp::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	switch(m_type)
	{
		case NODE_FUNCS:
		{
			std::vector<Node*> vecFuncs0;
			vecFuncs0 = this->flatten(NODE_FUNCS);
			for(Node *pNodeFunc : vecFuncs0)
				vecFuncs.push_back((NodeFunction*)pNodeFunc);

			for(NodeFunction* pFkt : vecFuncs)
			{
				if(pFkt->GetName() == "main")
					return pFkt->eval(pSym, vecFuncs);
			}

			std::cerr << "Error: No main function defined." << std::endl;
			return 0;
		}

		case NODE_STMTS:
		case NODE_ARGS:
		{
			if(m_pLeft)
			{
				//std::cout << "left: " << m_pLeft->m_type << std::endl;
				Symbol *pSymbol = m_pLeft->eval(pSym, vecFuncs);
				safe_delete(pSymbol, pSym);
			}
			if(m_pRight)
			{
				//std::cout << "right: " << m_pRight->m_type << std::endl;
				Symbol *pSymbol = m_pRight->eval(pSym, vecFuncs);
				safe_delete(pSymbol, pSym);
			}
			return 0;
		}

		case NODE_ASSIGN:
		{
			Symbol *pSymbol = m_pRight->eval(pSym, vecFuncs);
			const std::string& strIdent = ((NodeIdent*)m_pLeft)->m_strIdent;

			pSym->InsertSymbol(strIdent, pSymbol);
			return pSymbol;
		}
	};

	Symbol *pSymbolLeft = m_pLeft->eval(pSym, vecFuncs);
	Symbol *pSymbolRight = m_pRight->eval(pSym, vecFuncs);
	Symbol *pSymbol = Op(pSymbolLeft, pSymbolRight, m_type);
	safe_delete(pSymbolLeft, pSym);
	safe_delete(pSymbolRight, pSym);

	return pSymbol;
}


Symbol* NodeFunction::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	std::string strName = GetName();
	//std::cout << "in fkt " << strName << std::endl;

	std::vector<Node*> vecArgs;
	if(m_pArgs && (m_pArgs->m_type==NODE_IDENTS || m_pArgs->m_type==NODE_IDENT))
	{
		NodeBinaryOp* pArgs = (NodeBinaryOp*) m_pArgs;
		vecArgs = pArgs->flatten(NODE_IDENTS);
	}


	SymbolTable *pLocalSym = new SymbolTable;
	for(Node* pNode : vecArgs)
	{
		NodeIdent* pIdent = (NodeIdent*)pNode;
		//std::cout << "arg: " << pIdent->m_strIdent << std::endl;

		Symbol *pSymbol = pSym->GetSymbol(pIdent->m_strIdent);
		if(!pSymbol)
		{
			std::cerr << "Error: Symbol \"" << pIdent->m_strIdent << "\" not found."
						<< std::endl;
		}

		pLocalSym->InsertSymbol(pIdent->m_strIdent, pSymbol->clone());
	}


	Symbol *pRet = 0;
	if(m_pStmts)
	{
		pRet = m_pStmts->eval(pLocalSym, vecFuncs);
		if(!pRet)
			pRet = pLocalSym->GetSymbol("_ret_");
	}

	std::cout << "Local symbols for \"" << strName << "\":\n";
	pLocalSym->print();

	delete pLocalSym;
	return pRet;
}




std::vector<Node*> NodeBinaryOp::flatten(NodeType ntype) const
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
		vecNodes.push_back(const_cast<NodeBinaryOp*>(this));
	}

	return vecNodes;
}


std::string NodeFunction::GetName() const
{
	if(!m_pIdent || m_pIdent->m_type != NODE_IDENT)
		return "<null>";

	NodeIdent *pIdent = (NodeIdent*)m_pIdent;
	return pIdent->m_strIdent;
}
