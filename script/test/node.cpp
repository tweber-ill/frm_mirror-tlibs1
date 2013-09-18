/*
 * Simple Script
 * @author tweber
 */

#include "node.h"

template<typename T> T plus_op(T a, T b) { return a+b; }
template<typename T> T minus_op(T a, T b) { return a-b; }
template<typename T> T mult_op(T a, T b) { return a*b; }
template<typename T> T div_op(T a, T b) { return a/b; }

template<typename T> bool log_or_op(T a, T b) { return a||b; }
template<typename T> bool log_and_op(T a, T b) { return a&&b; }
template<typename T> bool log_eq_op(T a, T b) { return a==b; }
template<typename T> bool log_neq_op(T a, T b) { return a!=b; }
template<typename T> bool log_leq_op(T a, T b) { return a<=b; }
template<typename T> bool log_geq_op(T a, T b) { return a>=b; }
template<typename T> bool log_less_op(T a, T b) { return a<b; }
template<typename T> bool log_greater_op(T a, T b) { return a>b; }

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
	std::pair<NodeType, double (*)(double,double)>(NODE_POW, pow),
};

std::map<NodeType, int (*)(int,int)> g_mapBinOps_i =
{
	std::pair<NodeType, int (*)(int,int)>(NODE_PLUS, plus_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_MINUS, minus_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_MULT, mult_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_DIV, div_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_POW, pow_int)
};


std::map<NodeType, bool (*)(double,double)> g_mapBinLogOps_d =
{
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_OR, log_or_op),
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_AND, log_and_op),
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_EQ, log_eq_op),
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_NEQ, log_neq_op),
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_LEQ, log_leq_op),
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_GEQ, log_geq_op),
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_LESS, log_less_op),
	std::pair<NodeType, bool (*)(double,double)>(NODE_LOG_GREATER, log_greater_op)
};

std::map<NodeType, bool (*)(int,int)> g_mapBinLogOps_i =
{
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_OR, log_or_op),
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_AND, log_and_op),
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_EQ, log_eq_op),
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_NEQ, log_neq_op),
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_LEQ, log_leq_op),
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_GEQ, log_geq_op),
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_LESS, log_less_op),
	std::pair<NodeType, bool (*)(int,int)>(NODE_LOG_GREATER, log_greater_op)
};


static inline bool IsLogicalOp(NodeType op)
{
	switch(op)
	{
		case NODE_LOG_AND:
		case NODE_LOG_OR:
		case NODE_LOG_NOT:
		case NODE_LOG_EQ:
		case NODE_LOG_NEQ:
		case NODE_LOG_LESS:
		case NODE_LOG_GREATER:
		case NODE_LOG_LEQ:
		case NODE_LOG_GEQ:
			return 1;
	}

	return 0;
}

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
		Symbol *pRes = 0;
		if(IsLogicalOp(op))
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->strName = "<op_logical_result>";
			pResult->iVal = g_mapBinLogOps_d[op](((SymbolDouble*)pLeft)->dVal,
												((SymbolDouble*)pRight)->dVal);
			pRes = pResult;
		}
		else
		{
			SymbolDouble *pResult = new SymbolDouble();
			pResult->strName = "<op_result>";
			pResult->dVal = g_mapBinOps_d[op](((SymbolDouble*)pLeft)->dVal,
											((SymbolDouble*)pRight)->dVal);
			pRes = pResult;
		}

		return pRes;
	}
	else if(pLeft->GetType() == SYMBOL_INT)
	{
		Symbol *pRes = 0;
		if(IsLogicalOp(op))
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->strName = "<op_logical_result>";
			pResult->iVal = g_mapBinLogOps_i[op](((SymbolInt*)pLeft)->iVal,
												((SymbolInt*)pRight)->iVal);
			pRes = pResult;
		}
		else
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->strName = "<op_result>";
			pResult->iVal = g_mapBinOps_i[op](((SymbolInt*)pLeft)->iVal,
											((SymbolInt*)pRight)->iVal);
			pRes = pResult;
		}

		return pRes;
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


	std::vector<Symbol*> vecArgSyms;
	for(Node* pNode : vecArgs)
	{
		Symbol *pSymbol = pNode->eval(pSym, vecFuncs);
		//std::cout << "argument: " << pSymbol->print() << std::endl;

		vecArgSyms.push_back(pSymbol);
	}

	pFkt->SetArgSyms(&vecArgSyms);
	Symbol* pFktRet = pFkt->eval(pSym, vecFuncs);

	for(Symbol *pArgSym : vecArgSyms)
		safe_delete(pArgSym, pSym);
	return pFktRet;
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
			else if(pSymbol->GetType() == SYMBOL_INT)
				((SymbolInt*)pSymbol)->iVal = -((SymbolInt*)pSymbol)->iVal;

			return pSymbol;
		}

		case NODE_LOG_NOT:
		{
			Symbol *pSymbolEval = m_pChild->eval(pSym, vecFuncs);
			SymbolInt *pSymbolInt = new SymbolInt();

			if(pSymbolEval->GetType() == SYMBOL_DOUBLE)
				pSymbolInt->iVal = !((SymbolDouble*)pSymbolEval)->dVal;
			else if(pSymbolEval->GetType() == SYMBOL_INT)
				pSymbolInt->iVal = !((SymbolInt*)pSymbolEval)->iVal;

			safe_delete(pSymbolEval, pSym);
			return pSymbolInt;
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
	if(m_pVecArgSyms)
	{
		if(vecArgs.size() != m_pVecArgSyms->size())
		{
			std::cerr << "Error: Function \"" << strName << "\"" << " takes "
					 << vecArgs.size() << " arguments, but "
					 << m_pVecArgSyms->size() << " given."
					 << std::endl;
		}

		for(unsigned int iArg=0; iArg<vecArgs.size(); ++iArg)
		{
			Node* pNode = vecArgs[iArg];
			Symbol *pSymbol = (*m_pVecArgSyms)[iArg];

			NodeIdent* pIdent = (NodeIdent*)pNode;
			//std::cout << "arg: " << pIdent->m_strIdent << std::endl;

			/*Symbol *pSymbol = pSym->GetSymbol(pIdent->m_strIdent);
			if(!pSymbol)
			{
				std::cerr << "Error: Symbol \"" << pIdent->m_strIdent << "\" not found."
							<< std::endl;
			}*/

			pLocalSym->InsertSymbol(pIdent->m_strIdent, pSymbol->clone());
		}
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


Symbol* NodeIf::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	Symbol *pSymExpr = 0;
	Symbol *pSymRet = 0;
	if(m_pExpr)
		pSymExpr = m_pExpr->eval(pSym, vecFuncs);

	if(pSymExpr && pSymExpr->IsNotZero())
		pSymRet = (m_pIf ? m_pIf->eval(pSym, vecFuncs) : 0);
	else
		pSymRet = (m_pElse ? m_pElse->eval(pSym, vecFuncs) : 0);

	safe_delete(pSymExpr, pSym);
	return pSymRet;
}


Symbol* NodeWhile::eval(SymbolTable *pSym, std::vector<NodeFunction*>& vecFuncs) const
{
	if(!m_pExpr) return 0;
	if(!m_pStmt) return 0;

	Symbol *pSymRet = 0;

	while(1)
	{
		safe_delete(pSymRet, pSym);

		Symbol *pSymExpr = m_pExpr->eval(pSym, vecFuncs);

		if(pSymExpr && pSymExpr->IsNotZero())
			pSymRet = m_pStmt->eval(pSym, vecFuncs);

		safe_delete(pSymExpr, pSym);
	}
	return pSymRet;
}

//--------------------------------------------------------------------------------


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
