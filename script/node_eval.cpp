/*
 * Script interpreter
 * Node evaluation
 * @author tweber
 * @date 10 oct 2013
 */

#include "node.h"
#include "calls.h"


void safe_delete(Symbol *pSym, const SymbolTable* pSymTab)
{
	if(!pSym) return;

	// don't delete constants
	if(pSym->m_strName == "<const>")
		return;

	// don't delete symbols in table
	bool bIsInTable = pSymTab->IsPtrInMap(pSym);
	if(!bIsInTable)
		delete pSym;
}


Symbol* NodeReturn::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	Symbol *pRet = 0;

	if(m_pExpr)
		pRet = m_pExpr->eval(info, pSym)->clone();
	pSym->InsertSymbol("<ret>", pRet ? pRet : 0);

	info.bWantReturn = 1;
	return pRet;
}

Symbol* NodeIdent::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	Symbol *pSymbol = pSym->GetSymbol(m_strIdent);
	if(!pSymbol)
	{
		std::cerr << "Error: Symbol \"" << m_strIdent << "\" not in symbol table." << std::endl;
		return 0;
	}

	pSymbol->m_strIdent = m_strIdent;

	return pSymbol;
}

Symbol* NodeCall::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	if(m_pIdent->m_type != NODE_IDENT)
		return 0;
	//if(m_pArgs->m_type != NODE_ARGS)
	//	return 0;

	NodeIdent* pIdent = (NodeIdent*) m_pIdent;
	std::string strFkt = pIdent->m_strIdent;
	//std::cout << "call to " << strFkt << " with " << m_vecArgs.size() << " arguments." << std::endl;


	bool bCallUserFkt = 0;
	// user-defined function
	NodeFunction *pFkt = 0;
	for(NodeFunction *pFktIter : vecFuncs)
	{
		if(pFktIter && pFktIter->GetName()==strFkt)
			pFkt = pFktIter;
	}

	if(pFkt)
		bCallUserFkt = 1;

	/*if(!bCallUserFkt)
	{
		std::cerr << "Error: Trying to call unknown function \" << strFkt << \"."
					<< std::endl;
		return 0;
	}*/


	std::vector<Symbol*> vecArgSyms;
	for(Node* pNode : m_vecArgs)
	{
		Symbol *pSymbol = pNode->eval(info, pSym);
		//std::cout << "argument: " << pSymbol->print() << std::endl;

		vecArgSyms.push_back(pSymbol);
	}

	Symbol* pFktRet = 0;
	if(bCallUserFkt)	// call user-defined function
	{
		pFkt->SetArgSyms(&vecArgSyms);
		pFktRet = pFkt->eval(info, pSym);
		if(info.bWantReturn)
		{
			//std::cout << "returned" << std::endl;
			info.bWantReturn = 0;
		}
	}
	else				// call system function
	{
		pFktRet = ext_call(strFkt, vecArgSyms, info, pSym);
	}

	for(Symbol *pArgSym : vecArgSyms)
		safe_delete(pArgSym, pSym);
	return pFktRet;
}

Symbol* NodeDouble::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;
	return m_pSymbol;
}

Symbol* NodeInt::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;
	return m_pSymbol;
}

Symbol* NodeString::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;
	return m_pSymbol;
}


Symbol* NodeArray::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	NodeBinaryOp *pArr = (NodeBinaryOp*)m_pArr;
	std::vector<Node*> vecNodes = pArr->flatten(NODE_ARGS);

	SymbolArray *pSymArr = new SymbolArray;
	pSymArr->m_arr.reserve(vecNodes.size());

	for(Node* pNode : vecNodes)
	{
		if(!pNode) continue;

		Symbol *pSymbol = pNode->eval(info, pSym);
		pSymArr->m_arr.push_back(pSymbol);
	}

	return pSymArr;
}

Symbol* NodeArrayAccess::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(!m_pIdent || m_pIdent->m_type != NODE_IDENT)
	{
		std::cerr << "Error: Tried to access non-array." << std::endl;
		return 0;
	}

	std::string strIdent = ((NodeIdent*)m_pIdent)->m_strIdent;
	Symbol *pSymbol = pSym->GetSymbol(strIdent);

	if(!pSymbol)
	{
		std::cerr << "Error: Symbol \"" << strIdent
				<< "\" not found." << std::endl;
		return 0;
	}

	if(pSymbol->GetType() != SYMBOL_ARRAY)
	{
		std::cerr << "Error: Symbol \"" << strIdent
				<< "\" is no array." << std::endl;
		return 0;
	}



	for(Node *pIndices : m_vecIndices)
	{
		Symbol *pSymExpr = pIndices->eval(info, pSym);
		if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_INT)
		{
			std::cerr << "Error: Array index has to be of integer type."
						<< std::endl;
			return 0;
		}

		int iIdx = ((SymbolInt*)pSymExpr)->m_iVal;
		safe_delete(pSymExpr, pSym);

		if(iIdx<0 || iIdx>=((SymbolArray*)pSymbol)->m_arr.size())
		{
			std::cerr << "Error: Array index (" << iIdx
						<< ") out of bounds (array size: "
						<< ((SymbolArray*)pSymbol)->m_arr.size() << ")." << std::endl;
			return 0;
		}

		if(pSymbol->GetType() != SYMBOL_ARRAY)
		{
			std::cerr << "Error: Cannot take index of non-array." << std::endl;
			return 0;
		}

		pSymbol = ((SymbolArray*)pSymbol)->m_arr[iIdx];
	}

	return pSymbol;
}

Symbol* NodeUnaryOp::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	switch(m_type)
	{
		case NODE_UMINUS:
		{
			Symbol *pSymbolEval = m_pChild->eval(info, pSym);
			Symbol *pSymbol = pSymbolEval->clone();
			safe_delete(pSymbolEval, pSym);

			if(pSymbol->GetType() == SYMBOL_DOUBLE)
				((SymbolDouble*)pSymbol)->m_dVal = -((SymbolDouble*)pSymbol)->m_dVal;
			else if(pSymbol->GetType() == SYMBOL_INT)
				((SymbolInt*)pSymbol)->m_iVal = -((SymbolInt*)pSymbol)->m_iVal;

			return pSymbol;
		}

		case NODE_LOG_NOT:
		{
			Symbol *pSymbolEval = m_pChild->eval(info, pSym);
			SymbolInt *pSymbolInt = new SymbolInt();

			if(pSymbolEval->GetType() == SYMBOL_DOUBLE)
				pSymbolInt->m_iVal = !((SymbolDouble*)pSymbolEval)->m_dVal;
			else if(pSymbolEval->GetType() == SYMBOL_INT)
				pSymbolInt->m_iVal = !((SymbolInt*)pSymbolEval)->m_iVal;

			safe_delete(pSymbolEval, pSym);
			return pSymbolInt;
		}

		case NODE_STMTS:
		{
			if(m_pChild)
			{
				Symbol *pSymbol = m_pChild->eval(info, pSym);
				safe_delete(pSymbol, pSym);
			}
			return 0;
		}
	}
}

Symbol* NodeBinaryOp::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

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
					return pFkt->eval(info, pSym);
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
				Symbol *pSymbol = m_pLeft->eval(info, pSym);
				safe_delete(pSymbol, pSym);
			}
			if(m_pRight)
			{
				//std::cout << "right: " << m_pRight->m_type << std::endl;
				Symbol *pSymbol = m_pRight->eval(info, pSym);
				safe_delete(pSymbol, pSym);
			}
			return 0;
		}

		case NODE_ASSIGN:
		{
			Symbol *pSymbolOrg = m_pRight->eval(info, pSym);
			Symbol *pSymbol = pSymbolOrg->clone();
			safe_delete(pSymbolOrg, pSym);

			if(m_pLeft->m_type == NODE_IDENT)		// single variable
			{
				const std::string& strIdent = ((NodeIdent*)m_pLeft)->m_strIdent;
				pSym->InsertSymbol(strIdent, pSymbol);
			}
			else									// array
			{
				Symbol *pSymLeft = m_pLeft->eval(info, pSym);
				pSymLeft->assign(pSymbol);
			}

			return pSymbol;
		}
	};

	Symbol *pSymbolLeft = m_pLeft->eval(info, pSym);
	Symbol *pSymbolRight = m_pRight->eval(info, pSym);
	Symbol *pSymbol = Op(pSymbolLeft, pSymbolRight, m_type);
	safe_delete(pSymbolLeft, pSym);
	safe_delete(pSymbolRight, pSym);

	return pSymbol;
}


Symbol* NodeFunction::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;
	std::string strName = GetName();
	//std::cout << "in fkt " << strName << std::endl;


	SymbolTable *pLocalSym = new SymbolTable;
	if(m_pVecArgSyms)
	{
		if(m_vecArgs.size() != m_pVecArgSyms->size())
		{
			std::cerr << "Error: Function \"" << strName << "\"" << " takes "
					 << m_vecArgs.size() << " arguments, but "
					 << m_pVecArgSyms->size() << " given."
					 << std::endl;
		}

		for(unsigned int iArg=0; iArg<m_vecArgs.size(); ++iArg)
		{
			Node* pNode = m_vecArgs[iArg];
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
		pRet = m_pStmts->eval(info, pLocalSym);
		if(!pRet)
			pRet = pLocalSym->GetSymbol("<ret>");
	}

	//std::cout << "Local symbols for \"" << strName << "\":\n";
	//pLocalSym->print();

	delete pLocalSym;
	return pRet;
}


Symbol* NodeIf::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	Symbol *pSymExpr = 0;
	Symbol *pSymRet = 0;
	if(m_pExpr)
		pSymExpr = m_pExpr->eval(info, pSym);

	if(pSymExpr && pSymExpr->IsNotZero())
		pSymRet = (m_pIf ? m_pIf->eval(info, pSym) : 0);
	else
		pSymRet = (m_pElse ? m_pElse->eval(info, pSym) : 0);

	safe_delete(pSymExpr, pSym);
	safe_delete(pSymRet, pSym);

	return 0;
}


Symbol* NodeWhile::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	if(!m_pExpr) return 0;
	if(!m_pStmt) return 0;

	while(1)
	{
		Symbol *pSymRet = 0;
		Symbol *pSymExpr = m_pExpr->eval(info, pSym);

		if(pSymExpr && pSymExpr->IsNotZero())
			pSymRet = m_pStmt->eval(info, pSym);
		else
			break;

		safe_delete(pSymRet, pSym);
		safe_delete(pSymExpr, pSym);
	}

	return 0;
}

Symbol* NodeRangedFor::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.bWantReturn) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;
	if(!m_pIdent || !m_pExpr || !m_pStmt) return 0;

	if(m_pIdent->m_type != NODE_IDENT)
	{
		std::cerr << "Error: Range-based for loop needs identifier."
					<< std::endl;
		return 0;
	}

	Symbol *pSymRet = 0;
	Symbol *_pArr = m_pExpr->eval(info, pSym);
	if(_pArr->GetType() != SYMBOL_ARRAY)
	{
		std::cerr << "Error: Range-based for loop needs array." << std::endl;
		safe_delete(_pArr, pSym);
		return 0;
	}

	SymbolArray *pArr = (SymbolArray*)_pArr;


	const std::string& strIdent = ((NodeIdent*)m_pIdent)->m_strIdent;

	SymbolInt *pSymIter = new SymbolInt(0);
	std::string strIter = "<cur_iter_" + strIdent + ">";
	pSym->InsertSymbol(strIter, pSymIter);

	for(Symbol *pSymInArr : pArr->m_arr)
	{
		pSym->InsertSymbol(strIdent, pSymInArr);

		Symbol *pBodyRet = m_pStmt->eval(info, pSym);
		safe_delete(pBodyRet, pSym);

		++pSymIter->m_iVal;
	}

	pSym->RemoveSymbol(strIter);
	delete pSymIter;

	return 0;
}
