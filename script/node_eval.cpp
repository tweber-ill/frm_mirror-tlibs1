/*
 * Script interpreter
 * Node evaluation
 * @author tweber
 * @date 10 oct 2013
 */

#include "node.h"
#include "calls.h"


void safe_delete(Symbol *&pSym, const SymbolTable* pSymTab, const SymbolTable* pSymTabGlob)
{
	if(!pSym) return;

	// don't delete constants
	if(pSym->m_strName == "<const>")
		return;

	// don't delete array or map members
	if(pSym->m_pArr || pSym->m_pMap)
		return;

	// don't delete symbols in table
	bool bIsInTable = pSymTab->IsPtrInMap(pSym);
	bool bIsInGlobTable = pSymTabGlob->IsPtrInMap(pSym);
	if(!bIsInTable && !bIsInGlobTable)
	{
		//std::cout << "deleting " << (void*)pSym << ": " << pSym->GetType() << std::endl;
		delete pSym;
		pSym = 0;
	}
}


Symbol* NodeReturn::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	Symbol *pRet = 0;

	if(m_pExpr)
	{
		Symbol *pEval = m_pExpr->eval(info, pSym);
		pRet = pEval->clone();

		safe_delete(pEval, pSym, info.pGlobalSyms);
	}
	pSym->InsertSymbol("<ret>", pRet ? pRet : 0);

	info.bWantReturn = 1;
	return pRet;
}

Symbol* NodeBreak::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	if(info.pCurLoop==0)
	{
		std::cerr << linenr("Error", info)
				<< "Cannot use break outside loop."
				<< std::endl;
		return 0;
	}

	info.bWantBreak = 1;
	return 0;
}

Symbol* NodeContinue::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	if(info.pCurLoop==0)
	{
		std::cerr << linenr("Error", info)
				<< "Cannot use continue outside loop."
				<< std::endl;
		return 0;
	}

	info.bWantContinue = 1;
	return 0;
}

Symbol* NodeIdent::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	// local symbol
	Symbol *pSymbol = pSym->GetSymbol(m_strIdent);
	// global symbol
	Symbol *pSymbolGlob = info.pGlobalSyms->GetSymbol(m_strIdent);

	if(pSymbol && pSymbolGlob)
	{
		std::cerr << linenr("Warning", info) << "Symbol \"" << m_strIdent
				<< "\" exists in local and global scope, using local one." << std::endl;
	}

	// if no local symbol is available, use global symbol instead
	if(!pSymbol)
		pSymbol = pSymbolGlob;

	if(!pSymbol)
	{
		std::cerr << linenr("Error", info) << "Symbol \"" << m_strIdent
				<< "\" not in symbol table." << std::endl;
		return 0;
	}

	if(pSymbol->m_strIdent.size()==0)
		pSymbol->m_strIdent = m_strIdent;
	return pSymbol;
}

Symbol* NodeCall::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	info.pCurCaller = this;

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
		safe_delete(pArgSym, pSym, info.pGlobalSyms);
	return pFktRet;
}

Symbol* NodeDouble::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	return m_pSymbol;
}

Symbol* NodeInt::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	return m_pSymbol;
}

Symbol* NodeString::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	return m_pSymbol;
}

Symbol* NodePair::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::cerr << linenr("Error", info) 
		<< "Pairs should not be evaluated directly." 
		<< std::endl;
	return 0;
}

Symbol* NodeMap::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	NodeBinaryOp *pMap = (NodeBinaryOp*)m_pMap;
	std::vector<Node*> vecNodes;
	if(pMap) vecNodes = pMap->flatten(NODE_ARGS);

	SymbolMap *pSymMap = new SymbolMap;

	for(Node* pNode : vecNodes)
	{
		if(!pNode) continue;
		if(pNode->m_type != NODE_PAIR)
		{
			std::cerr << linenr("Error", info) 
				<< "Maps have to consist of key-value pairs." 
				<< std::endl;
			continue;
		}

		Symbol* pSymFirst = 0;
		Symbol* pSymSecond = 0;

		if(((NodePair*)pNode)->m_pFirst)
			pSymFirst = ((NodePair*)pNode)->m_pFirst->eval(info, pSym);
		if(((NodePair*)pNode)->m_pSecond)
			pSymSecond = ((NodePair*)pNode)->m_pSecond->eval(info, pSym);

		if(pSymFirst && pSymSecond)
		{
			if(pSymFirst->GetType() != SYMBOL_STRING)
			{
				std::cerr << linenr("Error", info) 
					<< "Only string keys are supported at the moment."
					<< std::endl;
				continue; 
			}

			const std::string& strKey = ((SymbolString*)pSymFirst)->m_strVal;

			pSymMap->m_map.insert(
				SymbolMap::t_map::value_type(strKey, pSymSecond->clone()));
		}

		safe_delete(pSymFirst, pSym, info.pGlobalSyms);
		safe_delete(pSymSecond, pSym, info.pGlobalSyms);
	}

	return pSymMap;
}

Symbol* NodeArray::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	NodeBinaryOp *pArr = (NodeBinaryOp*)m_pArr;
	std::vector<Node*> vecNodes;
	if(pArr) vecNodes = pArr->flatten(NODE_ARGS);

	SymbolArray *pSymArr = new SymbolArray;
	pSymArr->m_arr.reserve(vecNodes.size());

	for(Node* pNode : vecNodes)
	{
		if(!pNode) continue;

		Symbol *pSymbol = pNode->eval(info, pSym);
		pSymArr->m_arr.push_back(pSymbol->clone());
		pSymArr->UpdateIndices();

		safe_delete(pSymbol, pSym, info.pGlobalSyms);
	}

	return pSymArr;
}

Symbol* NodeArrayAccess::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	if(!m_pIdent /*|| m_pIdent->m_type != NODE_IDENT*/)
	{
		std::cerr << linenr("Error", info) << "Tried to access non-array." << std::endl;
		return 0;
	}

	std::string strIdent;
	Symbol* pSymbol = m_pIdent->eval(info, pSym);
	if(m_pIdent->m_type == NODE_IDENT)
		strIdent = pSymbol->m_strIdent;
	else
		strIdent = "<tmp_sym>";

	if(!pSymbol)
	{
		std::cerr << linenr("Error", info) << "Symbol \"" << strIdent
				<< "\" not found." << std::endl;
		return 0;
	}

	if(pSymbol->GetType() == SYMBOL_ARRAY)
	{
		for(Node *pIndices : m_vecIndices)
		{
			Symbol *pSymExpr = pIndices->eval(info, pSym);
			if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_INT)
			{
				std::cerr << linenr("Error", info) 
						<< "Array index has to be of integer type."
						<< std::endl;
				return 0;
			}

			int iIdx = pSymExpr->GetValInt();
			safe_delete(pSymExpr, pSym, info.pGlobalSyms);

			SymbolArray *pArr = (SymbolArray*)pSymbol;

			// convert negative indices
			if(iIdx < 0)
				iIdx = pArr->m_arr.size()  + iIdx;

			if(iIdx < 0)
			{
				std::cerr << linenr("Error", info) 
					<< "Invalid array index." 
					<< std::endl;

				return 0;
			}

			// index too high -> fill up with zeroes
			if(iIdx>=pArr->m_arr.size())
			{
	/*			std::cerr << linenr("Warning", info)
						<<  "Array index (" << iIdx
							<< ") out of bounds (array size: "
							<< pArr->m_arr.size() << ")."
							<< " Resizing."<< std::endl;
	*/

				unsigned int iOldSize = pArr->m_arr.size();
				for(unsigned int iRem=0; iRem<iIdx+1-iOldSize; ++iRem)
				{
					SymbolDouble *pNewSym = new SymbolDouble(0.);
					pNewSym->m_strName = "<const>";
					pArr->m_arr.push_back(pNewSym);
					//std::cout << "Inserting: " << iRem << std::endl;
				}
			}

			pSymbol = pArr->m_arr[iIdx];
			pArr->UpdateIndex(iIdx);
		}
	}
	else if(pSymbol->GetType() == SYMBOL_MAP)
	{
		if(m_vecIndices.size()==0)
		{
			std::cerr << linenr("Error", info) << "No key given for map." << std::endl;
			return 0;
		}
		else if(m_vecIndices.size()>1)
		{
			std::cerr << linenr("Warning", info) << "Multiple keys given for map, using first one."
						<< std::endl;
		}

		Node *pNodeKey = m_vecIndices[0];
		Symbol *pSymExpr = pNodeKey->eval(info, pSym);
		if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_STRING)
		{
			std::cerr << linenr("Error", info) << "Map key has to be of string type."
						<< std::endl;
			return 0;
		}

		const std::string& strKey = ((SymbolString*)pSymExpr)->m_strVal;

		SymbolMap *pMap = (SymbolMap*)pSymbol;
		SymbolMap::t_map::iterator iterMap = pMap->m_map.find(strKey);

		// key not yet in map -> insert it
		if(iterMap == pMap->m_map.end())
		{
			SymbolString *pNewSym = new SymbolString();
			pNewSym->m_strName = "<const>";
			iterMap = pMap->m_map.insert(SymbolMap::t_map::value_type(strKey, pNewSym)).first;
		}

		pSymbol = iterMap->second;
		pMap->UpdateIndex(strKey);
		safe_delete(pSymExpr, pSym, info.pGlobalSyms);
	}
	else
	{
		std::cerr << linenr("Error", info) << "Symbol \"" << strIdent
				<< "\" is neither an array nor a map." << std::endl;
		return 0;
	}

	return pSymbol;
}

static void uminus_inplace(Symbol* pSym, ParseInfo& info)
{
	if(!pSym) return;

	if(pSym->GetType() == SYMBOL_DOUBLE)
		((SymbolDouble*)pSym)->m_dVal = -((SymbolDouble*)pSym)->m_dVal;
	else if(pSym->GetType() == SYMBOL_INT)
		((SymbolInt*)pSym)->m_iVal = -((SymbolInt*)pSym)->m_iVal;
	else if(pSym->GetType() == SYMBOL_ARRAY)
	{
		for(Symbol* pElem : ((SymbolArray*)pSym)->m_arr)
			uminus_inplace(pElem, info);
	}
	/*else if(pSym->GetType() == SYMBOL_MAP)
	{
		for(SymbolMap::t_map::value_type& pair : ((SymbolMap*)pSym)->m_map)
			uminus_inplace(pair.second);
	}*/
	else
	{
		std::cerr << linenr("Error", info) 
			<< "Unary minus not defined for " 
			<< pSym->GetTypeName() << "."
			<< std::endl;
	}
}

Symbol* NodeUnaryOp::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	switch(m_type)
	{
		case NODE_UMINUS:
		{
			Symbol *pSymbolEval = m_pChild->eval(info, pSym);
			Symbol *pSymbol = pSymbolEval->clone();
			safe_delete(pSymbolEval, pSym, info.pGlobalSyms);

			uminus_inplace(pSymbol, info);
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

			safe_delete(pSymbolEval, pSym, info.pGlobalSyms);
			return pSymbolInt;
		}

		case NODE_STMTS:
		{
			if(m_pChild)
			{
				Symbol *pSymbol = m_pChild->eval(info, pSym);
				safe_delete(pSymbol, pSym, info.pGlobalSyms);
			}
			return 0;
		}
	}
}

Symbol* NodeBinaryOp::eval_assign(ParseInfo &info, SymbolTable *pSym) const
{
	if(m_pLeft==0 || m_pRight==0)
	{
		std::cerr << linenr("Error", info) << "NULL assignment." << std::endl;
		return 0;
	}

	Symbol *pSymbolOrg = m_pRight->eval(info, pSym);
	if(!pSymbolOrg)
	{
		std::cerr << linenr("Error", info)
				<< "Invalid rhs expression in assignment."
				<< std::endl;
		return 0;
	}

	Symbol *pSymbol = pSymbolOrg->clone();
	safe_delete(pSymbolOrg, pSym, info.pGlobalSyms);

	if(m_pLeft->m_type == NODE_IDENT)		// single variable
	{
		const std::string& strIdent = ((NodeIdent*)m_pLeft)->m_strIdent;

		Symbol* pSymGlob = info.pGlobalSyms->GetSymbol(strIdent);
		Symbol* pSymLoc = 0;
		if(!m_bGlobal)
			pSymLoc = pSym->GetSymbol(strIdent);

		if(pSymLoc && pSymGlob)
		{
			std::cerr << linenr("Warning", info) << "Symbol \"" << strIdent
					  << "\" exists in local and global scope, using local one." << std::endl;
		}

		if(pSymGlob && !pSymLoc && !m_bGlobal)
		{
			std::cerr << linenr("Warning", info) << "Overwriting global symbol \""
					<< strIdent << "\"." << std::endl;
			info.pGlobalSyms->InsertSymbol(strIdent, pSymbol);
		}
		else
		{
			if(m_bGlobal)
				info.pGlobalSyms->InsertSymbol(strIdent, pSymbol);
			else
				pSym->InsertSymbol(strIdent, pSymbol);
		}

		return pSymbol;
	}
	else								// array or map
	{
		Symbol *pSymLeft = m_pLeft->eval(info, pSym);
		if(!pSymLeft)
		{
			std::cerr << linenr("Error", info)
					<< "No array element found." << std::endl;
			return 0;
		}

		if(pSymLeft->GetType() == pSymbol->GetType())
		{
			pSymLeft->assign(pSymbol);
		}
		else if(pSymLeft->m_pArr)		// array
		{
			int iArrIdx = pSymLeft->m_iArrIdx;
			SymbolArray* pArr = pSymLeft->m_pArr;

			//std::cout << "Array: " << (void*) pArr << ", Index: " << iArrIdx << std::endl;

			if(pArr->m_arr.size() <= iArrIdx)
			{
/*						std::cerr << "Warning: Array index (" << iArrIdx
						<< ") out of bounds (array size: "
						<< pArr->m_arr.size() << ")."
						<< " Resizing."<< std::endl;
*/

				unsigned int iOldSize = pArr->m_arr.size();
				for(unsigned int iRem=0; iRem<iArrIdx+1-iOldSize; ++iRem)
				{
					SymbolDouble *pNewSym = new SymbolDouble(0.);
					pNewSym->m_strName = "<const>";
					pArr->m_arr.push_back(pNewSym);
				}
			}


			Symbol* pSymOld = pArr->m_arr[iArrIdx];
			if((void*)pSymOld != (void*)pSymLeft)
			{
				std::cerr << linenr("Error", info)
						<< "Array member mismatch." << std::endl;
				return 0;
			}


			//std::cout << pSymbol->GetType() << std::endl;
			pArr->m_arr[iArrIdx] = pSymbol;
			pSymbol->m_pArr = pArr;
			pSymbol->m_iArrIdx = iArrIdx;

			pSymOld->m_pArr = 0;
			safe_delete(pSymOld, pSym, info.pGlobalSyms);
		}
		else if(pSymLeft->m_pMap)		// map
		{
			const std::string& strMapKey = pSymLeft->m_strMapKey;
			SymbolMap* pMap = pSymLeft->m_pMap;

			SymbolMap::t_map::iterator iterMap = pMap->m_map.find(strMapKey);

			// symbol not in map -> insert a zero
			if(iterMap == pMap->m_map.end())
			{
				SymbolDouble *pNewSym = new SymbolDouble(0.);
				pNewSym->m_strName = "<const>";
				iterMap = pMap->m_map.insert(SymbolMap::t_map::value_type(strMapKey, pNewSym)).first;
			}

			Symbol* pSymOld = iterMap->second;
			if((void*)pSymOld != (void*)pSymLeft)
			{
				std::cerr << linenr("Error", info) <<
							"Map member mismatch." << std::endl;
				return 0;
			}

			pSymbol->m_pMap = pMap;
			pSymbol->m_strMapKey = strMapKey;
			iterMap->second = pSymbol;

			pSymOld->m_pMap = 0;
			safe_delete(pSymOld, pSym, info.pGlobalSyms);
		}
		else
		{
			std::cerr << linenr("Error", info)
					<< "Trying to access array/map member with no associated array/map."
					<< std::endl;
			return 0;
		}

		//safe_delete(pSymLeft, pSym, info.pGlobalSyms);
	}
	return pSymbol;
}

Symbol* NodeBinaryOp::eval_funcinit(ParseInfo &info, SymbolTable *pSym) const
{
	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	//std::cout << "Executing NODE_FUNCS for " << info.strInitScrFile << std::endl;
	NodeFunction *pToRun = 0;

	std::vector<Node*> vecFuncs0;
	vecFuncs0 = this->flatten(NODE_FUNCS);
	for(Node *_pNodeFunc : vecFuncs0)
	{
		if(!_pNodeFunc)
			continue;

		NodeFunction *pNodeFunc = (NodeFunction*)_pNodeFunc;
		pNodeFunc->m_strScrFile = info.strInitScrFile;

		if(pNodeFunc->GetName() == "__init__")
			pToRun = pNodeFunc;
		else
			vecFuncs.push_back(pNodeFunc);
	}

	// execute general entry point function
	if(pToRun)
	{
		std::string strExecBck = info.strExecFkt;

		Symbol *pSymInitRet = pToRun->eval(info, pSym);
		safe_delete(pSymInitRet, pSym, info.pGlobalSyms);

		info.strExecFkt = strExecBck;
	}

	// execute named entry point function
	if(info.strExecFkt != "")
	{
		for(NodeFunction* pFkt : vecFuncs)
		{
			if(pFkt->GetName() == info.strExecFkt)
			{
				// argument counts have to match
				if(info.pvecExecArg && pFkt->m_vecArgs.size() == info.pvecExecArg->size())
					pFkt->SetArgSyms(info.pvecExecArg);
				Symbol *pSymRet = pFkt->eval(info, pSym);
				info.pvecExecArg = 0;
				return pSymRet;
			}
		}

		std::cerr << linenr("Error", info) << "Function \"" << info.strExecFkt
				<< "\" not defined." << std::endl;
	}
	return 0;
}

Symbol* NodeBinaryOp::eval_sequential(ParseInfo &info, SymbolTable *pSym) const
{
	if(m_pLeft)
	{
		//std::cout << "left: " << m_pLeft->m_type << std::endl;
		Symbol *pSymbol = m_pLeft->eval(info, pSym);
		safe_delete(pSymbol, pSym, info.pGlobalSyms);
	}
	if(m_pRight)
	{
		//std::cout << "right: " << m_pRight->m_type << std::endl;
		Symbol *pSymbol = m_pRight->eval(info, pSym);
		safe_delete(pSymbol, pSym, info.pGlobalSyms);
	}
	return 0;
}

Symbol* NodeBinaryOp::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	switch(m_type)
	{
		case NODE_STMTS:
		case NODE_ARGS:
			return eval_sequential(info, pSym);

		case NODE_ASSIGN: return eval_assign(info, pSym);

		// should only be called once per module
		case NODE_FUNCS: return eval_funcinit(info, pSym);
	};

	Symbol *pSymbolLeft = m_pLeft->eval(info, pSym);
	Symbol *pSymbolRight = m_pRight->eval(info, pSym);
	Symbol *pSymbol = Op(pSymbolLeft, pSymbolRight, m_type);
	safe_delete(pSymbolLeft, pSym, info.pGlobalSyms);
	safe_delete(pSymbolRight, pSym, info.pGlobalSyms);

	return pSymbol;
}


Symbol* NodeFunction::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	info.pCurFunction = this;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;
	std::string strName = GetName();
	//std::cout << "in fkt " << strName << std::endl;


	SymbolTable *pLocalSym = new SymbolTable;
	if(m_pVecArgSyms)
	{
		if(m_vecArgs.size() != m_pVecArgSyms->size())
		{
			std::cerr << linenr("Error", info) << "Function \""
					<< strName << "\"" << " takes "
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
		{
			pRet = pLocalSym->GetSymbol("<ret>");
			pLocalSym->RemoveSymbolNoDelete("<ret>");
		}
	}

	//std::cout << "Local symbols for \"" << strName << "\":\n";
	//pLocalSym->print();

	delete pLocalSym;
	return pRet;
}


Symbol* NodeIf::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	Symbol *pSymExpr = 0;
	Symbol *pSymRet = 0;
	if(m_pExpr)
		pSymExpr = m_pExpr->eval(info, pSym);

	if(pSymExpr && pSymExpr->IsNotZero())
		pSymRet = (m_pIf ? m_pIf->eval(info, pSym) : 0);
	else
		pSymRet = (m_pElse ? m_pElse->eval(info, pSym) : 0);

	safe_delete(pSymExpr, pSym, info.pGlobalSyms);
	safe_delete(pSymRet, pSym, info.pGlobalSyms);

	return 0;
}


Symbol* NodeWhile::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	if(!m_pExpr) return 0;
	if(!m_pStmt) return 0;

	info.pCurLoop = this;
	while(1)
	{
		Symbol *pSymRet = 0;
		Symbol *pSymExpr = m_pExpr->eval(info, pSym);

		if(pSymExpr && pSymExpr->IsNotZero())
			pSymRet = m_pStmt->eval(info, pSym);
		else
			break;

		safe_delete(pSymRet, pSym, info.pGlobalSyms);
		safe_delete(pSymExpr, pSym, info.pGlobalSyms);

		if(info.bWantBreak)
		{
			info.bWantBreak = 0;
			break;
		}
		if(info.bWantContinue)
		{
			info.bWantContinue = 0;
			continue;
		}
	}
	info.pCurLoop = 0;

	return 0;
}

Symbol* NodeRangedFor::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;
	if(!m_pIdent || !m_pExpr || !m_pStmt) return 0;

	if(m_pIdent->m_type != NODE_IDENT)
	{
		std::cerr << linenr("Error", info) << "Range-based for loop needs identifier."
					<< std::endl;
		return 0;
	}

	Symbol *pSymRet = 0;
	Symbol *_pArr = m_pExpr->eval(info, pSym);
	if(!_pArr)
	{
		std::cerr << linenr("Error", info)
				<< "Invalid array for loop." << std::endl;
		return 0;
	}

	if(_pArr->GetType() != SYMBOL_ARRAY)
	{
		std::cerr << linenr("Error", info)
				<< "Range-based for loop needs array." << std::endl;
		safe_delete(_pArr, pSym, info.pGlobalSyms);
		return 0;
	}

	SymbolArray *pArr = (SymbolArray*)_pArr;


	const std::string& strIdent = ((NodeIdent*)m_pIdent)->m_strIdent;

	SymbolInt *pSymIter = new SymbolInt(0);
	std::string strIter = "<cur_iter_" + strIdent + ">";
	pSym->InsertSymbol(strIter, pSymIter);

	info.pCurLoop = this;
	for(unsigned int iArr=0; iArr<pArr->m_arr.size(); ++iArr)
	{
		Symbol *pSymInArr = pArr->m_arr[iArr];
		pSym->InsertSymbol(strIdent, pSymInArr);

		Symbol *pBodyRet = m_pStmt->eval(info, pSym);
		safe_delete(pBodyRet, pSym, info.pGlobalSyms);


		// write back symbol in case an assignment has taken place
		Symbol *pNewSym = pSym->GetSymbol(strIdent);
		if(pSymInArr != pNewSym)
		{
			pArr->m_arr[iArr] = pNewSym;
			//delete pSymInArr;
			pSymInArr = pNewSym;
		}
		pSym->RemoveSymbolNoDelete(strIdent);


		++pSymIter->m_iVal;

		if(info.bWantBreak)
		{
			info.bWantBreak = 0;
			break;
		}
		if(info.bWantContinue)
		{
			info.bWantContinue = 0;
			continue;
		}
	}
	info.pCurLoop = 0;

	pSym->RemoveSymbol(strIter);

	//std::cout << "ranged for:" << pArr->m_strName << ", " << pArr->m_strIdent << std::endl;
	safe_delete(_pArr, pSym, info.pGlobalSyms);
	return 0;
}
