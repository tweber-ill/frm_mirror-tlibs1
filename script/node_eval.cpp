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
	if(pSym->GetName() == T_STR"<const>")
		return;

	// don't delete array or map members
	if(pSym->GetArrPtr() || pSym->GetMapPtr())
		return;

	// don't delete symbols in table
	bool bIsInTable = pSymTab->IsPtrInMap(pSym);
	bool bIsInGlobTable = pSymTabGlob->IsPtrInMap(pSym);
	if(!bIsInTable && !bIsInGlobTable)
	{
		//G_COUT << "deleting " << (void*)pSym << ": " << pSym->GetType() << std::endl;
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
	pSym->InsertSymbol(T_STR"<ret>", pRet ? pRet : 0);

	info.bWantReturn = 1;
	return pRet;
}

Symbol* NodeBreak::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	if(info.pCurLoop==0)
	{
		G_CERR << linenr(T_STR"Error", info)
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
		G_CERR << linenr(T_STR"Error", info)
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
		G_CERR << linenr(T_STR"Warning", info) << "Symbol \"" << m_strIdent
				<< "\" exists in local and global scope, using local one." << std::endl;
	}

	// if no local symbol is available, use global symbol instead
	if(!pSymbol)
		pSymbol = pSymbolGlob;

	if(!pSymbol)
	{
		G_CERR << linenr(T_STR"Error", info) << "Symbol \"" << m_strIdent
				<< "\" not in symbol table." << std::endl;
		return 0;
	}

	//if(pSymbol->GetIdent().size()==0)	// !! cur_iter needs this symbol's identifier !!
		pSymbol->SetIdent(m_strIdent);
	return pSymbol;
}

Symbol* NodeCall::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	info.pCurCaller = this;

	if(m_pIdent->GetType() != NODE_IDENT)
		return 0;
	//if(m_pArgs->GetType() != NODE_ARGS)
	//	return 0;

	NodeIdent* pIdent = (NodeIdent*) m_pIdent;
	const t_string& strFkt = pIdent->GetIdent();
	//G_COUT << "call to " << strFkt << " with " << m_vecArgs.size() << " arguments." << std::endl;


	bool bCallUserFkt = 0;
	// user-defined function
	NodeFunction *pFkt = info.GetFunction(strFkt);;
	if(pFkt)
		bCallUserFkt = 1;

	/*if(!bCallUserFkt)
	{
		G_CERR << "Error: Trying to call unknown function \" << strFkt << \"."
					<< std::endl;
		return 0;
	}*/


	SymbolArray arrArgs;
	arrArgs.SetDontDel(1);
	std::vector<Symbol*> &vecArgSyms = arrArgs.GetArr();
	for(Node* pNode : m_vecArgs)
	{
		// TODO: Unpack operation for vector.
		if(pNode->GetType() == NODE_UNPACK)
		{
			Node *pChild = ((NodeUnaryOp*)pNode)->GetChild();
			if(!pChild)
			{
				G_CERR << linenr(T_STR"Error", info)
						<< "Invalid symbol to unpack." << std::endl;
				continue;
			}
			Symbol *pSymChild = pChild->eval(info, pSym);
			if(pSymChild->GetType() != SYMBOL_MAP)
			{
				G_CERR << linenr(T_STR"Error", info)
						<< "Cannot unpack non-map." << std::endl;
				continue;
			}

			SymbolMap::t_map& mapSyms = ((SymbolMap*)pSymChild)->GetMap();

			std::vector<t_string> vecParamNames = pFkt->GetParamNames();
			for(unsigned int iParam=vecArgSyms.size(); iParam<vecParamNames.size(); ++iParam)
			{
				const t_string& strParamName = vecParamNames[iParam];
				SymbolMap::t_map::iterator iter = mapSyms.find(strParamName);
				if(iter == mapSyms.end())
				{
					G_CERR << linenr(T_STR"Error", info)
							<< "Parameter \"" << strParamName
							<< "\" not in map." << std::endl;

					vecArgSyms.push_back(new SymbolDouble(0.));
					continue;
				}

				vecArgSyms.push_back(iter->second->clone());
			}

			safe_delete(pSymChild, pSym, info.pGlobalSyms);
		}
		else
		{
			Symbol *pSymbol = pNode->eval(info, pSym);
			//G_COUT << "argument: " << pSymbol->print() << std::endl;

			vecArgSyms.push_back(pSymbol);
		}
	}

	Symbol* pFktRet = 0;
	if(bCallUserFkt)	// call user-defined function
	{
		//pFkt->SetArgSyms(&vecArgSyms);
		pSym->InsertSymbol(T_STR"<args>", &arrArgs);
		pFktRet = pFkt->eval(info, pSym);
		pSym->RemoveSymbolNoDelete(T_STR"<args>");
		info.bWantReturn = 0;
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

	G_CERR << linenr(T_STR"Error", info)
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
		if(pNode->GetType() != NODE_PAIR)
		{
			G_CERR << linenr(T_STR"Error", info)
				<< "Maps have to consist of key-value pairs." 
				<< std::endl;
			continue;
		}

		Symbol* pSymFirst = 0;
		Symbol* pSymSecond = 0;

		if(((NodePair*)pNode)->GetFirst())
			pSymFirst = ((NodePair*)pNode)->GetFirst()->eval(info, pSym);
		if(((NodePair*)pNode)->GetSecond())
			pSymSecond = ((NodePair*)pNode)->GetSecond()->eval(info, pSym);

		if(pSymFirst && pSymSecond)
		{
			if(pSymFirst->GetType() != SYMBOL_STRING)
			{
				G_CERR << linenr(T_STR"Error", info)
					<< "Only string keys are supported at the moment."
					<< std::endl;
				continue; 
			}

			const t_string& strKey = ((SymbolString*)pSymFirst)->GetVal();

			pSymMap->GetMap().insert(
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
	pSymArr->GetArr().reserve(vecNodes.size());

	for(Node* pNode : vecNodes)
	{
		if(!pNode) continue;

		Symbol *pSymbol = pNode->eval(info, pSym);
		pSymArr->GetArr().push_back(pSymbol->clone());
		pSymArr->UpdateIndices();

		safe_delete(pSymbol, pSym, info.pGlobalSyms);
	}

	return pSymArr;
}

Symbol* NodeRange::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	G_CERR << linenr(T_STR"Error", info)
				<< "Range should not be evaluated directly."
				<< std::endl;
	return 0;
}

void NodeRange::GetRangeIndices(ParseInfo &info, SymbolTable *pSym,
					t_int iMaxLen, t_int& iBeginIdx, t_int& iEndIdx)
{
	/*if(!pArr)
	{
		G_CERR << linenr(T_STR"Error", info)
					<< "Range operation needs an array."
					<< std::endl;
		return;
	}*/

	if(m_rangetype == RANGE_FULL)
	{
		iBeginIdx = 0;
		iEndIdx = iMaxLen;
	}
	else if(m_rangetype == RANGE_BEGINEND)
	{
		if(m_pBegin==0 || m_pEnd==0)
		{
			G_CERR << linenr(T_STR"Error", info)
						<< "Invalid range."
						<< std::endl;
			return;
		}

		Symbol *pSymBeg = m_pBegin->eval(info, pSym);
		Symbol *pSymEnd = m_pEnd->eval(info, pSym);

		iBeginIdx = pSymBeg->GetValInt();
		iEndIdx = pSymEnd->GetValInt();

		// convert negative indices
		if(iBeginIdx < 0) iBeginIdx = iMaxLen + iBeginIdx;
		if(iEndIdx < 0) iEndIdx = iMaxLen + iEndIdx;


		if(iBeginIdx<0 || iBeginIdx>=iMaxLen)
		{
			G_CERR << linenr(T_STR"Error", info)
						<< "Lower array index out of bounds: "
						<< iBeginIdx << ". Adjusting to lower limit."
						<< std::endl;
			iBeginIdx = 0;
		}
		if(iEndIdx<-1 || iEndIdx>iMaxLen)
		{
			G_CERR << linenr(T_STR"Error", info)
						<< "Upper array index out of bounds: "
						<< iEndIdx << ". Adjusting to upper limit."
						<< std::endl;
			iEndIdx = iMaxLen;
		}

		//G_COUT << "begin: " << iBeginIdx << ", end: " << iEndIdx << std::endl;

		safe_delete(pSymBeg, pSym, info.pGlobalSyms);
		safe_delete(pSymEnd, pSym, info.pGlobalSyms);
	}
	else
	{
		G_CERR << linenr(T_STR"Error", info)
					<< "Invalid range operation."
					<< std::endl;
	}
}

unsigned int get_max_cols(SymbolArray* pArr, std::vector<unsigned int>* pvecCols=0)
{
	unsigned int iCols = 0;
	for(Symbol* pSym : pArr->GetArr())
	{
		if(!pSym)
		{
			if(pvecCols) pvecCols->push_back(0);
			continue;
		}
		if(pSym->GetType() == SYMBOL_ARRAY)
			iCols = std::max<unsigned int>(iCols, ((SymbolArray*)pSym)->GetArr().size());
		else
			iCols = std::max<unsigned int>(iCols, 1);

		if(pvecCols) pvecCols->push_back(iCols);
	}

	return iCols;
}

Symbol* get_mat_elem(SymbolArray* pArr, unsigned int iLine, unsigned int iCol)
{
	if(iLine >= pArr->GetArr().size())
		return 0;

	Symbol* pLine = pArr->GetArr()[iLine];
	if(pLine->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrLine = (SymbolArray*)pLine;
		if(iCol >= pArrLine->GetArr().size())
			return 0;

		return pArrLine->GetArr()[iCol]->clone();
	}
	else	// single scalar element
	{
		if(iCol == 0)
			return pLine;
		else
			return 0;
	}
}

SymbolArray* transpose(SymbolArray* pArr, std::vector<unsigned int>* pvecCols=0)
{
	unsigned int iLines = pArr->GetArr().size();
	unsigned int iCols = get_max_cols(pArr, pvecCols);

	SymbolArray *pArrNew = new SymbolArray();
	pArrNew->GetArr().reserve(iCols);

	for(unsigned int iCol=0; iCol<iCols; ++iCol)
	{
		SymbolArray* pArrCol = new SymbolArray();
		pArrCol->GetArr().reserve(iLines);

		for(unsigned int iLine=0; iLine<iLines; ++iLine)
		{
			Symbol *pElem = get_mat_elem(pArr, iLine, iCol);
			if(!pElem) pElem = new SymbolInt(0);
			pArrCol->GetArr().push_back(pElem);
		}

		pArrCol->UpdateIndices();
		pArrNew->GetArr().push_back(pArrCol);
	}

	pArrNew->UpdateIndices();
	return pArrNew;
}

Symbol* NodeArrayAccess::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	if(!m_pIdent /*|| m_pIdent->GetType() != NODE_IDENT*/)
	{
		G_CERR << linenr(T_STR"Error", info) << "Tried to access non-array." << std::endl;
		return 0;
	}

	t_string strIdent;
	Symbol* pSymbol = m_pIdent->eval(info, pSym);
	if(m_pIdent->GetType() == NODE_IDENT)
		strIdent = pSymbol->GetIdent();
	else
		strIdent = T_STR"<tmp_sym>";

	if(!pSymbol)
	{
		G_CERR << linenr(T_STR"Error", info) << "Symbol \"" << strIdent
				<< "\" not found." << std::endl;
		return 0;
	}

	if(pSymbol->GetType() == SYMBOL_ARRAY)
	{
		bool bCreatedSym = 0;
		bool bEvenIndex = 0;
		bool bAlreadyTransposed = 0;

		for(Node *pIndices : m_vecIndices)
		{
			SymbolArray *pArr = (SymbolArray*)pSymbol;

			if(bCreatedSym /*&& bEvenIndex*/)
			{
				pArr = transpose(pArr);
				delete pSymbol;
				pSymbol = pArr;

				bAlreadyTransposed = 1;
			}

			// TODO: assigment for ranged access
			if(pIndices->GetType() == NODE_RANGE)	// range index
			{
				NodeRange *pRange = (NodeRange*)pIndices;
				t_int iBeginIdx = 0, iEndIdx = 0;
				pRange->GetRangeIndices(info, pSym, pArr->GetArr().size(), iBeginIdx, iEndIdx);

				t_int iStep = 1;
				t_int iSize = iEndIdx - iBeginIdx;
				if(iEndIdx < iBeginIdx)
				{
					iSize = -iSize;
					iStep = -1;
				}

				SymbolArray *pSubArr = new SymbolArray();
				pSubArr->GetArr().reserve(iSize);

				for(t_int iIdx=iBeginIdx, iNewIdx=0; iIdx!=iEndIdx && iNewIdx<iSize; iIdx+=iStep, ++iNewIdx)
				{
					Symbol *pElemClone = pArr->GetArr()[iIdx]->clone();
					pSubArr->GetArr().push_back(pElemClone);
					pSubArr->UpdateIndex(iNewIdx);
				}

				if(bAlreadyTransposed)
				{
					SymbolArray *pOrgSubArr = pSubArr;
					pSubArr = transpose(pSubArr);
					delete pOrgSubArr;

					bAlreadyTransposed = 0;
				}

				if(bCreatedSym)
				{
					delete pSymbol;
					pSymbol = 0;
					bCreatedSym = 0;
				}

				safe_delete(pSymbol, pSym, info.pGlobalSyms);
				pSymbol = pSubArr;
				bCreatedSym = 1;
			}
			else								// integer index
			{
				Symbol *pSymExpr = pIndices->eval(info, pSym);
				if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_INT)
				{
					G_CERR << linenr(T_STR"Error", info)
							<< "Array index has to be of integer type."
							<< std::endl;
					return 0;
				}

				t_int iIdx = pSymExpr->GetValInt();
				safe_delete(pSymExpr, pSym, info.pGlobalSyms);

				// convert negative indices
				if(iIdx < 0)
					iIdx = pArr->GetArr().size()  + iIdx;

				if(iIdx < 0)
				{
					G_CERR << linenr(T_STR"Error", info)
						<< "Invalid array index."
						<< std::endl;

					return 0;
				}

				// index too high -> fill up with zeroes
				if(iIdx>=pArr->GetArr().size())
				{
					unsigned int iOldSize = pArr->GetArr().size();
					for(unsigned int iRem=0; iRem<iIdx+1-iOldSize; ++iRem)
					{
						SymbolDouble *pNewSym = new SymbolDouble(0.);
						pNewSym->SetName(T_STR"<const>");
						pArr->GetArr().push_back(pNewSym);
						//G_COUT << "Inserting: " << iRem << std::endl;
					}
				}

				if((void*)pSymbol != (void*)pArr)
					safe_delete(pSymbol, pSym, info.pGlobalSyms);

				pSymbol = pArr->GetArr()[iIdx];
				pArr->UpdateIndex(iIdx);
			}

			bEvenIndex = !bEvenIndex;
		}
	}
	else if(pSymbol->GetType() == SYMBOL_MAP)
	{
		if(m_vecIndices.size()==0)
		{
			G_CERR << linenr(T_STR"Error", info) << "No key given for map." << std::endl;
			return 0;
		}
		else if(m_vecIndices.size()>1)
		{
			G_CERR << linenr(T_STR"Warning", info) << "Multiple keys given for map, using first one."
						<< std::endl;
		}

		Node *pNodeKey = m_vecIndices[0];
		Symbol *pSymExpr = pNodeKey->eval(info, pSym);
		if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_STRING)
		{
			G_CERR << linenr(T_STR"Error", info) << "Map key has to be of string type."
						<< std::endl;
			return 0;
		}

		const t_string& strKey = ((SymbolString*)pSymExpr)->GetVal();

		SymbolMap *pMap = (SymbolMap*)pSymbol;
		SymbolMap::t_map::iterator iterMap = pMap->GetMap().find(strKey);

		// key not yet in map -> insert it
		if(iterMap == pMap->GetMap().end())
		{
			SymbolString *pNewSym = new SymbolString();
			pNewSym->SetName(T_STR"<const>");
			iterMap = pMap->GetMap().insert(SymbolMap::t_map::value_type(strKey, pNewSym)).first;
		}

		if((void*)pSymbol != (void*)iterMap->second)
			safe_delete(pSymbol, pSym, info.pGlobalSyms);

		pSymbol = iterMap->second;
		pMap->UpdateIndex(strKey);
		safe_delete(pSymExpr, pSym, info.pGlobalSyms);
	}
	else if(pSymbol->GetType() == SYMBOL_STRING)
	{
		SymbolString *pStr = (SymbolString*)pSymbol;
		t_int iStrLen = pStr->GetVal().length();

		if(m_vecIndices.size()!=1)
		{
			G_CERR << linenr(T_STR"Error", info)
					<< "Need exactly one string index."
					<< std::endl;
			return 0;
		}

		Node *pIndices = m_vecIndices[0];

		if(pIndices->GetType() == NODE_RANGE)	// range index
		{
			NodeRange *pRange = (NodeRange*)pIndices;
			t_int iBeginIdx = 0, iEndIdx = 0;
			pRange->GetRangeIndices(info, pSym, iStrLen, iBeginIdx, iEndIdx);

			t_int iStep = 1;
			t_int iSize = iEndIdx - iBeginIdx;
			if(iEndIdx < iBeginIdx)
			{
				iSize = -iSize;
				iStep = -1;
			}

			t_string strNew;
			strNew.resize(iSize);

			for(t_int iIdx=iBeginIdx, iNewIdx=0; iIdx!=iEndIdx && iNewIdx<iSize; iIdx+=iStep, ++iNewIdx)
				strNew[iNewIdx] = pStr->GetVal()[iIdx];

			SymbolString *pSubStr = new SymbolString(strNew);
			safe_delete(pSymbol, pSym, info.pGlobalSyms);
			pSymbol = pSubStr;
		}
		else								// integer index
		{
			Symbol *pSymExpr = pIndices->eval(info, pSym);
			if(pSymExpr==0 || pSymExpr->GetType()!=SYMBOL_INT)
			{
				G_CERR << linenr(T_STR"Error", info)
						<< "String index has to be of integer type."
						<< std::endl;
				return 0;
			}

			t_int iIdx = pSymExpr->GetValInt();
			safe_delete(pSymExpr, pSym, info.pGlobalSyms);

			// convert negative indices
			if(iIdx < 0)
				iIdx = iStrLen  + iIdx;

			if(iIdx < 0 || iIdx >= iStrLen)
			{
				G_CERR << linenr(T_STR"Error", info)
					<< "String index out of bounds."
					<< std::endl;

				return 0;
			}

			t_string strNew;
			strNew += pStr->GetVal()[iIdx];
			safe_delete(pSymbol, pSym, info.pGlobalSyms);
			pSymbol = new SymbolString(strNew);
		}
	}
	else
	{
		G_CERR << linenr(T_STR"Error", info) << "Symbol \"" << strIdent
				<< "\" is neither an array nor a map." << std::endl;
		return 0;
	}

	return pSymbol;
}

static void uminus_inplace(Symbol* pSym, ParseInfo& info)
{
	if(!pSym) return;

	if(pSym->GetType() == SYMBOL_DOUBLE)
		((SymbolDouble*)pSym)->SetVal(-((SymbolDouble*)pSym)->GetVal());
	else if(pSym->GetType() == SYMBOL_INT)
		((SymbolInt*)pSym)->SetVal(-((SymbolInt*)pSym)->GetVal());
	else if(pSym->GetType() == SYMBOL_ARRAY)
	{
		for(Symbol* pElem : ((SymbolArray*)pSym)->GetArr())
			uminus_inplace(pElem, info);
	}
	/*else if(pSym->GetType() == SYMBOL_MAP)
	{
		for(SymbolMap::t_map::value_type& pair : ((SymbolMap*)pSym)->GetMap())
			uminus_inplace(pair.second);
	}*/
	else
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "Unary minus not defined for " 
			<< pSym->GetTypeName() << "."
			<< std::endl;
	}
}

Symbol* NodeUnaryOp::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	switch(GetType())
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
				pSymbolInt->SetVal(!((SymbolDouble*)pSymbolEval)->GetVal());
			else if(pSymbolEval->GetType() == SYMBOL_INT)
				pSymbolInt->SetVal(!((SymbolInt*)pSymbolEval)->GetVal());

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

		case NODE_UNPACK:
		{
			G_CERR << linenr(T_STR"Warning", info)
					<< "Unpack operation only allowed in function call. Ignoring."
					<< std::endl;

			if(m_pChild)
				return m_pChild->eval(info, pSym);
		}
		default:
			G_CERR << linenr(T_STR"Warning", info)
					<< "Unknown node type: " << GetType()
					<< std::endl;
			break;
	}

	return 0;
}

Symbol* NodeBinaryOp::eval_assign(ParseInfo &info, SymbolTable *pSym,
 				Node *pLeft, Node *pRight, Symbol *pSymRightAlt,
				const bool *pbGlob) const
{
	if(pLeft==0) pLeft = m_pLeft;
	if(pRight==0) pRight = m_pRight;
	if(pbGlob==0) pbGlob = &m_bGlobal;

	if(pLeft==0 || pRight==0)
	{
		G_CERR << linenr(T_STR"Error", info) << "NULL assignment." << std::endl;
		return 0;
	}

//	std::cout << pLeft->GetType() << std::endl;

	Symbol *pSymbolOrg = 0;
	if(pSymRightAlt)	// use RHS symbol if given instead of RHS node
		pSymbolOrg = pSymRightAlt;
	else
		pSymbolOrg = pRight->eval(info, pSym);

	if(!pSymbolOrg)
	{
		G_CERR << linenr(T_STR"Error", info)
				<< "Invalid rhs expression in assignment."
				<< std::endl;
		return 0;
	}

	Symbol *pSymbol = pSymbol = pSymbolOrg->clone();
	safe_delete(pSymbolOrg, pSym, info.pGlobalSyms);

	if(pLeft->GetType() == NODE_IDENT)		// single variable
	{
		const t_string& strIdent = ((NodeIdent*)pLeft)->GetIdent();

		Symbol* pSymGlob = info.pGlobalSyms->GetSymbol(strIdent);
		Symbol* pSymLoc = 0;
		if(!*pbGlob)
			pSymLoc = pSym->GetSymbol(strIdent);

		if(pSymLoc && pSymGlob)
		{
			G_CERR << linenr(T_STR"Warning", info) << "Symbol \"" << strIdent
					  << "\" exists in local and global scope, using local one." << std::endl;
		}

		if(pSymGlob && !pSymLoc && !*pbGlob)
		{
			G_CERR << linenr(T_STR"Warning", info) << "Overwriting global symbol \""
					<< strIdent << "\"." << std::endl;
			info.pGlobalSyms->InsertSymbol(strIdent, pSymbol);
		}
		else
		{
			if(*pbGlob)
				info.pGlobalSyms->InsertSymbol(strIdent, pSymbol);
			else
				pSym->InsertSymbol(strIdent, pSymbol);
		}

		return pSymbol;
	}
	else if(pLeft->GetType() == NODE_ARRAY)				// e.g. [a,b] = [1,2];
	{	// TODO: check that LHS does not want eval!
		if(pSymbol->GetType() != SYMBOL_ARRAY)
		{
			G_CERR << linenr(T_STR"Error", info)
				<< "Assignment needs an array on the right-hand side."
				<< std::endl;
			return pSymbol;
		}
		SymbolArray *pArrRight = (SymbolArray*)pSymbol;

		std::vector<Node*> vecLeftArgs = ((NodeBinaryOp*)((NodeArray*)pLeft)->GetArr())->flatten(NODE_ARGS);
		unsigned int iArrSize = vecLeftArgs.size();
		if(vecLeftArgs.size() != pArrRight->GetArr().size())
		{
                        G_CERR << linenr(T_STR"Warning", info)
				<< "Size mismatch between assigned and returned array: "
				<< vecLeftArgs.size() << " != " << pArrRight->GetArr().size() << "."
				<< std::endl;

			iArrSize = std::min(vecLeftArgs.size(), pArrRight->GetArr().size());
		}

		for(unsigned int iArr=0; iArr<vecLeftArgs.size(); ++iArr)
		{
			Node *pNodeLeft = vecLeftArgs[iArr];
			Symbol *pSymRight = pArrRight->GetArr()[iArr];

			//std::cout << ((NodeIdent*)pNode)->GetIdent() << std::endl;
			eval_assign(info, pSym, pNodeLeft, 0, pSymRight, pbGlob);
		}

		return pSymbol;
	}
	else								// array or map
	{
		//std::cout << "in assign -> else" << std::endl;
		Symbol *pSymLeft = pLeft->eval(info, pSym);
		if(!pSymLeft)
		{
			G_CERR << linenr(T_STR"Error", info)
					<< "No array element found." << std::endl;
			return 0;
		}

		if(pSymLeft->GetType() == pSymbol->GetType())
		{
			pSymLeft->assign(pSymbol);
		}
		else if(pSymLeft->GetArrPtr())		// array
		{
			t_int iArrIdx = pSymLeft->GetArrIdx();
			SymbolArray* pArr = pSymLeft->GetArrPtr();

			//G_COUT << "Array: " << (void*) pArr << ", Index: " << iArrIdx << std::endl;

			if(pArr->GetArr().size() <= iArrIdx)
			{
/*						G_CERR << "Warning: Array index (" << iArrIdx
						<< ") out of bounds (array size: "
						<< pArr->GetArr().size() << ")."
						<< " Resizing."<< std::endl;
*/

				unsigned int iOldSize = pArr->GetArr().size();
				for(unsigned int iRem=0; iRem<iArrIdx+1-iOldSize; ++iRem)
				{
					SymbolDouble *pNewSym = new SymbolDouble(0.);
					pNewSym->SetName(T_STR"<const>");
					pArr->GetArr().push_back(pNewSym);
				}
			}


			Symbol* pSymOld = pArr->GetArr()[iArrIdx];
			if((void*)pSymOld != (void*)pSymLeft)
			{
				G_CERR << linenr(T_STR"Error", info)
						<< "Array member mismatch." << std::endl;
				return 0;
			}


			//G_COUT << pSymbol->GetType() << std::endl;
			pArr->GetArr()[iArrIdx] = pSymbol;
			pSymbol->SetArrPtr(pArr);
			pSymbol->SetArrIdx(iArrIdx);

			pSymOld->SetArrPtr(0);
			safe_delete(pSymOld, pSym, info.pGlobalSyms);
		}
		else if(pSymLeft->GetMapPtr())		// map
		{
			const t_string& strMapKey = pSymLeft->GetMapKey();
			SymbolMap* pMap = pSymLeft->GetMapPtr();

			SymbolMap::t_map::iterator iterMap = pMap->GetMap().find(strMapKey);

			// symbol not in map -> insert a zero
			if(iterMap == pMap->GetMap().end())
			{
				SymbolDouble *pNewSym = new SymbolDouble(0.);
				pNewSym->SetName(T_STR"<const>");
				iterMap = pMap->GetMap().insert(SymbolMap::t_map::value_type(strMapKey, pNewSym)).first;
			}

			Symbol* pSymOld = iterMap->second;
			if((void*)pSymOld != (void*)pSymLeft)
			{
				G_CERR << linenr(T_STR"Error", info) <<
							"Map member mismatch." << std::endl;
				return 0;
			}

			pSymbol->SetMapPtr(pMap);
			pSymbol->SetMapKey(strMapKey);
			iterMap->second = pSymbol;

			pSymOld->SetMapPtr(0);
			safe_delete(pSymOld, pSym, info.pGlobalSyms);
		}
		else
		{
			G_CERR << linenr(T_STR"Error", info)
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

	//G_COUT << "Executing NODE_FUNCS for " << info.strInitScrFile << std::endl;
	NodeFunction *pToRun = 0;

	std::vector<Node*> vecFuncs0;
	vecFuncs0 = this->flatten(NODE_FUNCS);
	for(Node *_pNodeFunc : vecFuncs0)
	{
		if(!_pNodeFunc)
			continue;

		NodeFunction *pNodeFunc = (NodeFunction*)_pNodeFunc;
		pNodeFunc->SetScrFile(info.strInitScrFile);

		const t_string& strFktName = pNodeFunc->GetName();

		if(strFktName == T_STR"__init__" || strFktName == T_STR"module_init")
		{
			pToRun = pNodeFunc;
		}
		else if(info.GetFunction(strFktName))
		{
			if(strFktName != "main")
				G_CERR << linenr(T_STR"Warning", info)
					<< "Function \"" << strFktName
					<< "\" redefined in \"" << info.strInitScrFile << "\"."
					<< " Ignoring." << std::endl;
		}
		else
		{
			if(has_ext_call(strFktName))
			{
				G_CERR << linenr(T_STR"Warning", info)
					<< "Function \"" << strFktName
					<< "\" in \"" << info.strInitScrFile
					<< "\" overwrites a system function."
					<< std::endl;
			}
			vecFuncs.push_back(pNodeFunc);
		}
	}

	// execute general entry point function
	if(pToRun)
	{
		t_string strExecBck = info.strExecFkt;

		Symbol *pSymInitRet = pToRun->eval(info, /*pSym*/0);
		safe_delete(pSymInitRet, pSym, info.pGlobalSyms);

		info.strExecFkt = strExecBck;
	}

	// execute named entry point function
	if(info.strExecFkt != T_STR"")
	{
		for(NodeFunction* pFkt : vecFuncs)
		{
			if(pFkt->GetName() == info.strExecFkt)
			{
				if(pFkt->GetArgVec().size() == 0)
					pSym->RemoveSymbolNoDelete(T_STR"<args>");

				Symbol *pSymRet = pFkt->eval(info, pSym);
				if(pSym) pSym->RemoveSymbolNoDelete(T_STR"<args>");

				return pSymRet;
			}
		}

		G_CERR << linenr(T_STR"Error", info) << "Function \"" << info.strExecFkt
				<< "\" not defined." << std::endl;
	}
	return 0;
}

Symbol* NodeBinaryOp::eval_recursive(ParseInfo &info, SymbolTable *pSym) const
{
	if(m_pLeft)
	{
		//G_COUT << "left: " << m_pLeft->GetType() << std::endl;
		Symbol *pSymbol = m_pLeft->eval(info, pSym);
		safe_delete(pSymbol, pSym, info.pGlobalSyms);
	}
	if(m_pRight)
	{
		//G_COUT << "right: " << m_pRight->GetType() << std::endl;
		Symbol *pSymbol = m_pRight->eval(info, pSym);
		safe_delete(pSymbol, pSym, info.pGlobalSyms);
	}
	return 0;
}

Symbol* NodeBinaryOp::eval_sequential(ParseInfo &info, SymbolTable *pSym) const
{
	for(Node* pNode : m_vecNodesFlat)
	{
		if(info.IsExecDisabled()) break;

		Symbol *pSymbol = pNode->eval(info, pSym);
		safe_delete(pSymbol, pSym, info.pGlobalSyms);
	}
	return 0;
}

Symbol* NodeBinaryOp::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	if(m_vecNodesFlat.size() != 0)
	{
		return eval_sequential(info, pSym);
	}
	else
	{
		switch(GetType())
		{
			case NODE_STMTS:
			//case NODE_ARGS:
			{
				G_CERR << linenr(T_STR"Info", info)
						<< "Should better be evaluated sequentially in interpreter."
						<< std::endl;
				return eval_recursive(info, pSym);
			}

			case NODE_ASSIGN: return eval_assign(info, pSym);

			// should only be called once per module
			case NODE_FUNCS: return eval_funcinit(info, pSym);

			default:
				break;
		};

		Symbol *pSymbolLeft = m_pLeft->eval(info, pSym);
		Symbol *pSymbolRight = 0;
		Symbol *pSymbol = 0;

		// optimisation: 0 && x == 0
		if(GetType() == NODE_LOG_AND && pSymbolLeft->GetValInt()==0)
			pSymbol = new SymbolInt(0);
		// optimisation: 1 || x == 1
		else if(GetType() == NODE_LOG_OR && pSymbolLeft->GetValInt()==1)
			pSymbol = new SymbolInt(1);
		else
		{
			pSymbolRight = m_pRight->eval(info, pSym);
			pSymbol = Op(pSymbolLeft, pSymbolRight, GetType());
		}

		safe_delete(pSymbolLeft, pSym, info.pGlobalSyms);
		safe_delete(pSymbolRight, pSym, info.pGlobalSyms);

		return pSymbol;
	}
}


Symbol* NodeFunction::eval(ParseInfo &info, SymbolTable* pTableSup) const
{
	if(info.IsExecDisabled()) return 0;
	info.pCurFunction = this;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;
	const t_string& strName = GetName();
	//G_COUT << "in fkt " << strName << std::endl;

	SymbolTable *pLocalSym = new SymbolTable;

	SymbolArray* pArgs = 0;
	if(pTableSup)
		pArgs = (SymbolArray*)pTableSup->GetSymbol(T_STR"<args>");
	if(pArgs)
	{
		const std::vector<Symbol*> *pVecArgSyms = &pArgs->GetArr();

		/*if(m_vecArgs.size() != pVecArgSyms->size())
		{
			G_CERR << linenr(T_STR"Error", info) << "Function \""
					<< strName << "\"" << " takes "
					<< m_vecArgs.size() << " arguments, but "
					<< pVecArgSyms->size() << " given."
					<< std::endl;
		}*/

		unsigned int iArgSize = /*std::min(*/m_vecArgs.size()/*, pVecArgSyms->size())*/;
		for(unsigned int iArg=0; iArg<iArgSize; ++iArg)
		{
			NodeIdent* pIdent = (NodeIdent*)m_vecArgs[iArg];
			const Node* pDefArg = pIdent->GetDefArg();

			Symbol *pSymbol = 0;

			if(iArg < pVecArgSyms->size())		// argument given by caller
				pSymbol = (*pVecArgSyms)[iArg];
			if(pSymbol==0 && pDefArg)		// default argument
			{
				pSymbol = pDefArg->eval(info, pTableSup);

				if(iArg < pVecArgSyms->size() && pSymbol)
					G_CERR << linenr(T_STR"Warning", info)
						<< "Given argument \"" << pIdent->GetIdent() 
						<< "\" for function \"" 
						<< strName << "\" not valid. " 
						<< "Using default argument." << std::endl;
			}
			if(pSymbol==0)
			{
				G_CERR << linenr(T_STR"Error", info) << "Argument \""
					<< pIdent->GetIdent() << "\" for function \""
					<< strName << "\" not given." << std::endl;
				continue;
			}

			//G_COUT << "arg: " << pIdent->GetIdent() << std::endl;
			pLocalSym->InsertSymbol(pIdent->GetIdent(), pSymbol?pSymbol->clone():0);
		}
	}


	Symbol *pRet = 0;
	if(m_pStmts)
	{
		pRet = m_pStmts->eval(info, pLocalSym);
		if(!pRet)
		{
			pRet = pLocalSym->GetSymbol(T_STR"<ret>");
			pLocalSym->RemoveSymbolNoDelete(T_STR"<ret>");
		}
	}

	//G_COUT << "Local symbols for \"" << strName << "\":\n";
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

Symbol* NodeFor::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::vector<NodeFunction*>& vecFuncs = info.vecFuncs;

	if(!m_pExprCond) return 0;
	if(!m_pStmt) return 0;

	info.pCurLoop = this;

	Symbol *pSymInit = 0;
	if(m_pExprInit)
	{
		pSymInit = m_pExprInit->eval(info, pSym);
		safe_delete(pSymInit, pSym, info.pGlobalSyms);
	}

	while(1)
	{
		Symbol *pSymRet = 0;
		Symbol *pSymExpr = m_pExprCond->eval(info, pSym);

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

		if(m_pExprEnd)
		{
			Symbol *pSymEnd = m_pExprEnd->eval(info, pSym);
			safe_delete(pSymEnd, pSym, info.pGlobalSyms);
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

	if(m_pIdent->GetType() != NODE_IDENT)
	{
		G_CERR << linenr(T_STR"Error", info) << "Range-based for loop needs identifier."
					<< std::endl;
		return 0;
	}

	Symbol *pSymRet = 0;
	Symbol *_pArr = m_pExpr->eval(info, pSym);
	if(!_pArr)
	{
		G_CERR << linenr(T_STR"Error", info)
				<< "Invalid array for loop." << std::endl;
		return 0;
	}

	if(_pArr->GetType() != SYMBOL_ARRAY)
	{
		G_CERR << linenr(T_STR"Error", info)
				<< "Range-based for loop needs array." << std::endl;
		safe_delete(_pArr, pSym, info.pGlobalSyms);
		return 0;
	}

	SymbolArray *pArr = (SymbolArray*)_pArr;


	const t_string& strIdent = ((NodeIdent*)m_pIdent)->GetIdent();
	//std::cout << "ranged ident: " << strIdent << std::endl;

	SymbolInt *pSymIter = new SymbolInt(0);
	t_string strIter = T_STR"<cur_iter_" + strIdent + T_STR">";
	pSym->InsertSymbol(strIter, pSymIter);

	info.pCurLoop = this;
	for(unsigned int iArr=0; iArr<pArr->GetArr().size(); ++iArr)
	{
		Symbol *pSymInArr = pArr->GetArr()[iArr];
		pSym->InsertSymbol(strIdent, pSymInArr);

		Symbol *pBodyRet = m_pStmt->eval(info, pSym);
		safe_delete(pBodyRet, pSym, info.pGlobalSyms);


		// write back symbol in case an assignment has taken place
		Symbol *pNewSym = pSym->GetSymbol(strIdent);
		if(pSymInArr != pNewSym)
		{
			pArr->GetArr()[iArr] = pNewSym;
			//delete pSymInArr;
			pSymInArr = pNewSym;
		}
		pSym->RemoveSymbolNoDelete(strIdent);


		++pSymIter->GetVal();

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

	//G_COUT << "ranged for:" << pArr->GetName() << ", " << pArr->GetIdent() << std::endl;
	safe_delete(_pArr, pSym, info.pGlobalSyms);
	return 0;
}
