/*
 * Script interpreter
 * Node evaluation
 * @author tweber
 * @date 10 oct 2013
 */

#include "node.h"
#include "calls.h"
#include "helper/log.h"


Symbol* NodeReturn::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;
	Symbol *pRet = 0;

	if(m_pExpr)
	{
		Symbol *pEval = m_pExpr->eval(info, pSym);
		if(clone_if_needed(pEval, pRet))
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
		std::ostringstream ostrErr;
		ostrErr << linenr(info)
				<< "Cannot use break outside loop."
				<< std::endl;
		throw Err(ostrErr.str(),0);
	}

	info.bWantBreak = 1;
	return 0;
}

Symbol* NodeContinue::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	if(info.pCurLoop==0)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(info)
				<< "Cannot use continue outside loop."
				<< std::endl;
		throw Err(ostrErr.str(),0);
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
		log_warn(linenr(info), "Symbol \"", m_strIdent,
			"\" exists in local and global scope, using local one.");
	}

	// if no local symbol is available, use global symbol instead
	if(!pSymbol)
		pSymbol = pSymbolGlob;

	if(!pSymbol)
	{
		// no throw, so null symbol or symbol validity query can be used in script
		log_err(linenr(info), "Symbol \"", m_strIdent,
			"\" not in symbol table.");
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
				std::ostringstream ostrErr;
				ostrErr << linenr(info)
					<< "Invalid symbol to unpack." << std::endl;
				throw Err(ostrErr.str(),0);
			}
			Symbol *pSymChild = pChild->eval(info, pSym);
			if(pSymChild->GetType() != SYMBOL_MAP)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(info)
					<< "Cannot unpack non-map." << std::endl;
				throw Err(ostrErr.str(),0);
			}

			SymbolMap::t_map& mapSyms = ((SymbolMap*)pSymChild)->GetMap();

			std::vector<t_string> vecParamNames = pFkt->GetParamNames();
			for(unsigned int iParam=vecArgSyms.size(); iParam<vecParamNames.size(); ++iParam)
			{
				const t_string& strParamName = vecParamNames[iParam];
				SymbolMap::t_map::iterator iter = mapSyms.find(SymbolMapKey(strParamName));
				if(iter == mapSyms.end())
				{
					log_err(linenr(info), "Parameter \"", strParamName,
						"\" not in map. Using 0.");

					vecArgSyms.push_back(new SymbolDouble(0.));
					continue;
				}

				Symbol *pSymClone;
				clone_if_needed(iter->second, pSymClone);
				vecArgSyms.push_back(pSymClone);
			}

			safe_delete(pSymChild, pSym, info.pGlobalSyms);
		}
		else
		{
			//log_debug("type: ", pNode->GetType());
			Symbol *pSymbol = pNode->eval(info, pSym);
			vecArgSyms.push_back(pSymbol);
		}
	}

	arrArgs.UpdateIndices();
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
		if(info.bEnableDebug)
		{
			std::string strTrace = "syscall: " + strFkt + ", " 
					+ std::to_string(vecArgSyms.size()) + " args";
			info.PushTraceback(std::move(strTrace));
		}

		pFktRet = ext_call(strFkt, vecArgSyms, info, pSym);

		if(info.bEnableDebug)
		{
			info.PopTraceback();
		}
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

	std::ostringstream ostrErr;
	ostrErr << linenr(info)
		<< "Pairs should not be evaluated directly." 
		<< std::endl;
	throw Err(ostrErr.str(),0);
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
			std::ostringstream ostrErr;
			ostrErr << linenr(info)
				<< "Maps have to consist of key-value pairs." 
				<< std::endl;
			throw Err(ostrErr.str(),0);
			//continue;
		}

		Symbol* pSymFirst = 0;
		Symbol* pSymSecond = 0;

		if(((NodePair*)pNode)->GetFirst())
			pSymFirst = ((NodePair*)pNode)->GetFirst()->eval(info, pSym);
		if(((NodePair*)pNode)->GetSecond())
			pSymSecond = ((NodePair*)pNode)->GetSecond()->eval(info, pSym);

		bool bSecondCloned = 0;
		if(pSymFirst && pSymSecond)
		{
			Symbol *pSymClone;
			bSecondCloned = clone_if_needed(pSymSecond, pSymClone);

			pSymMap->GetMap().insert(
				SymbolMap::t_map::value_type(SymbolMapKey(pSymFirst),
				pSymClone));
		}

		safe_delete(pSymFirst, pSym, info.pGlobalSyms);
		if(bSecondCloned)
			safe_delete(pSymSecond, pSym, info.pGlobalSyms);
	}

	pSymMap->UpdateIndices();
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
		Symbol *pClone;
		if(clone_if_needed(pSymbol, pClone))
			safe_delete(pSymbol, pSym, info.pGlobalSyms);

		pSymArr->GetArr().push_back(pClone);
		pSymArr->UpdateLastNIndices(1);
	}
	return pSymArr;
}

Symbol* NodeRange::eval(ParseInfo &info, SymbolTable *pSym) const
{
	if(info.IsExecDisabled()) return 0;

	std::ostringstream ostrErr;
	ostrErr << linenr(info) << "Range should not be evaluated directly." << std::endl;
	throw Err(ostrErr.str(), 0);
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
			std::ostringstream ostrErr;
			ostrErr << linenr(info) << "Invalid range." << std::endl;
			throw Err(ostrErr.str(),0);
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
			log_err(linenr(info), "Lower array index out of bounds: ",
				iBeginIdx, ". Adjusting to lower limit.");
			iBeginIdx = 0;
		}
		if(iEndIdx<-1 || iEndIdx>iMaxLen)
		{
			log_err(linenr(info), "Upper array index out of bounds: ",
				iEndIdx, ". Adjusting to upper limit.");
			iEndIdx = iMaxLen;
		}

		//G_COUT << "begin: " << iBeginIdx << ", end: " << iEndIdx << std::endl;

		safe_delete(pSymBeg, pSym, info.pGlobalSyms);
		safe_delete(pSymEnd, pSym, info.pGlobalSyms);
	}
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(info) << "Invalid range operation." << std::endl;
		throw Err(ostrErr.str(),0);
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
		std::ostringstream ostrErr;
		ostrErr << linenr(info) << "Tried to access non-array." << std::endl;
		throw Err(ostrErr.str(),0);
	}

	t_string strIdent;
	Symbol* pSymbol = m_pIdent->eval(info, pSym);

	if(!pSymbol)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(info) 
			<< "Symbol for array not found." << std::endl;
		throw Err(ostrErr.str(),0);
	}

	if(m_pIdent->GetType() == NODE_IDENT)
		strIdent = pSymbol->GetIdent();
	else
		strIdent = T_STR"<tmp_sym>";


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
				safe_delete(pSymbol, pSym, info.pGlobalSyms);
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
					std::ostringstream ostrErr;
					ostrErr << linenr(info)
							<< "Array index has to be of integer type."
							<< std::endl;
					throw Err(ostrErr.str(),0);
				}

				t_int iIdx = pSymExpr->GetValInt();
				safe_delete(pSymExpr, pSym, info.pGlobalSyms);

				// convert negative indices
				if(iIdx < 0)
					iIdx = pArr->GetArr().size()  + iIdx;

				if(iIdx < 0)
				{
					std::ostringstream ostrErr;
					ostrErr << linenr(info)
						<< "Invalid array index."
						<< std::endl;
					throw Err(ostrErr.str(),0);
				}

				// index too high -> fill up with zeroes
				if(iIdx>=pArr->GetArr().size())
				{
					unsigned int iOldSize = pArr->GetArr().size();
					for(unsigned int iRem=0; iRem<iIdx+1-iOldSize; ++iRem)
					{
						SymbolDouble *pNewSym = new SymbolDouble(0.);
						pNewSym->SetConst(1);
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
			std::ostringstream ostrErr;
			ostrErr << linenr(info) << "No key given for map." << std::endl;
			throw Err(ostrErr.str(),0);
		}
		else if(m_vecIndices.size()>1)
		{
			log_warn(linenr(info), "Multiple keys given for map, using first one.");
		}

		Node *pNodeKey = m_vecIndices[0];
		Symbol *pSymExpr = pNodeKey->eval(info, pSym);
		if(pSymExpr==0)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(info) << "Map key is invalid."
						<< std::endl;
			throw Err(ostrErr.str(),0);
		}

		SymbolMap *pMap = (SymbolMap*)pSymbol;
		SymbolMapKey key = SymbolMapKey(pSymExpr);
		SymbolMap::t_map::iterator iterMap = pMap->GetMap().find(key);

		// key not yet in map -> insert it
		if(iterMap == pMap->GetMap().end())
		{
			SymbolString *pNewSym = new SymbolString();
			pNewSym->SetConst(1);
			iterMap = pMap->GetMap().insert(SymbolMap::t_map::value_type(key, pNewSym)).first;
		}

		if((void*)pSymbol != (void*)iterMap->second)
			safe_delete(pSymbol, pSym, info.pGlobalSyms);

		pSymbol = iterMap->second;
		pMap->UpdateIndex(key);
		safe_delete(pSymExpr, pSym, info.pGlobalSyms);
	}
	else if(pSymbol->GetType() == SYMBOL_STRING)
	{
		SymbolString *pStr = (SymbolString*)pSymbol;
		t_int iStrLen = pStr->GetVal().length();

		if(m_vecIndices.size()!=1)
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(info)
					<< "Need exactly one string index."
					<< std::endl;
			throw(Err(ostrErr.str()),0);
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
				std::ostringstream ostrErr;
				ostrErr << linenr(info)
						<< "String index has to be of integer type."
						<< std::endl;
				throw Err(ostrErr.str(),0);
			}

			t_int iIdx = pSymExpr->GetValInt();
			safe_delete(pSymExpr, pSym, info.pGlobalSyms);

			// convert negative indices
			if(iIdx < 0)
				iIdx = iStrLen  + iIdx;

			if(iIdx < 0 || iIdx >= iStrLen)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(info) << "String index out of bounds."
					<< std::endl;

				throw Err(ostrErr.str(),0);
			}

			t_string strNew;
			strNew += pStr->GetVal()[iIdx];
			safe_delete(pSymbol, pSym, info.pGlobalSyms);
			pSymbol = new SymbolString(strNew);
		}
	}
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(info) << "Symbol \"" << strIdent
				<< "\" is neither an array nor a map." << std::endl;
		throw Err(ostrErr.str(),0);
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
		log_err(linenr(info), "Unary minus not defined for ", pSym->GetTypeName(), ".");
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
			Symbol *pClone;
			if(clone_if_needed(pSymbolEval, pClone))
				safe_delete(pSymbolEval, pSym, info.pGlobalSyms);

			uminus_inplace(pClone, info);
			return pClone;
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
			log_warn(linenr(info), "Unpack operation only allowed in function call. Ignoring.");

			if(m_pChild)
				return m_pChild->eval(info, pSym);
		}
		default:
			log_warn(linenr(info), "Unknown node type: ", GetType());
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
		log_err(linenr(info), "NULL assignment.");
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
		log_err(linenr(info), "Invalid rhs expression in assignment.");
		return 0;
	}

	Symbol *pSymbol;
	if(clone_if_needed(pSymbolOrg, pSymbol))
		safe_delete(pSymbolOrg, pSym, info.pGlobalSyms);
	pSymbol->SetRval(0);

	if(pLeft->GetType() == NODE_IDENT)		// single variable
	{
		const t_string& strIdent = ((NodeIdent*)pLeft)->GetIdent();
		//log_debug("Assigning ", strIdent, " = ", pSymbol);

		Symbol* pSymGlob = info.pGlobalSyms->GetSymbol(strIdent);
		Symbol* pSymLoc = 0;
		if(!*pbGlob)
			pSymLoc = pSym->GetSymbol(strIdent);

		if(pSymLoc && pSymGlob)
		{
			log_warn(linenr(info), "Symbol \"", strIdent,
				"\" exists in local and global scope, using local one.");
		}

		if(pSymGlob && !pSymLoc && !*pbGlob)
		{
			log_warn(linenr(info), "Overwriting global symbol \"", strIdent, "\".");
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
			log_err(linenr(info), "Assignment needs an array on the right-hand side.");
			return pSymbol;
		}
		SymbolArray *pArrRight = (SymbolArray*)pSymbol;

		std::vector<Node*> vecLeftArgs = ((NodeBinaryOp*)((NodeArray*)pLeft)->GetArr())->flatten(NODE_ARGS);
		unsigned int iArrSize = vecLeftArgs.size();
		if(vecLeftArgs.size() != pArrRight->GetArr().size())
		{
                        log_warn(linenr(info),
				"Size mismatch between assigned and returned array: ",
				vecLeftArgs.size(), " != ", pArrRight->GetArr().size(), ".");

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
			std::ostringstream ostrErr;
			ostrErr << linenr(info)
				<< "No array element found." << std::endl;
			throw Err(ostrErr.str(),0);
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
					pNewSym->SetConst(1);
					pArr->GetArr().push_back(pNewSym);
					pArr->UpdateLastNIndices(1);
				}
			}


			Symbol* pSymOld = pArr->GetArr()[iArrIdx];
			if((void*)pSymOld != (void*)pSymLeft)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(info)
					<< "Array member mismatch." << std::endl;
				throw Err(ostrErr.str(),0);
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
			const SymbolMapKey& MapKey = pSymLeft->GetMapKey();
			SymbolMap* pMap = pSymLeft->GetMapPtr();

			SymbolMap::t_map::iterator iterMap = pMap->GetMap().find(MapKey);

			// symbol not in map -> insert a zero
			if(iterMap == pMap->GetMap().end())
			{
				SymbolDouble *pNewSym = new SymbolDouble(0.);
				pNewSym->SetConst(1);
				iterMap = pMap->GetMap().insert(SymbolMap::t_map::value_type(MapKey, pNewSym)).first;
			}

			Symbol* pSymOld = iterMap->second;
			if((void*)pSymOld != (void*)pSymLeft)
			{
				std::ostringstream ostrErr;
				ostrErr << linenr(info) << "Map member mismatch." << std::endl;
				throw Err(ostrErr.str(),0);
			}

			pSymbol->SetMapPtr(pMap);
			pSymbol->SetMapKey(MapKey);
			iterMap->second = pSymbol;

			pSymOld->SetMapPtr(0);
			safe_delete(pSymOld, pSym, info.pGlobalSyms);
		}
		else
		{
			std::ostringstream ostrErr;
			ostrErr << linenr(info)
				<< "Trying to access array/map member with no associated array/map."
				<< std::endl;
			throw Err(ostrErr.str(),0);
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
				log_warn(linenr(info),
					"Function \"", strFktName,
					"\" redefined in \"", info.strInitScrFile, "\".",
					" Ignoring.");
		}
		else
		{
			if(has_ext_call(strFktName))
			{
				log_warn(linenr(info),
					"Function \"", strFktName,
					"\" in \"", info.strInitScrFile,
					"\" overwrites a system function.");
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

		log_err(linenr(info), "Function \"", info.strExecFkt, "\" not defined.");
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
				log_info(linenr(info), "Should better be evaluated sequentially in interpreter.");
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

		if(pSymbol!=pSymbolLeft) safe_delete(pSymbolLeft, pSym, info.pGlobalSyms);
		if(pSymbol!=pSymbolRight) safe_delete(pSymbolRight, pSym, info.pGlobalSyms);

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

	std::unique_ptr<SymbolTable> ptrLocalSym(new SymbolTable);
	SymbolTable *pLocalSym = ptrLocalSym.get();

	SymbolArray* pArgs = 0;
	if(pTableSup)
		pArgs = (SymbolArray*)pTableSup->GetSymbol(T_STR"<args>");

	unsigned int iArgSize = 0;
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

		iArgSize = /*std::min(*/m_vecArgs.size()/*, pVecArgSyms->size())*/;
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
					log_warn(linenr(info),
						"Given argument \"", pIdent->GetIdent(), 
						"\" for function \"", 
						strName, "\" not valid. ", 
						"Using default argument.");
			}
			if(pSymbol==0)
			{
				log_err(linenr(info), "Argument \"",
					pIdent->GetIdent(), "\" for function \"",
					strName, "\" not given. Ignoring.");
				continue;
			}

			//G_COUT << "arg: " << pIdent->GetIdent() << std::endl;
			Symbol* pSymToInsert;
			if(clone_if_needed(pSymbol, pSymToInsert))
				/*safe_delete(pSymbol, pLocalSym, info.pGlobalSyms)*/;
			pSymToInsert->SetRval(0);
			pLocalSym->InsertSymbol(pIdent->GetIdent(), pSymToInsert);
		}
	}


	if(info.bEnableDebug)
	{
		std::string strTrace = "call: " + GetName() + ", " 
					+ std::to_string(iArgSize) + " args";
		info.PushTraceback(std::move(strTrace));
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

	if(info.bEnableDebug)
	{
		info.PopTraceback();
	}

	//G_COUT << "Local symbols for \"" << strName << "\":\n";
	//pLocalSym->print();

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
		std::ostringstream ostrErr;
		ostrErr << linenr(info) 
			<< "Range-based for loop needs identifier." << std::endl;
		throw Err(ostrErr.str(),0);
	}

	Symbol *pSymRet = 0;
	Symbol *_pArr = m_pExpr->eval(info, pSym);
	if(!_pArr)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(info)
			<< "Invalid array for loop." << std::endl;
		return 0;
	}

	if(_pArr->GetType() != SYMBOL_ARRAY)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(info)
				<< "Range-based for loop needs array." << std::endl;
		safe_delete(_pArr, pSym, info.pGlobalSyms);

		throw Err(ostrErr.str(),0);
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
