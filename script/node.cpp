/*
 * Script interpreter
 * @author tweber
 * @date 2013-2014
 */

#include "node.h"
#include "calls.h"
#include "helper/math.h"

ParseInfo::ParseInfo() : pmapModules(0), phandles(0),
			pmutexGlobal(0), pmutexInterpreter(0),
			pGlobalSyms(0),
			pCurFunction(0), pCurCaller(0), bWantReturn(0),
			pCurLoop(0), bWantBreak(0), bWantContinue(0),
			bDestroyParseInfo(1)
{
	pmapModules = new t_mods();
	pGlobalSyms = new SymbolTable();
	phandles = new HandleManager();
	pmutexGlobal = new std::mutex();
	pmutexInterpreter = new std::mutex();
}

ParseInfo::~ParseInfo()
{
	if(!bDestroyParseInfo)
		return;

	if(phandles) delete phandles;
	if(pGlobalSyms) delete pGlobalSyms;
	if(pmutexGlobal) delete pmutexGlobal;
	if(pmutexInterpreter) delete pmutexInterpreter;

	phandles = 0;
	pGlobalSyms = 0;

	if(pmapModules)
	{
		for(ParseInfo::t_mods::value_type vals : *pmapModules)
		{
			if(vals.second)
				delete vals.second;
		}

		pmapModules->clear();
		delete pmapModules;
		pmapModules = 0;
	}
}


template<typename T> T plus_op(T a, T b) { return a+b; }
template<typename T> T minus_op(T a, T b) { return a-b; }
template<typename T> T mult_op(T a, T b) { return a*b; }
template<typename T> T div_op(T a, T b) { return a/b; }
template<typename T> T mod_op(T a, T b) { return a%b; }
template<> t_real mod_op(t_real a, t_real b) { return ::fmod(a,b); }

template<typename T> bool log_or_op(T a, T b) { return a||b; }
template<typename T> bool log_and_op(T a, T b) { return a&&b; }
template<typename T> bool log_eq_op(T a, T b) { return a==b; }
template<> bool log_eq_op(t_real a, t_real b) { return ::float_equal(a,b); }
template<typename T> bool log_neq_op(T a, T b) { return !log_eq_op(a,b); }
template<typename T> bool log_leq_op(T a, T b) { return a<=b; }
template<typename T> bool log_geq_op(T a, T b) { return a>=b; }
template<typename T> bool log_less_op(T a, T b) { return a<b; }
template<typename T> bool log_greater_op(T a, T b) { return a>b; }

static t_int pow_int(t_int a, t_int b)
{
	return t_int(pow(a,b));
}

std::unordered_map<NodeType, t_real (*)(t_real, t_real), EnumDirectHash<NodeType>> g_mapBinOps_d = 
{
	std::pair<NodeType, t_real (*)(t_real, t_real)>(NODE_PLUS, plus_op),
	std::pair<NodeType, t_real (*)(t_real, t_real)>(NODE_MINUS, minus_op),
	std::pair<NodeType, t_real (*)(t_real, t_real)>(NODE_MULT, mult_op),
	std::pair<NodeType, t_real (*)(t_real, t_real)>(NODE_DIV, div_op),
	std::pair<NodeType, t_real (*)(t_real, t_real)>(NODE_MOD, mod_op),
	std::pair<NodeType, t_real (*)(t_real, t_real)>(NODE_POW, pow),
};

std::unordered_map<NodeType, t_int (*)(t_int, t_int), EnumDirectHash<NodeType>> g_mapBinOps_i =
{
	std::pair<NodeType, t_int (*)(t_int, t_int)>(NODE_PLUS, plus_op),
	std::pair<NodeType, t_int (*)(t_int, t_int)>(NODE_MINUS, minus_op),
	std::pair<NodeType, t_int (*)(t_int, t_int)>(NODE_MULT, mult_op),
	std::pair<NodeType, t_int (*)(t_int, t_int)>(NODE_DIV, div_op),
	std::pair<NodeType, t_int (*)(t_int, t_int)>(NODE_MOD, mod_op),
	std::pair<NodeType, t_int (*)(t_int, t_int)>(NODE_POW, pow_int)
};

std::unordered_map<NodeType, t_string (*)(t_string, t_string), EnumDirectHash<NodeType>> g_mapBinOps_s =
{
	std::pair<NodeType, t_string (*)(t_string,t_string)>(NODE_PLUS, plus_op),
};


std::unordered_map<NodeType, bool (*)(t_real, t_real), EnumDirectHash<NodeType>> g_mapBinLogOps_d =
{
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_OR, log_or_op),
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_AND, log_and_op),
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_EQ, log_eq_op),
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_NEQ, log_neq_op),
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_LEQ, log_leq_op),
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_GEQ, log_geq_op),
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_LESS, log_less_op),
	std::pair<NodeType, bool (*)(t_real, t_real)>(NODE_LOG_GREATER, log_greater_op)
};

std::unordered_map<NodeType, bool (*)(t_int, t_int), EnumDirectHash<NodeType>> g_mapBinLogOps_i =
{
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_OR, log_or_op),
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_AND, log_and_op),
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_EQ, log_eq_op),
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_NEQ, log_neq_op),
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_LEQ, log_leq_op),
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_GEQ, log_geq_op),
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_LESS, log_less_op),
	std::pair<NodeType, bool (*)(t_int, t_int)>(NODE_LOG_GREATER, log_greater_op)
};

std::unordered_map<NodeType, bool (*)(t_string,t_string), EnumDirectHash<NodeType>> g_mapBinLogOps_s =
{
	std::pair<NodeType, bool (*)(t_string,t_string)>(NODE_LOG_EQ, log_eq_op),
	std::pair<NodeType, bool (*)(t_string,t_string)>(NODE_LOG_NEQ, log_neq_op),
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

Symbol* Node::Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op)
{
	if(!pSymLeft || !pSymRight) return 0;

	const Symbol *pLeft = pSymLeft;
	const Symbol *pRight = pSymRight;
	bool bCleanLeft = 0;
	bool bCleanRight = 0;


	// map operations
	if(pSymLeft->GetType()==SYMBOL_MAP && pSymRight->GetType()==SYMBOL_MAP)
	{
		// merge maps
		if(op == NODE_PLUS)
		{
			SymbolMap *pSymMap = new SymbolMap();

			const SymbolMap::t_map& mapLeft = ((SymbolMap*)pSymLeft)->GetMap();
			const SymbolMap::t_map& mapRight = ((SymbolMap*)pSymRight)->GetMap();
			SymbolMap::t_map& mapRes = pSymMap->GetMap();

			for(const SymbolMap::t_map::value_type& pair : mapLeft)
			{
				mapRes.insert(SymbolMap::t_map::value_type(
						pair.first, pair.second->clone()));
			}

			for(const SymbolMap::t_map::value_type& pair : mapRight)
			{
				mapRes.insert(SymbolMap::t_map::value_type(
						pair.first, pair.second->clone()));
			}

			return pSymMap;
		}
	}


	// vector operations
	if(pSymLeft->GetType()==SYMBOL_ARRAY || pSymRight->GetType()==SYMBOL_ARRAY)
	{
		SymbolArray* pSymResult = new SymbolArray();

		// (vector op vector) piecewise operation
		if(pSymLeft->GetType()==SYMBOL_ARRAY && pSymRight->GetType()==SYMBOL_ARRAY)
		{
			const std::vector<Symbol*>& vecLeft = ((SymbolArray*)pSymLeft)->GetArr();
			const std::vector<Symbol*>& vecRight = ((SymbolArray*)pSymRight)->GetArr();

			unsigned int iArrSize = std::min(vecLeft.size(), vecRight.size());
			pSymResult->GetArr().reserve(iArrSize);

			for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
			{
				Symbol *pOpResult = Op(vecLeft[iElem], vecRight[iElem], op);
				pSymResult->GetArr().push_back(pOpResult);
			}
		}
		// (vector op scalar) operation
		else if(pSymLeft->GetType()==SYMBOL_ARRAY && pSymRight->GetType()!=SYMBOL_ARRAY)
		{
			const std::vector<Symbol*>& vecLeft = ((SymbolArray*)pSymLeft)->GetArr();

			unsigned int iArrSize = vecLeft.size();
			pSymResult->GetArr().reserve(iArrSize);

			for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
			{
				Symbol *pOpResult = Op(vecLeft[iElem], pSymRight, op);
				pSymResult->GetArr().push_back(pOpResult);
			}

		}
		// (scalar op vector) operation
		else if(pSymLeft->GetType()!=SYMBOL_ARRAY && pSymRight->GetType()==SYMBOL_ARRAY)
		{
			const std::vector<Symbol*>& vecRight = ((SymbolArray*)pSymRight)->GetArr();

			unsigned int iArrSize = vecRight.size();
			pSymResult->GetArr().reserve(iArrSize);

			for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
			{
				Symbol *pOpResult = Op(pSymLeft, vecRight[iElem], op);
				pSymResult->GetArr().push_back(pOpResult);
			}
		}
		return pSymResult;
	}

	// cast to equal symbol types
	if(pSymLeft->GetType() != pSymRight->GetType())
	{
		//int * string
		if(op==NODE_MULT && (pSymLeft->GetType()==SYMBOL_STRING||pSymRight->GetType()==SYMBOL_STRING))
		{
			t_int iVal = 0;
			t_string strVal;
			bool bValidIntStr = 0;

			if(pSymLeft->GetType()==SYMBOL_STRING && pSymRight->GetType()==SYMBOL_INT)
			{
				iVal = ((SymbolInt*)pSymRight)->GetVal();
				strVal = ((SymbolString*)pSymLeft)->GetVal();
				bValidIntStr = 1;
			}
			else if(pSymRight->GetType()==SYMBOL_STRING && pSymLeft->GetType()==SYMBOL_INT)
			{
				iVal = ((SymbolInt*)pSymLeft)->GetVal();
				strVal = ((SymbolString*)pSymRight)->GetVal();
				bValidIntStr = 1;
			}

			if(bValidIntStr)
			{
				SymbolString *pSymStrRet = new SymbolString();
				for(t_int iCopy=0; iCopy<iVal; ++iCopy)
					pSymStrRet->GetVal() += strVal;
				return pSymStrRet;
			}
		}

		if(pLeft->GetType()==SYMBOL_INT && pRight->GetType()==SYMBOL_DOUBLE)
		{
			pLeft = pLeft->ToType(SYMBOL_DOUBLE);
			bCleanLeft = 1;
		}

		if(pLeft->GetType()==SYMBOL_DOUBLE && pRight->GetType()==SYMBOL_INT)
		{
			pRight = pRight->ToType(SYMBOL_DOUBLE);
			bCleanRight = 1;
		}

		if(pLeft->GetType()==SYMBOL_STRING)
		{
			pRight = pRight->ToType(SYMBOL_STRING);
			bCleanRight = 1;
		}
		else if(pRight->GetType()==SYMBOL_STRING)
		{
			pLeft = pLeft->ToType(SYMBOL_STRING);
			bCleanLeft = 1;
		}
	}


	Symbol *pRes = 0;
	if(pLeft->GetType() == SYMBOL_DOUBLE)
	{
		if(IsLogicalOp(op))
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->SetName(T_STR"<op_logical_result>");

			if(g_mapBinLogOps_d.find(op) != g_mapBinLogOps_d.end())
			{
				pResult->SetVal(g_mapBinLogOps_d[op](((SymbolDouble*)pLeft)->GetVal(),
								((SymbolDouble*)pRight)->GetVal()));
			}
			else
			{
				G_CERR << "Error: Operator \"" << op 
						<< "\" not defined for real type."
						<< std::endl;
				pResult->SetVal(0);
			}
			pRes = pResult;
		}
		else
		{
			SymbolDouble *pResult = new SymbolDouble();
			pResult->SetName(T_STR"<op_result>");
			if(g_mapBinOps_d.find(op) != g_mapBinOps_d.end())
			{
				pResult->SetVal(g_mapBinOps_d[op](((SymbolDouble*)pLeft)->GetVal(),
						((SymbolDouble*)pRight)->GetVal()));
			}
			else
			{
				G_CERR << "Error: Operator \"" << op 
						<< "\" not defined for real type."
						<< std::endl;
			}
			pRes = pResult;
		}
	}
	else if(pLeft->GetType() == SYMBOL_INT)
	{
		if(IsLogicalOp(op))
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->SetName(T_STR"<op_logical_result>");
			if(g_mapBinLogOps_i.find(op) != g_mapBinLogOps_i.end())
			{
				pResult->SetVal(g_mapBinLogOps_i[op](((SymbolInt*)pLeft)->GetVal(),
								((SymbolInt*)pRight)->GetVal()));
			}
			else
			{
				G_CERR << "Error: Operator \"" << op 
						<< "\" not defined for int type."
						<< std::endl;
			}
			pRes = pResult;
		}
		else
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->SetName(T_STR"<op_result>");
			if(g_mapBinOps_i.find(op) != g_mapBinOps_i.end())
			{
				pResult->SetVal(g_mapBinOps_i[op](((SymbolInt*)pLeft)->GetVal(),
								((SymbolInt*)pRight)->GetVal()));
			}
			else
			{
				G_CERR << "Error: Operator \"" << op 
						<< "\" not defined for int type."
						<< std::endl;
			}
			pRes = pResult;
		}
	}
	else if(pLeft->GetType() == SYMBOL_STRING)
	{
		if(IsLogicalOp(op))
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->SetName(T_STR"<op_logical_result>");
			if(g_mapBinLogOps_s.find(op) != g_mapBinLogOps_s.end())
			{
				pResult->SetVal(g_mapBinLogOps_s[op](((SymbolString*)pLeft)->GetVal(),
									((SymbolString*)pRight)->GetVal()));
			}
			else
			{
				G_CERR << "Error: Operator \"" << op 
						<< "\" not defined for string type."
						<< std::endl;
			}
			pRes = pResult;
		}
		else
		{
			SymbolString *pResult = new SymbolString();
			pResult->SetName(T_STR"<op_result>");
			if(g_mapBinOps_s.find(op) != g_mapBinOps_s.end())
			{
				pResult->SetVal(g_mapBinOps_s[op](((SymbolString*)pLeft)->GetVal(),
									((SymbolString*)pRight)->GetVal()));
			}
			else
			{
				G_CERR << "Error: Operator \"" << op 
						<< "\" not defined for string type."
						<< std::endl;
			}
			pRes = pResult;
		}
	}

	if(bCleanLeft) delete pLeft;
	if(bCleanRight) delete pRight;

	return pRes;
}



NodeFunction* ParseInfo::GetFunction(const t_string& strName)
{
	for(NodeFunction* pFunc : vecFuncs)
	{
		if(pFunc && pFunc->GetName()==strName)
			return pFunc;
	}

	return 0;
}



t_string Node::linenr(const t_string& strErr, const ParseInfo &info) const
{
	t_ostringstream ostr;
	ostr << strErr;

	t_ostringstream ostrDetail;

	bool bHasLine = 0;
	if(m_iLine>0)
	{
		ostrDetail << "line " << m_iLine;
		bHasLine = 1;
	}

	if(info.pCurFunction)
	{
		if(bHasLine)
			ostrDetail << " ";

		const t_string& strFile = info.pCurFunction->GetScrFile();
		const t_string& strFkt = info.pCurFunction->GetName();

		if(strFile == T_STR"")
			ostrDetail << "in unknown file";
		else
			ostrDetail << "in \"" << strFile << "\"";

		if(strFkt != T_STR"")
			ostrDetail << " -> \"" << strFkt << "\"";
	}

	if(ostrDetail.str().length() > 0)
	{
		ostr << " (";
		ostr << ostrDetail.str();
		ostr << ")";
	}
	ostr << ": ";

	return ostr.str();
}

//--------------------------------------------------------------------------------

NodeDouble::NodeDouble(t_real dVal)
	: Node(NODE_DOUBLE), m_dVal(dVal)
{
	m_pSymbol = new SymbolDouble;
	m_pSymbol->SetVal(m_dVal);
	m_pSymbol->SetName(T_STR"<const>");
}

NodeInt::NodeInt(t_int iVal)
	: Node(NODE_INT), m_iVal(iVal)
{
	m_pSymbol = new SymbolInt;
	m_pSymbol->SetVal(m_iVal);
	m_pSymbol->SetName(T_STR"<const>");
}

NodeString::NodeString(t_string strVal)
	: Node(NODE_STRING), m_strVal(strVal)
{
	m_pSymbol = new SymbolString;
	m_pSymbol->SetVal(m_strVal);
	m_pSymbol->SetName(T_STR"<const>");
}

NodeArray::NodeArray(NodeArray* pArr)
	: Node(NODE_ARRAY), m_pArr(pArr)
{
	// don't evaluate here because these are not
	// necessarily constant and need to evaluate
	// sub-expression
}

NodeMap::NodeMap(NodeMap* pMap)
	: Node(NODE_MAP), m_pMap(pMap)
{}

NodePair::NodePair(Node* pFirst, Node* pSecond)
	: Node(NODE_PAIR), m_pFirst(pFirst), m_pSecond(pSecond)
{}

NodeRange::NodeRange(RangeType rt)
	: Node(NODE_RANGE), m_rangetype(RANGE_FULL), m_pBegin(0), m_pEnd(0)
{}

NodeRange::NodeRange(Node *pBegin, Node *pEnd)
	: Node(NODE_RANGE), m_rangetype(RANGE_BEGINEND), m_pBegin(pBegin), m_pEnd(pEnd)
{}

NodeArrayAccess::NodeArrayAccess(Node* pIdent, Node* pExpr)
	: Node(NODE_ARRAY_ACCESS), m_pIdent(pIdent), m_pExpr(pExpr)
{
	NodeBinaryOp* pIndices = (NodeBinaryOp*) m_pExpr;
	if(pIndices)
		m_vecIndices = pIndices->flatten(NODE_ARGS);
}


NodeBinaryOp::NodeBinaryOp(Node* pLeft, Node* pRight, NodeType ntype)
	: Node(ntype), m_pLeft(pLeft), m_pRight(pRight), m_bGlobal(0)
{}

void NodeBinaryOp::FlattenNodes(NodeType ntype)
{
	if(GetType() == ntype)
	{
		m_vecNodesFlat = this->flatten(ntype);
		//G_COUT << "Node type: " << ntype << ", flattened: " << m_vecNodesFlat.size() << std::endl;
	}
}


NodeCall::NodeCall(Node* _pIdent, Node* _pArgs)
	: Node(NODE_CALL), m_pIdent(_pIdent), m_pArgs(_pArgs)
{
	NodeBinaryOp* pArgs = (NodeBinaryOp*) m_pArgs;

	if(pArgs)
		m_vecArgs = pArgs->flatten(NODE_ARGS);
}


NodeFunction::NodeFunction(Node* pLeft, Node* pMiddle, Node* pRight)
	: Node(NODE_FUNC), m_pIdent(pLeft), m_pArgs(pMiddle), m_pStmts(pRight)
{
	if(m_pArgs && (m_pArgs->GetType()==NODE_IDENTS || m_pArgs->GetType()==NODE_IDENT))
	{
		NodeBinaryOp* pArgs = (NodeBinaryOp*) m_pArgs;
		m_vecArgs = pArgs->flatten(NODE_IDENTS);
	}
}

std::vector<t_string> NodeFunction::GetParamNames() const
{
	std::vector<t_string> vecNames;

	for(Node *pNode : m_vecArgs)
		vecNames.push_back(((NodeIdent*)pNode)->GetIdent());

	return vecNames;
}

//--------------------------------------------------------------------------------

Node* NodeReturn::clone() const
{
	NodeReturn* pNode = new NodeReturn(m_pExpr?m_pExpr->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeBreak::clone() const
{
	NodeBreak* pNode = new NodeBreak();
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeContinue::clone() const
{
	NodeContinue* pNode = new NodeContinue();
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeIdent::clone() const
{
	NodeIdent* pNode = new NodeIdent(m_strIdent);
	pNode->m_pDefArg = this->m_pDefArg->clone();
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeCall::clone() const
{
	NodeCall* pNode = new NodeCall(m_pIdent?m_pIdent->clone():0,
						m_pArgs?m_pArgs->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeIf::clone() const
{
	NodeIf* pNode = new NodeIf(m_pExpr?m_pExpr->clone():0,
					m_pIf?m_pIf->clone():0,
					m_pElse?m_pElse->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeWhile::clone() const
{
	NodeWhile* pNode = new NodeWhile(m_pExpr?m_pExpr->clone():0,
						m_pStmt?m_pStmt->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeFor::clone() const
{
	NodeFor* pNode = new NodeFor(m_pExprInit?m_pExprInit->clone():0,
								 m_pExprCond?m_pExprCond->clone():0,
								 m_pExprEnd?m_pExprEnd->clone():0,
								m_pStmt?m_pStmt->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeRangedFor::clone() const
{
	NodeRangedFor* pNode = new NodeRangedFor(m_pIdent?m_pIdent->clone():0,
							m_pExpr?m_pExpr->clone():0,
							m_pStmt?m_pStmt->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeDouble::clone() const
{
	NodeDouble* pNode = new NodeDouble(m_dVal);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeInt::clone() const
{
	NodeInt* pNode = new NodeInt(m_iVal);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeString::clone() const
{
	NodeString* pNode = new NodeString(m_strVal);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeArray::clone() const
{
	NodeArray* pNode = new NodeArray(m_pArr?m_pArr->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeMap::clone() const
{
	NodeMap* pNode = new NodeMap(m_pMap?m_pMap->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodePair::clone() const
{
	NodePair* pNode = new NodePair(m_pFirst?m_pFirst->clone():0, 
					m_pSecond?m_pSecond->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeRange::clone() const
{
	NodeRange *pNode = new NodeRange(m_rangetype);
	if(m_pBegin) pNode->m_pBegin = m_pBegin->clone();
	if(m_pEnd) pNode->m_pEnd = m_pEnd->clone();

	return pNode;
}

Node* NodeArrayAccess::clone() const
{
	NodeArrayAccess *pNode = new NodeArrayAccess(m_pIdent?m_pIdent->clone():0,
										m_pExpr?m_pExpr->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeUnaryOp::clone() const
{
	NodeUnaryOp* pNode = new NodeUnaryOp(m_pChild?m_pChild->clone():0, GetType());
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeBinaryOp::clone() const
{
	NodeBinaryOp* pNode = new NodeBinaryOp(m_pLeft?m_pLeft->clone():0,
							m_pRight?m_pRight->clone():0,
							GetType());
	pNode->m_vecNodesFlat = this->m_vecNodesFlat;
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeFunction::clone() const
{
	NodeFunction* pFkt = new NodeFunction(m_pIdent?m_pIdent->clone():0,
							m_pArgs?m_pArgs->clone():0,
							m_pStmts?m_pStmts->clone():0);
	pFkt->m_strScrFile = this->m_strScrFile;

	*((Node*)pFkt) = *((Node*)this);
	return pFkt;
}


//--------------------------------------------------------------------------------

// TODO: Fix: Gets called for non casted NodeBinaryOps which are of other type!
std::vector<Node*> NodeBinaryOp::flatten(NodeType ntype) const
{
	//G_COUT << GetType() << ", " << ntype << std::endl;

	std::vector<Node*> vecNodes;

	if(GetType() == ntype)
	{
		NodeBinaryOp *pLeft = (NodeBinaryOp*) m_pLeft;
		NodeBinaryOp *pRight = (NodeBinaryOp*) m_pRight;

		if(m_pLeft)
		{
			if(m_pLeft->GetType() == ntype)
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
			if(m_pRight->GetType() == ntype)
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


const t_string& NodeFunction::GetName() const
{
	static t_string strNull = T_STR"<null>";
	if(!m_pIdent || m_pIdent->GetType() != NODE_IDENT)
		return strNull;

	NodeIdent *pIdent = (NodeIdent*)m_pIdent;
	return pIdent->GetIdent();
}
