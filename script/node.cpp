/*
 * Script interpreter
 * @author tweber
 */

#include "node.h"
#include "calls.h"

ParseInfo::ParseInfo() : pmapModules(0), phandles(0),
			pmutexGlobal(0), pmutexInterpreter(0),
			pvecExecArg(0), pGlobalSyms(0),
			pCurFunction(0), pCurCaller(0), bWantReturn(0),
			pCurLoop(0), bWantBreak(0), bWantContinue(0),
			bDestroyParseInfo(1)
{
	pmapModules = new t_mods();
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
template<> double mod_op(double a, double b) { return ::fmod(a,b); }

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
	std::pair<NodeType, double (*)(double,double)>(NODE_MOD, mod_op),
	std::pair<NodeType, double (*)(double,double)>(NODE_POW, pow),
};

std::map<NodeType, int (*)(int,int)> g_mapBinOps_i =
{
	std::pair<NodeType, int (*)(int,int)>(NODE_PLUS, plus_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_MINUS, minus_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_MULT, mult_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_DIV, div_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_MOD, mod_op),
	std::pair<NodeType, int (*)(int,int)>(NODE_POW, pow_int)
};

std::map<NodeType, std::string (*)(std::string, std::string)> g_mapBinOps_s =
{
	std::pair<NodeType, std::string (*)(std::string,std::string)>(NODE_PLUS, plus_op),
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

std::map<NodeType, bool (*)(std::string,std::string)> g_mapBinLogOps_s =
{
	std::pair<NodeType, bool (*)(std::string,std::string)>(NODE_LOG_EQ, log_eq_op),
	std::pair<NodeType, bool (*)(std::string,std::string)>(NODE_LOG_NEQ, log_neq_op),
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

	// vector operations
	if(pSymLeft->GetType()==SYMBOL_ARRAY || pSymRight->GetType()==SYMBOL_ARRAY)
	{
		SymbolArray* pSymResult = new SymbolArray();

		// (vector op vector) piecewise operation
		if(pSymLeft->GetType()==SYMBOL_ARRAY && pSymRight->GetType()==SYMBOL_ARRAY)
		{
			const std::vector<Symbol*>& vecLeft = ((SymbolArray*)pSymLeft)->m_arr;
			const std::vector<Symbol*>& vecRight = ((SymbolArray*)pSymRight)->m_arr;

			unsigned int iArrSize = std::min(vecLeft.size(), vecRight.size());
			pSymResult->m_arr.reserve(iArrSize);

			for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
			{
				Symbol *pOpResult = Op(vecLeft[iElem], vecRight[iElem], op);
				pSymResult->m_arr.push_back(pOpResult);
			}
		}
		// (vector op scalar) operation
		else if(pSymLeft->GetType()==SYMBOL_ARRAY && pSymRight->GetType()!=SYMBOL_ARRAY)
		{
			const std::vector<Symbol*>& vecLeft = ((SymbolArray*)pSymLeft)->m_arr;

			unsigned int iArrSize = vecLeft.size();
			pSymResult->m_arr.reserve(iArrSize);

			for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
			{
				Symbol *pOpResult = Op(vecLeft[iElem], pSymRight, op);
				pSymResult->m_arr.push_back(pOpResult);
			}

		}
		// (scalar op vector) operation
		else if(pSymLeft->GetType()!=SYMBOL_ARRAY && pSymRight->GetType()==SYMBOL_ARRAY)
		{
			const std::vector<Symbol*>& vecRight = ((SymbolArray*)pSymRight)->m_arr;

			unsigned int iArrSize = vecRight.size();
			pSymResult->m_arr.reserve(iArrSize);

			for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
			{
				Symbol *pOpResult = Op(pSymLeft, vecRight[iElem], op);
				pSymResult->m_arr.push_back(pOpResult);
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
			int iVal = 0;
			std::string strVal;
			bool bValidIntStr = 0;

			if(pSymLeft->GetType()==SYMBOL_STRING && pSymRight->GetType()==SYMBOL_INT)
			{
				iVal = ((SymbolInt*)pSymRight)->m_iVal;
				strVal = ((SymbolString*)pSymLeft)->m_strVal;
				bValidIntStr = 1;
			}
			else if(pSymRight->GetType()==SYMBOL_STRING && pSymLeft->GetType()==SYMBOL_INT)
			{
				iVal = ((SymbolInt*)pSymLeft)->m_iVal;
				strVal = ((SymbolString*)pSymRight)->m_strVal;
				bValidIntStr = 1;
			}

			if(bValidIntStr)
			{
				SymbolString *pSymStrRet = new SymbolString();
				for(int iCopy=0; iCopy<iVal; ++iCopy)
					pSymStrRet->m_strVal += strVal;
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
			pResult->m_strName = "<op_logical_result>";
			
			if(g_mapBinLogOps_d.find(op) != g_mapBinLogOps_d.end())
			{
				pResult->m_iVal = g_mapBinLogOps_d[op](((SymbolDouble*)pLeft)->m_dVal,
													((SymbolDouble*)pRight)->m_dVal);
			}
			else
			{
				std::cerr << "Error: Operator \"" << op << "\" not defined for double type."
						 << std::endl;
				pResult->m_iVal = 0;		 
			}
			pRes = pResult;
		}
		else
		{
			SymbolDouble *pResult = new SymbolDouble();
			pResult->m_strName = "<op_result>";
			if(g_mapBinOps_d.find(op) != g_mapBinOps_d.end())
			{
				pResult->m_dVal = g_mapBinOps_d[op](((SymbolDouble*)pLeft)->m_dVal,
												((SymbolDouble*)pRight)->m_dVal);
			}
			else
			{
				std::cerr << "Error: Operator \"" << op << "\" not defined for double type."
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
			pResult->m_strName = "<op_logical_result>";
			if(g_mapBinLogOps_i.find(op) != g_mapBinLogOps_i.end())
			{
				pResult->m_iVal = g_mapBinLogOps_i[op](((SymbolInt*)pLeft)->m_iVal,
													((SymbolInt*)pRight)->m_iVal);
			}
			else
			{
				std::cerr << "Error: Operator \"" << op << "\" not defined for int type."
						 << std::endl;
			}
			pRes = pResult;
		}
		else
		{
			SymbolInt *pResult = new SymbolInt();
			pResult->m_strName = "<op_result>";
			if(g_mapBinOps_i.find(op) != g_mapBinOps_i.end())
			{
				pResult->m_iVal = g_mapBinOps_i[op](((SymbolInt*)pLeft)->m_iVal,
												((SymbolInt*)pRight)->m_iVal);
			}
			else
			{
				std::cerr << "Error: Operator \"" << op << "\" not defined for int type."
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
			pResult->m_strName = "<op_logical_result>";
			if(g_mapBinLogOps_s.find(op) != g_mapBinLogOps_s.end())
			{
				pResult->m_iVal = g_mapBinLogOps_s[op](((SymbolString*)pLeft)->m_strVal,
													((SymbolString*)pRight)->m_strVal);
			}
			else
			{
				std::cerr << "Error: Operator \"" << op << "\" not defined for string type."
						 << std::endl;
			}
			pRes = pResult;
		}
		else
		{
			SymbolString *pResult = new SymbolString();
			pResult->m_strName = "<op_result>";
			if(g_mapBinOps_s.find(op) != g_mapBinOps_s.end())
			{
				pResult->m_strVal = g_mapBinOps_s[op](((SymbolString*)pLeft)->m_strVal,
												((SymbolString*)pRight)->m_strVal);
			}
			else
			{
				std::cerr << "Error: Operator \"" << op << "\" not defined for string type."
						 << std::endl;
			}
			pRes = pResult;
		}
	}
	
	if(bCleanLeft) delete pLeft;
	if(bCleanRight) delete pRight;

	return pRes;
}



NodeFunction* ParseInfo::GetFunction(const std::string& strName)
{
	for(NodeFunction* pFunc : vecFuncs)
	{
		if(pFunc && pFunc->GetName()==strName)
			return pFunc;
	}

	return 0;
}



std::string Node::linenr(const std::string& strErr, const ParseInfo &info) const
{
	std::ostringstream ostr;

	ostr << strErr << " ";
	ostr << "(";

	if(m_iLine==0)
		ostr << "unknown line";
	else
		ostr << "line " << m_iLine;

	if(info.pCurFunction)
	{
		const std::string& strFile = info.pCurFunction->m_strScrFile;
		std::string strFkt = info.pCurFunction->GetName();

		if(strFile == "")
			ostr << " in unknown file";
		else
			ostr << " in \"" << strFile << "\"";

		if(strFkt != "")
			ostr << " -> \"" << strFkt << "\"";
	}

	ostr << "): ";

	return ostr.str();
}

//--------------------------------------------------------------------------------

NodeDouble::NodeDouble(double dVal)
	: Node(NODE_DOUBLE), m_dVal(dVal)
{
	m_pSymbol = new SymbolDouble;
	m_pSymbol->m_dVal = m_dVal;
	m_pSymbol->m_strName = "<const>";
}

NodeInt::NodeInt(int iVal)
	: Node(NODE_INT), m_iVal(iVal)
{
	m_pSymbol = new SymbolInt;
	m_pSymbol->m_iVal = m_iVal;
	m_pSymbol->m_strName = "<const>";
}

NodeString::NodeString(std::string strVal)
	: Node(NODE_STRING), m_strVal(strVal)
{
	m_pSymbol = new SymbolString;
	m_pSymbol->m_strVal = m_strVal;
	m_pSymbol->m_strName = "<const>";
}

NodeArray::NodeArray(NodeArray* pArr)
	: Node(NODE_ARRAY), m_pArr(pArr)
{
	// don't evaluate here because these are not
	// necessarily constant and need to evaluate
	// sub-expression
}

NodeArrayAccess::NodeArrayAccess(Node* pIdent, Node* pExpr)
	: Node(NODE_ARRAY_ACCESS), m_pIdent(pIdent), m_pExpr(pExpr)
{
	NodeBinaryOp* pIndices = (NodeBinaryOp*) m_pExpr;
	if(pIndices)
		m_vecIndices = pIndices->flatten(NODE_ARGS);
}


NodeCall::NodeCall(Node* _pIdent, Node* _pArgs)
	: Node(NODE_CALL), m_pIdent(_pIdent), m_pArgs(_pArgs)
{
	NodeBinaryOp* pArgs = (NodeBinaryOp*) m_pArgs;
	
	if(pArgs)
		m_vecArgs = pArgs->flatten(NODE_ARGS);
}


NodeFunction::NodeFunction(Node* pLeft, Node* pMiddle, Node* pRight)
	: Node(NODE_FUNC), m_pIdent(pLeft), m_pArgs(pMiddle), m_pStmts(pRight),
		m_pVecArgSyms(0)
{
	if(m_pArgs && (m_pArgs->m_type==NODE_IDENTS || m_pArgs->m_type==NODE_IDENT))
	{
		NodeBinaryOp* pArgs = (NodeBinaryOp*) m_pArgs;
		m_vecArgs = pArgs->flatten(NODE_IDENTS);
	}
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

Node* NodeArrayAccess::clone() const
{
	NodeArrayAccess *pNode = new NodeArrayAccess(m_pIdent?m_pIdent->clone():0,
										m_pExpr?m_pExpr->clone():0);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeUnaryOp::clone() const
{
	NodeUnaryOp* pNode = new NodeUnaryOp(m_pChild?m_pChild->clone():0,
										m_type);
	*((Node*)pNode) = *((Node*)this);
	return pNode;
}

Node* NodeBinaryOp::clone() const
{
	NodeBinaryOp* pNode = new NodeBinaryOp(m_pLeft?m_pLeft->clone():0,
							m_pRight?m_pRight->clone():0,
							m_type);
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
	//std::cout << m_type << ", " << ntype << std::endl;

	std::vector<Node*> vecNodes;

	if(m_type==ntype)
	{
		NodeBinaryOp *pLeft = (NodeBinaryOp*) m_pLeft;
		NodeBinaryOp *pRight = (NodeBinaryOp*) m_pRight;
		
		if(m_pLeft)
		{
			if(m_pLeft->m_type==ntype)
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
			if(m_pRight->m_type==ntype)
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
