/*
 * Script interpreter
 * @author tweber
 */

#include "node.h"
#include "calls.h"

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

	if(pLeft->GetType() == SYMBOL_DOUBLE)
	{
		Symbol *pRes = 0;
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

		return pRes;
	}
	else if(pLeft->GetType() == SYMBOL_INT)
	{
		Symbol *pRes = 0;
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

		return pRes;
	}
	else if(pLeft->GetType() == SYMBOL_STRING)
	{
		Symbol *pRes = 0;
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

		return pRes;
	}
	
	if(bCleanLeft) delete pLeft;
	if(bCleanRight) delete pRight;

	return 0;
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
	return new NodeReturn(m_pExpr?m_pExpr->clone():0);
}

Node* NodeIdent::clone() const
{
	return new NodeIdent(m_strIdent);
}

Node* NodeCall::clone() const
{
	return new NodeCall(m_pIdent?m_pIdent->clone():0,
						m_pArgs?m_pArgs->clone():0);
}

Node* NodeIf::clone() const
{
	return new NodeIf(m_pExpr?m_pExpr->clone():0,
					m_pIf?m_pIf->clone():0,
					m_pElse?m_pElse->clone():0);
}

Node* NodeWhile::clone() const
{
	return new NodeWhile(m_pExpr?m_pExpr->clone():0,
						m_pStmt?m_pStmt->clone():0);
}

Node* NodeRangedFor::clone() const
{
	return new NodeRangedFor(m_pIdent?m_pIdent->clone():0,
							m_pExpr?m_pExpr->clone():0,
							m_pStmt?m_pStmt->clone():0);
}

Node* NodeDouble::clone() const
{
	return new NodeDouble(m_dVal);
}

Node* NodeInt::clone() const
{
	return new NodeInt(m_iVal);
}

Node* NodeString::clone() const
{
	return new NodeString(m_strVal);
}

Node* NodeArray::clone() const
{
	return new NodeArray(m_pArr?m_pArr->clone():0);
}

Node* NodeArrayAccess::clone() const
{
	return new NodeArrayAccess(m_pIdent?m_pIdent->clone():0,
							m_pExpr?m_pExpr->clone():0);
}

Node* NodeUnaryOp::clone() const
{
	return new NodeUnaryOp(m_pChild?m_pChild->clone():0,
							m_type);
}

Node* NodeBinaryOp::clone() const
{
	return new NodeBinaryOp(m_pLeft?m_pLeft->clone():0,
							m_pRight?m_pRight->clone():0,
							m_type);
}

Node* NodeFunction::clone() const
{
	return new NodeFunction(m_pIdent?m_pIdent->clone():0,
							m_pArgs?m_pArgs->clone():0,
							m_pStmts?m_pStmts->clone():0);
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
