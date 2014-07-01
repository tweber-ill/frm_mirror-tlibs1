/*
 * Script interpreter
 * @author tweber
 * @date 2013-2014
 */

#ifndef __MIEZE_NODE__
#define __MIEZE_NODE__

#include "types.h"

#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <mutex>

#include "symbol.h"
#include "handles.h"

struct Node;
struct NodeFunction;
struct NodeCall;

struct ParseInfo
{
	// external imported modules
	typedef std::map<t_string, Node*> t_mods;
	t_mods *pmapModules;

	// function to execute, e.g. "main"
	t_string strExecFkt;
	t_string strInitScrFile;

	// all functions from all modules
	std::vector<NodeFunction*> vecFuncs;

	// global symbol table
	SymbolTable *pGlobalSyms;

	HandleManager *phandles;
	// mutex for script
	std::mutex *pmutexGlobal;
	// mutex for interpreter
	std::mutex *pmutexInterpreter;

	// currently active function
	const NodeFunction *pCurFunction;
	const NodeCall *pCurCaller;
	bool bWantReturn = 0;

	const Node* pCurLoop;
	bool bWantBreak;
	bool bWantContinue;

	bool bDestroyParseInfo;

	ParseInfo();
	virtual ~ParseInfo();

	NodeFunction* GetFunction(const t_string& strName);

	bool IsExecDisabled() const
	{
		return bWantReturn || bWantBreak || bWantContinue;
	}

protected:
	void init_global_syms();
};


enum NodeType
{
	NODE_NOP,

	NODE_PLUS,
	NODE_MINUS,
	NODE_DIV,
	NODE_MOD,
	NODE_MULT,
	NODE_POW,

	NODE_UMINUS,

	NODE_ASSIGN,
	NODE_CALL,

	NODE_STMTS,

	NODE_DOUBLE,
	NODE_INT,
	NODE_STRING,

	NODE_ARRAY,
	NODE_ARRAY_ACCESS,

	NODE_MAP,
	NODE_PAIR,

	NODE_IDENTS,
	NODE_IDENT,

	NODE_RANGE,
	NODE_ARGS,

	NODE_FUNCS,
	NODE_FUNC,

	NODE_LOG_AND,
	NODE_LOG_OR,
	NODE_LOG_NOT,
	NODE_LOG_EQ,
	NODE_LOG_NEQ,
	NODE_LOG_LESS,
	NODE_LOG_GREATER,
	NODE_LOG_LEQ,
	NODE_LOG_GEQ,

	NODE_IF,
	NODE_WHILE,
	NODE_FOR,
	NODE_RANGED_FOR,

	NODE_RETURN,
	NODE_CONTINUE,
	NODE_BREAK,

	NODE_UNPACK,

	NODE_INVALID
};


extern std::map<NodeType, t_real (*)(t_real, t_real)> g_mapBinOps_d;
//extern Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op);
extern void safe_delete(Symbol *&pSym, const SymbolTable* pSymTab, const SymbolTable* pGlobalSymTab);


struct Node
{
protected:
	NodeType m_type;
	unsigned int m_iLine;

public:
	Node(NodeType ntype) : m_type(ntype), m_iLine(0) {}
	virtual ~Node() {}
	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const = 0;

	virtual Node* clone() const = 0;

	NodeType GetType() const { return m_type; }
	unsigned int GetLine() const { return m_iLine; }
	void SetLine(unsigned int iLine) { m_iLine = iLine; }

	// create a string containing the line number for error output
	t_string linenr(const t_string& strErr, const ParseInfo &info) const;

	static Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op);
};

struct NodeReturn : public Node
{
protected:
	Node *m_pExpr;

public:
	NodeReturn(Node *pExpr=0)
		: Node(NODE_RETURN), m_pExpr(pExpr)
	{}

	NodeReturn(void *pExpr)
		: NodeReturn((Node*)pExpr)
	{}

	virtual ~NodeReturn()
	{
		if(m_pExpr) delete m_pExpr;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetExpr() const { return m_pExpr; }
};

struct NodeBreak : public Node
{
	NodeBreak()
		: Node(NODE_BREAK)
	{}

	virtual ~NodeBreak()
	{}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;
};

struct NodeContinue : public Node
{
	NodeContinue()
		: Node(NODE_CONTINUE)
	{}

	virtual ~NodeContinue()
	{}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;
};

struct NodeIdent : public Node
{
protected:
	t_string m_strIdent;

public:
	NodeIdent(const t_string& strIdent)
		: Node(NODE_IDENT), m_strIdent(strIdent)
	{}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	const t_string& GetIdent() const { return m_strIdent; }
};

struct NodeCall : public Node
{
protected:
	Node *m_pIdent;
	Node *m_pArgs;

public:
	NodeCall(Node* pIdent, Node* pArgs);

	NodeCall(void* pIdent, void* pArgs)
		: NodeCall((Node*)pIdent, (Node*)pArgs)
	{}

	virtual ~NodeCall()
	{
		if(m_pIdent) delete m_pIdent;
		if(m_pArgs) delete m_pArgs;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetIdent() const { return m_pIdent; }
	Node* GetArgs() const { return m_pArgs; }


protected:
	std::vector<Node*> m_vecArgs;
};

struct NodeIf : public Node
{
protected:
	Node *m_pExpr;
	Node *m_pIf;
	Node *m_pElse;

public:
	NodeIf(Node* pExpr, Node* pIf, Node* pElse=0)
		: Node(NODE_IF), m_pExpr(pExpr), m_pIf(pIf), m_pElse(pElse)
	{}

	NodeIf(void* pExpr, void* pIf, void* pElse=0)
		: NodeIf((Node*)pExpr, (Node*)pIf, (Node*)pElse)
	{}

	virtual ~NodeIf()
	{
		if(m_pExpr) delete m_pExpr;
		if(m_pIf) delete m_pIf;
		if(m_pElse) delete m_pElse;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetExpr() const { return m_pExpr; }
	Node* GetIf() const { return m_pIf; }
	Node* GetElse() const { return m_pElse; }
};

struct NodeWhile : public Node
{
protected:
	Node *m_pExpr;
	Node *m_pStmt;

public:
	NodeWhile(Node* pExpr, Node* pStmt)
		: Node(NODE_WHILE), m_pExpr(pExpr), m_pStmt(pStmt)
	{}

	NodeWhile(void* pExpr, void* pStmt)
		: NodeWhile((Node*)pExpr, (Node*)pStmt)
	{}

	virtual ~NodeWhile()
	{
		if(m_pExpr) delete m_pExpr;
		if(m_pStmt) delete m_pStmt;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetExpr() const { return m_pExpr; }
	Node* GetStmt() const { return m_pStmt; }
};

struct NodeFor : public Node
{
protected:
	Node *m_pExprInit;
	Node *m_pExprCond;
	Node *m_pExprEnd;
	Node *m_pStmt;

public:
	NodeFor(Node* pExprInit, Node* pExprCond, Node* pExprEnd, Node* pStmt)
		: Node(NODE_FOR), m_pExprInit(pExprInit),
			m_pExprCond(pExprCond), m_pExprEnd(pExprEnd),
			m_pStmt(pStmt)
	{}

	NodeFor(void* pExprInit, void* pExprCond, void* pExprEnd, void* pStmt)
		: NodeFor((Node*)pExprInit, (Node*)pExprCond, (Node*)pExprEnd, (Node*)pStmt)
	{}

	virtual ~NodeFor()
	{
		if(m_pExprInit) delete m_pExprInit;
		if(m_pExprCond) delete m_pExprCond;
		if(m_pExprEnd) delete m_pExprEnd;
		if(m_pStmt) delete m_pStmt;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetExprInit() const { return m_pExprInit; }
	Node* GetExprCond() const { return m_pExprCond; }
	Node* GetExprEnd() const { return m_pExprEnd; }
	Node* GetStmt() const { return m_pStmt; }
};

struct NodeRangedFor : public Node
{
protected:
	Node *m_pIdent;
	Node *m_pExpr;
	Node *m_pStmt;

public:
	NodeRangedFor(Node* pIdent, Node* pExpr, Node* pStmt)
		: Node(NODE_RANGED_FOR),
		  m_pIdent(pIdent), m_pExpr(pExpr), m_pStmt(pStmt)
	{}

	NodeRangedFor(void* pIdent, void* pExpr, void* pStmt)
		: NodeRangedFor((Node*)pIdent, (Node*)pExpr, (Node*)pStmt)
	{}

	virtual ~NodeRangedFor()
	{
		if(m_pIdent) delete m_pIdent;
		if(m_pExpr) delete m_pExpr;
		if(m_pStmt) delete m_pStmt;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetIdent() const { return m_pIdent; }
	Node* GetExpr() const { return m_pExpr; }
	Node* GetStmt() const { return m_pStmt; }
};


struct NodeDouble : public Node
{
protected:
	t_real m_dVal;

public:
	NodeDouble(t_real dVal);
	virtual ~NodeDouble()
	{
		if(m_pSymbol) delete m_pSymbol;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	t_real GetVal() const { return m_dVal; }

protected:
	SymbolDouble *m_pSymbol;
};

struct NodeInt : public Node
{
protected:
	t_int m_iVal;

public:
	NodeInt(t_int iVal);
	virtual ~NodeInt()
	{
		if(m_pSymbol) delete m_pSymbol;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	t_int GetVal() const { return m_iVal; }

protected:
	SymbolInt *m_pSymbol;
};

struct NodeString : public Node
{
protected:
	t_string m_strVal;

public:
	NodeString(t_string strVal);
	virtual ~NodeString()
	{
		if(m_pSymbol) delete m_pSymbol;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	const t_string& GetVal() const { return m_strVal; }

protected:
	SymbolString *m_pSymbol;
};

struct NodeArray : public Node
{
protected:
	Node *m_pArr;

public:
	NodeArray(NodeArray* pArr);
	NodeArray(void *pArr)
		: NodeArray((NodeArray*)pArr)
	{}

	virtual ~NodeArray()
	{
		if(m_pArr) delete m_pArr;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetArr() const { return m_pArr; }
};

struct NodeMap : public Node
{
protected:
	Node *m_pMap;

public:
	NodeMap(NodeMap* pMap);
	NodeMap(void* pMap)
		: NodeMap((NodeMap*)pMap)
	{}

	virtual ~NodeMap()
	{
		if(m_pMap) delete m_pMap;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetMap() const { return m_pMap; }
};

struct NodePair : public Node
{
protected:
	Node *m_pFirst, *m_pSecond;

public:
	NodePair(Node *pFirst, Node *pSecond);
	NodePair(void *pFirst, void *pSecond)
		: NodePair((NodePair*)pFirst, (NodePair*)pSecond)
	{}

	virtual ~NodePair()
	{
		if(m_pFirst) delete m_pFirst;
		if(m_pSecond) delete m_pSecond;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetFirst() const { return m_pFirst; }
	Node* GetSecond() const { return m_pSecond; }
};

enum RangeType
{
	RANGE_BEGINEND,
	RANGE_FULL,
};

struct NodeRange : public Node
{
protected:
	RangeType m_rangetype;
	Node *m_pBegin;
	Node *m_pEnd;

public:
	NodeRange(RangeType rt);
	NodeRange(Node* pBegin, Node* pEnd);
	NodeRange(void* pBegin, void* pEnd)
		: NodeRange((Node*)pBegin, (Node*)pEnd)
	{}

	virtual ~NodeRange()
	{
		if(m_pBegin) delete m_pBegin;
		if(m_pEnd) delete m_pEnd;
	}

	void GetRangeIndices(ParseInfo &info, SymbolTable *pSym,
				t_int iMaxLen, t_int& iBeginIdx, t_int& iEndIdx);

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	RangeType GetRangeType() const { return m_rangetype; }
	Node* GetBegin() const { return m_pBegin; }
	Node* GetEnd() const { return m_pEnd; }
};

struct NodeArrayAccess : public Node
{
protected:
	Node *m_pIdent;
	Node *m_pExpr;

public:
	NodeArrayAccess(Node* pIdent, Node* pExpr);
	NodeArrayAccess(void* pIdent, void* pExpr)
		: NodeArrayAccess((Node*)pIdent, (Node*)pExpr)
	{}

	virtual ~NodeArrayAccess()
	{
		if(m_pIdent) delete m_pIdent;
		if(m_pExpr) delete m_pExpr;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetIdent() const { return m_pIdent; }
	Node* GetExpr() const { return m_pExpr; }

protected:
	std::vector<Node*> m_vecIndices;
};

struct NodeUnaryOp : public Node
{
protected:
	Node *m_pChild;

public:
	NodeUnaryOp(Node *pChild, NodeType ntype) 
		: Node(ntype), m_pChild(pChild)
	{}

	NodeUnaryOp(void* pChild, NodeType ntype)
		: NodeUnaryOp((Node*)pChild, ntype)
	{}

	virtual ~NodeUnaryOp()
	{
		if(m_pChild) delete m_pChild;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetChild() const { return m_pChild; }
};

struct NodeBinaryOp : public Node
{
protected:
	Node *m_pLeft, *m_pRight;
	bool m_bOwnsLeft = 1;
	bool m_bGlobal = 0;
	std::vector<Node*> m_vecNodesFlat;

public:
	NodeBinaryOp(Node* pLeft, Node* pRight, NodeType ntype);
	NodeBinaryOp(void* pLeft, void* pRight, NodeType ntype)
		: NodeBinaryOp((Node*)pLeft, (Node*)pRight, ntype)
	{}

	virtual ~NodeBinaryOp()
	{
		if(m_pLeft && m_bOwnsLeft) delete m_pLeft;
		if(m_pRight && ((void*)m_pRight)!=((void*)m_pLeft)) delete m_pRight;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Symbol* eval_assign(ParseInfo &info, SymbolTable *pSym=0, 
					Node* pLeft=0, Node *pRight=0, Symbol* pSymRightAlt=0, 
					const bool *pbGlob=0) const;
	virtual Symbol* eval_funcinit(ParseInfo &info, SymbolTable *pSym=0) const;

	virtual Symbol* eval_recursive(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Symbol* eval_sequential(ParseInfo &info, SymbolTable *pSym=0) const;

	virtual Node* clone() const;

	void FlattenNodes(NodeType ntype);
	std::vector<Node*> flatten(NodeType ntype=NODE_ARGS) const;

	Node* GetLeft() const { return m_pLeft; }
	Node* GetRight() const { return m_pRight; }
	const std::vector<Node*>& GetNodesFlat() const { return m_vecNodesFlat; }

	void SetGlobal(bool bGlob) { m_bGlobal = bGlob; }
	void SetOwnsLeft(bool bOwns) { m_bOwnsLeft = bOwns; }
};

struct NodeFunction : public Node
{
protected:
	Node *m_pIdent, *m_pArgs, *m_pStmts;

	// arguments the function takes
	std::vector<Node*> m_vecArgs;

	// script file this function resides in
	t_string m_strScrFile;

public:
	NodeFunction(Node* pLeft, Node* pMiddle, Node* pRight);

	NodeFunction(void* pLeft, void *pMiddle, void* pRight)
		: NodeFunction((Node*)pLeft, (Node*)pMiddle, (Node*)pRight)
	{}

	virtual ~NodeFunction()
	{
		if(m_pIdent) delete m_pIdent;
		if(m_pArgs) delete m_pArgs;
		if(m_pStmts) delete m_pStmts;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	Node* GetIdent() const { return m_pIdent; }
	Node* GetArgs() const { return m_pArgs; }
	Node* GetStmts() const { return m_pStmts; }

	const std::vector<Node*>& GetArgVec() const { return m_vecArgs; }

	const t_string& GetScrFile() const { return m_strScrFile; }
	void SetScrFile(const t_string& str) { m_strScrFile = str; }

	const t_string& GetName() const;
	std::vector<t_string> GetParamNames() const;
};

#endif
