/*
 * Script interpreter
 * @author tweber
 * @date 2013-2014
 */

#ifndef __MIEZE_NODE__
#define __MIEZE_NODE__

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
	typedef std::map<std::string, Node*> t_mods;
	t_mods *pmapModules;

	// function to execute, e.g. "main"
	std::string strExecFkt;
	const std::vector<Symbol*>* pvecExecArg;
	std::string strInitScrFile;

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

	NodeFunction* GetFunction(const std::string& strName);

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
	NODE_RANGED_FOR,
	
	NODE_RETURN,
	NODE_CONTINUE,
	NODE_BREAK,
	
	
	NODE_INVALID
};


extern std::map<NodeType, double (*)(double,double)> g_mapBinOps_d;
extern Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op);
extern void safe_delete(Symbol *pSym, const SymbolTable* pSymTab);


struct Node
{
	NodeType m_type;
	unsigned int m_iLine;

	Node(NodeType ntype) : m_type(ntype), m_iLine(0) {}
	virtual ~Node() {}
	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const = 0;

	virtual Node* clone() const = 0;

	// create a string containing the line number for error output
	std::string linenr(const std::string& strErr, const ParseInfo &info) const;

	static Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op);
};

struct NodeReturn : public Node
{
	Node *m_pExpr;
	
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
	std::string m_strIdent;

	NodeIdent(const std::string& strIdent)
		: Node(NODE_IDENT), m_strIdent(strIdent)
	{}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;
};

struct NodeCall : public Node
{
	Node *m_pIdent;
	Node *m_pArgs;

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

	
protected:
	std::vector<Node*> m_vecArgs;
};

struct NodeIf : public Node
{
	Node *m_pExpr;
	Node *m_pIf;
	Node *m_pElse;

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
};

struct NodeWhile : public Node
{
	Node *m_pExpr;
	Node *m_pStmt;

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
};

struct NodeRangedFor : public Node
{
	Node *m_pIdent;
	Node *m_pExpr;
	Node *m_pStmt;

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
};


struct NodeDouble : public Node
{
	double m_dVal;
	
	NodeDouble(double dVal);
	virtual ~NodeDouble()
	{
		if(m_pSymbol) delete m_pSymbol;
	}
	
	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;
	
protected:
	SymbolDouble *m_pSymbol;
};

struct NodeInt : public Node
{
	double m_iVal;
	
	NodeInt(int iVal);
	virtual ~NodeInt()
	{
		if(m_pSymbol) delete m_pSymbol;
	}
	
	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

protected:
	SymbolInt *m_pSymbol;
};

struct NodeString : public Node
{
	std::string m_strVal;

	NodeString(std::string strVal);
	virtual ~NodeString()
	{
		if(m_pSymbol) delete m_pSymbol;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

protected:
	SymbolString *m_pSymbol;
};

struct NodeArray : public Node
{
	Node *m_pArr;
	
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
};

struct NodeMap : public Node
{
	Node *m_pMap;

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
};

struct NodePair : public Node
{
	Node *m_pFirst, *m_pSecond;

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
};

struct NodeArrayAccess : public Node
{
	Node *m_pIdent;
	Node *m_pExpr;
	
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

protected:
	std::vector<Node*> m_vecIndices;
};

struct NodeUnaryOp : public Node
{
	Node *m_pChild;

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
};

struct NodeBinaryOp : public Node
{
	Node *m_pLeft, *m_pRight;
	bool m_bGlobal;

	NodeBinaryOp(Node* pLeft, Node* pRight, NodeType ntype)
		: Node(ntype), m_pLeft(pLeft), m_pRight(pRight), m_bGlobal(0)
	{}

	NodeBinaryOp(void* pLeft, void* pRight, NodeType ntype)
		: NodeBinaryOp((Node*)pLeft, (Node*)pRight, ntype)
	{}

	virtual ~NodeBinaryOp()
	{
		if(m_pLeft) delete m_pLeft;
		if(m_pRight && ((void*)m_pRight)!=((void*)m_pLeft)) delete m_pRight;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Symbol* eval_assign(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Symbol* eval_funcinit(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Symbol* eval_sequential(ParseInfo &info, SymbolTable *pSym=0) const;
	virtual Node* clone() const;

	std::vector<Node*> flatten(NodeType ntype=NODE_ARGS) const;
};

struct NodeFunction : public Node
{
	Node *m_pIdent, *m_pArgs, *m_pStmts;

	// symbols for the arguments the function takes
	const std::vector<Symbol*>* m_pVecArgSyms;

	// script file this function resides in
	std::string m_strScrFile;

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

	std::string GetName() const;

	// arguments to function
	void SetArgSyms(const std::vector<Symbol*>* pvecArgSyms)
	{ this->m_pVecArgSyms = pvecArgSyms; }
	
	
//protected:
	// arguments the function takes
	std::vector<Node*> m_vecArgs;
};

#endif
