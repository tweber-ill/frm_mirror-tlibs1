/*
 * Script interpreter
 * @author tweber
 */

#ifndef __MIEZE_NODE__
#define __MIEZE_NODE__

#include <math.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "symbol.h"
#include "handles.h"


struct NodeFunction;
struct ParseInfo
{
	SymbolTable* pGlobalSyms;
	std::vector<NodeFunction*> vecFuncs;
	HandleManager handles;
	
	bool bWantReturn = 0;
	
	ParseInfo() : bWantReturn(0)
	{}


	NodeFunction* GetFunction(const std::string& strName);
};


enum NodeType
{
	NODE_NOP,

	NODE_PLUS,
	NODE_MINUS,
	NODE_DIV,
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
	
	
	NODE_INVALID
};


extern std::map<NodeType, double (*)(double,double)> g_mapBinOps_d;
extern Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op);
extern void safe_delete(Symbol *pSym, const SymbolTable* pSymTab);


struct NodeFunction;

struct Node
{
	NodeType m_type;

	Node(NodeType ntype) : m_type(ntype) {}
	virtual ~Node() {}
	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const = 0;
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
};

struct NodeIdent : public Node
{
	std::string m_strIdent;

	NodeIdent(const std::string& strIdent)
		: Node(NODE_IDENT), m_strIdent(strIdent)
	{}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;
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

protected:
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
};

struct NodeBinaryOp : public Node
{
	Node *m_pLeft, *m_pRight;

	NodeBinaryOp(Node* pLeft, Node* pRight, NodeType ntype)
		: Node(ntype), m_pLeft(pLeft), m_pRight(pRight)
	{}

	NodeBinaryOp(void* pLeft, void* pRight, NodeType ntype)
		: NodeBinaryOp((Node*)pLeft, (Node*)pRight, ntype)
	{}

	virtual ~NodeBinaryOp()
	{
		if(m_pLeft) delete m_pLeft;
		if(m_pRight) delete m_pRight;
	}

	virtual Symbol* eval(ParseInfo &info, SymbolTable *pSym=0) const;

	std::vector<Node*> flatten(NodeType ntype=NODE_ARGS) const;
};

struct NodeFunction : public Node
{
	Node *m_pIdent, *m_pArgs, *m_pStmts;
	const std::vector<Symbol*>* m_pVecArgSyms;

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

	std::string GetName() const;

	void SetArgSyms(const std::vector<Symbol*>* pvecArgSyms)
	{ this-> m_pVecArgSyms = pvecArgSyms; }
	
	
protected:
	std::vector<Node*> m_vecArgs;
};

#endif
