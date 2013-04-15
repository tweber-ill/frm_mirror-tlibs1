/*
 * Simple Script
 * @author tweber
 */

#ifndef __MIEZE_NODE__
#define __MIEZE_NODE__

#include <math.h>
#include <string>
#include <iostream>
#include <map>
#include "symbol.h"


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
	NODE_IDENT
};


extern std::map<NodeType, double (*)(double,double)> g_mapBinOps_d;
extern Symbol* Op(const Symbol *pSymLeft, const Symbol *pSymRight, NodeType op);
extern void safe_delete(Symbol *pSym, const SymbolTable* pSymTab);


struct Node
{
	NodeType m_type;

	Node(NodeType ntype) : m_type(ntype)
	{}

	virtual ~Node()
	{}
	virtual Symbol* eval(SymbolTable *pSym) const = 0;
};

struct NodeIdent : public Node
{
	std::string m_strIdent;

	NodeIdent(const std::string& strIdent)
		: Node(NODE_IDENT), m_strIdent(strIdent)
	{}

	virtual Symbol* eval(SymbolTable *pSym) const
	{
		return pSym->GetSymbol(m_strIdent);
	}
};

struct NodeCall : public Node
{
        Node *m_pIdent;

	NodeCall(Node* pIdent)
		: Node(NODE_CALL), m_pIdent(pIdent)
	{}

	NodeCall(void* pIdent)
		: NodeCall((Node*)pIdent)
	{}

	virtual ~NodeCall()
	{
		if(m_pIdent) delete m_pIdent;
	}
	
	virtual Symbol* eval(SymbolTable *pSym) const
	{
		return 0;
	}
};

struct NodeDouble : public Node
{
	double m_dVal;
	
	NodeDouble(double dVal)
		: Node(NODE_DOUBLE), m_dVal(dVal)
	{}
	
	virtual Symbol* eval(SymbolTable *pSym) const
	{
		SymbolDouble *pSymbol = new SymbolDouble;
		pSymbol->dVal = m_dVal;
		pSymbol->strName = "<const>";
		return pSymbol;
	}
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

	virtual Symbol* eval(SymbolTable *pSym) const
	{
		switch(m_type)
		{
			case NODE_UMINUS:
			{
				Symbol *pSymbolEval = m_pChild->eval(pSym);
				Symbol *pSymbol = pSymbolEval->clone();
				safe_delete(pSymbolEval, pSym);
				
				if(pSymbol->GetType() == SYMBOL_DOUBLE)
					((SymbolDouble*)pSymbol)->dVal = -((SymbolDouble*)pSymbol)->dVal;
				
				return pSymbol;
			}
		}
	}
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

        virtual Symbol* eval(SymbolTable *pSym) const
        {
        	switch(m_type)
        	{
        		case NODE_STMTS:
        			if(m_pLeft)
					{
						Symbol *pSymbol = m_pLeft->eval(pSym);
						safe_delete(pSymbol, pSym);
					}
        			if(m_pRight)
					{
						Symbol *pSymbol = m_pRight->eval(pSym);
						safe_delete(pSymbol, pSym);
					}
        			return 0;
        		case NODE_ASSIGN:
        			Symbol *pSymbol = m_pRight->eval(pSym);
        			const std::string& strIdent = ((NodeIdent*)m_pLeft)->m_strIdent;

					pSym->InsertSymbol(strIdent, pSymbol);
        			return pSymbol;
        	};

			Symbol *pSymbolLeft = m_pLeft->eval(pSym);
			Symbol *pSymbolRight = m_pRight->eval(pSym);
        	Symbol *pSymbol = Op(pSymbolLeft, pSymbolRight, m_type);
			safe_delete(pSymbolLeft, pSym);
			safe_delete(pSymbolRight, pSym);
			return pSymbol;
        }
};

#endif
