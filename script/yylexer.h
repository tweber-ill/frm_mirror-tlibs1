/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#ifndef __MIEZE_YYLEXER__
#define __MIEZE_YYLEXER__

#include "tokens.h"
#include "parseobj.h"

#ifdef __cplusplus
extern "C" {
#endif


int yylex(void*, void* _pParseObj=0);
void yyerror(void* pParseObj, const char*);
//void yyerror(const char*);


#ifdef __cplusplus
}
#endif


#endif
