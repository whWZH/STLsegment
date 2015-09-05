//#include "tin2.h"
#include "stdafx.h"
#include "T_line.h"

T_line::T_line()
{
	p0=-1;
	p1=-1;
	useCount=0;

}

T_line::~T_line()
{

}

bool operator==(T_line line0,T_line line1)
{
	if(line0.p0==line1.p0&&line0.p1==line1.p1)
		return true;
	if(line0.p0==line1.p1&&line0.p1==line1.p0)
		return true;
	return false;
}
