// residue.cpp
//
// Last Modified: 25, Feb 2007
//
// Copyright (c) 2004-2007 Shinsuke Yamada
//
// This file is part of PRIME.
//
// PRIME is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PRIME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PRIME; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "residue.h"
#include "params.h"
#include <iostream>
#include <algorithm>
#include <cctype>

using namespace std;

namespace prrn
{
    const RESIDUE* getResCode(rtype rt)
    {
	switch(rt)
	{
	    case DNA:
		return DNACode;
	    case RNA:
		return RNACode;
	    case AA:
	    default:
		return AACode;
	}
    }

    RESCODE getAANum(RESIDUE c)
    {
	switch(toupper(c))
	{
	    case '.':	// for not penalized gap
		return 0;
	    case '-':
	    case '*':
		return 1;
	    case 'X':
		return 2;
	    case 'A':
		return 3;
	    case 'R':
		return 4;
	    case 'N':
		return 5;
	    case 'D':
		return 6;
	    case 'C':
		return 7;
	    case 'Q':
		return 8;
	    case 'E':
		return 9;
	    case 'G':
		return 10;
	    case 'H':
		return 11;
	    case 'I':
		return 12;
	    case 'L':
		return 13;
	    case 'K':
		return 14;
	    case 'M':
		return 15;
	    case 'F':
		return 16;
	    case 'P':
		return 17;
	    case 'S':
		return 18;
	    case 'T':
		return 19;
	    case 'W':
		return 20;
	    case 'Y':
		return 21;
	    case 'V':
		return 22;
	    case 'B':
		return 23;
	    case 'Z':
		return 24;
	    default:
		std::cerr << "irregular residue `" << c << "' has been included; it is substituted by `X'\n";
		return 2;
	}
	return nAACode;
    }

    RESCODE getNANum(RESIDUE c)
    {
	switch(toupper(c))
	{
	    case '.':	// for not penalized gap
		return 0;
	    case '-':
	    case '*':
		return 1;
	    case 'N':
		return 2;
	    case 'G':
		return 3;
	    case 'A':
		return 4;
	    case 'T':
	    case 'U':
		return 5;
	    case 'C':
		return 6;
	    case 'R':
		return 7;
	    case 'Y':
		return 8;
	    case 'M':
		return 9;
	    case 'K':
		return 10;
	    case 'S':
		return 11;
	    case 'W':
		return 12;
	    case 'H':
		return 13;
	    case 'B':
		return 14;
	    case 'V':
		return 15;
	    case 'D':
		return 16;
	    default:
		std::cerr << "irregular residue `" << c << "' has been included; it is substituted by `N'\n";
		return 2;
	}
	return nNACode;
    }

    void convert2rescode(string::const_iterator begin,
	    string::const_iterator end,
	    vector<RESCODE>& target)
    {
	transform(begin, end, back_inserter(target),
		(nResType == nAACode) ? getAANum : getNANum);
    }

    void convert2string(vector<RESCODE>::const_iterator begin,
	    vector<RESCODE>::const_iterator end,
	    string::iterator target)
    {
	switch(res_type)
	{
	    case DNA:
		transform(begin, end, target, getDNACode());
		break;
	    case RNA:
		transform(begin, end, target, getRNACode());
		break;
	    case AA:
	    default:
		transform(begin, end, target, getAACode());
	}
    }
}
