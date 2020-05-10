// residue.h
//
// Last Modified: 11, Jul 2007
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

#ifndef __RESIDUE_H__
#define __RESIDUE_H__

#include <vector>
#include <string>
#include <functional>
#ifndef CLIMITS
#include <limits>
#else
#include <climits>
#endif

namespace prrn
{
    typedef char   RESIDUE;
    typedef unsigned char RESCODE;

    const RESCODE nonPenaltyGap = 0;
    const RESCODE PenaltyGap = 1;

    const RESCODE nAminoAcid = 20;
    const RESCODE nNucleicAcid = 4;

    const RESCODE nNonStdResidue = 3;
    const RESCODE nNonStdBase = 11;

    const RESCODE nAACode = 25;
    const RESCODE nNACode = 17;

#ifndef CLIMITS
    const RESCODE limitRescodeMax = std::numeric_limits<RESCODE>::max();
#else
    const RESCODE limitRescodeMax = UCHAR_MAX;
#endif

    enum rtype {AA, DNA, RNA};

    const RESIDUE AACode[nAACode] = {
	'.', '-', 'X', 'A', 'R',
	'N', 'D', 'C', 'Q', 'E',
	'G', 'H', 'I', 'L', 'K',
	'M', 'F', 'P', 'S', 'T',
	'W', 'Y', 'V', 'B', 'Z'};

    // A = 1 (+1)
    // C = 2 (+1)
    // G = 4 (+1)
    // T, U = 8 (+1)
    const RESIDUE DNACode[nNACode] = {
	'.', '-', 'N', 'G', 'A', 'T',
	'C', 'R', 'Y', 'M', 'K', 'S',
	'W', 'H', 'B', 'V', 'D'};

    const RESIDUE RNACode[nNACode] = {
	'.', '-', 'N', 'G', 'A', 'U',
	'C', 'R', 'Y', 'M', 'K', 'S',
	'W', 'H', 'B', 'V', 'D'};

    inline bool isResidue(RESCODE r)
    {
	return (r > PenaltyGap);
    }

    const RESIDUE* getResCode(rtype rt);

    struct getAACode : public std::unary_function<RESCODE, RESIDUE>
    {
	result_type operator()(argument_type i)
	{
	    return AACode[i];
	}
    };

    struct getDNACode : public std::unary_function<RESCODE, RESIDUE>
    {
	result_type operator()(argument_type i)
	{
	    return DNACode[i];
	}
    };

    struct getRNACode : public std::unary_function<RESCODE, RESIDUE>
    {
	result_type operator()(argument_type i)
	{
	    return RNACode[i];
	}
    };

    RESCODE getAANum(RESIDUE c);
    RESCODE getNANum(RESIDUE c);

    void convert2rescode(std::string::const_iterator,
	    std::string::const_iterator, std::vector<RESCODE>&);
    void convert2string(std::vector<RESCODE>::const_iterator,
	    std::vector<RESCODE>::const_iterator, std::string::iterator);
}

#endif
