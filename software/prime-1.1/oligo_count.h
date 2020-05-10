// oligo_count.h
//
// Last Modified: 9, May 2007
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

#ifndef __OLIGO_COUNT_H__
#define __OLIGO_COUNT_H__

#include "residue.h"
#include <vector>

class Sequence;

// Dayhoff(6): AGPST(X), C, DENQ(BZ), FWY, HKR, ILMV
const char dayhoff[prrn::nAACode] = {
    6, 6, 0, 0, 4,
    2, 2, 1, 2, 2,
    0, 4, 5, 5, 4,
    5, 3, 0, 0, 0,
    3, 3, 5, 2, 2};
const size_t nrc_dayhoff = 6;

// SE-B(10): AST(X), C, DN(B), EQ(Z), FY, G, HW, ILMV, KR, P
const char se_b[prrn::nAACode] = {
    10, 10, 0, 0, 8,
    2, 2, 1, 3, 3,
    5, 6, 7, 7, 8,
    7, 4, 9, 0, 0,
    6, 4, 7, 2, 3};
const size_t nrc_se_b = 10;

// WG(8): WFY, MLIV, C, G, P, ATS(X), NDEQ(BZ), HRK
// Note that this classification has been modified due to non-standard residue codes
const char wg[prrn::nAACode] = {
    8, 8, 5, 5, 7,
    6, 6, 2, 6, 6,
    3, 7, 1, 1, 7,
    1, 0, 4, 5, 5,
    0, 0, 1, 6, 6};
const size_t nrc_wg = 8;

/*
const RESIDUE AACode[nAACode] = {
    '.', '-', 'X', 'A', 'R',
    'N', 'D', 'C', 'Q', 'E',
    'G', 'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S', 'T',
    'W', 'Y', 'V', 'B', 'Z'};
    */

class Oligomer
{
    public:
	Oligomer();
	Oligomer(const Sequence&);
	~Oligomer(){}

	static void setResidueClass(const char*, size_t);
	static void setOligomerLength(size_t);
	static double getApproxSeqID(const Oligomer&, const Oligomer&);

    private:
	size_t nsgmt_;
	std::vector<size_t> oligo_;

	static size_t kmer_;

	static const char* pRC_;
	static size_t nrc_;
};

#endif
