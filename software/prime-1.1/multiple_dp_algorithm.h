// multiple_dp_algorithm.h
//
// Last Modified: 31, May 2007
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

#ifndef __MULTIPLE_DP_ALGORITHM_H__
#define __MULTIPLE_DP_ALGORITHM_H__

#include "substitution_matrix.h"
#include <vector>

class Sequence;
class Alignment;
class PTHolder;
class PhylogeneticTree;
class PairwiseDPAlgorithm;

class MultipleDPAlgorithm
{
    public:
	virtual ~MultipleDPAlgorithm(){}
	virtual void showAlignmentMode() const = 0;
	virtual Alignment getAlignment(const Alignment&, const Alignment&,
		const PTHolder* = 0) const = 0;
	virtual prrn::SCORE getPSAscore(const Alignment&, size_t, size_t) const = 0;
	virtual prrn::SCORE getSPScore(const Alignment&, const PTHolder* = 0) const = 0;

#ifdef MSA_TEST
	double getDPScore() const
	{
	    return score_;
	}
#endif

    protected:
	MultipleDPAlgorithm()
	    : sm_(SubstitutionMatrix::getInstance())
	    , go_(sm_.getOpenPenalty())
	    , ge_(sm_.getExtensionPenalty())
	    {}

	const SubstitutionMatrix& sm_;
	const prrn::SCORE go_;
	const prrn::SCORE ge_;

#ifdef MSA_TEST
	mutable double score_;
#endif
};

void getProgressiveAlignment(const std::vector<Sequence>&,
	const PairwiseDPAlgorithm&, const MultipleDPAlgorithm&,
	const PTHolder&, Alignment&, bool = false);

void getBranchOrder(int, std::vector<int>&);

void refineAlignment(const MultipleDPAlgorithm&, size_t,
	const PTHolder&, std::vector<int>&, Alignment&);

#endif
