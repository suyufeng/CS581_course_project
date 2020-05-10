// pairwise_dp_algorithm.h
//
// Last Modified: 14, Jul 2007
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

#ifndef __PAIRWISE_DP_ALGORITHM_H__
#define __PAIRWISE_DP_ALGORITHM_H__

#include "prrn.h"
#include "substitution_matrix.h"
#include <vector>

class Alignment;
class Sequence;

class PairwiseDPAlgorithm
{
    public:
	virtual ~PairwiseDPAlgorithm(){}

	virtual void showAlignmentMode() const = 0;
	virtual Alignment getAlignment(const Sequence&, const Sequence&) const = 0;
    protected:
	PairwiseDPAlgorithm()
	    : sub_(SubstitutionMatrix::getInstance())
	    , go_(sub_.getOpenPenalty())
	    {}

	const SubstitutionMatrix& sub_;
	const prrn::SCORE go_;
};

class PairwiseDPGlobal : public PairwiseDPAlgorithm
{
    public:
	PairwiseDPGlobal();
	~PairwiseDPGlobal(){}

	virtual void showAlignmentMode() const;
	virtual Alignment getAlignment(const Sequence&, const Sequence&) const;

    private:
	const prrn::SCORE ge_;
};

class PairwiseDPLoclGlobal : public PairwiseDPAlgorithm
{
    public:
	PairwiseDPLoclGlobal();
	~PairwiseDPLoclGlobal(){}

	virtual void showAlignmentMode() const;
	virtual Alignment getAlignment(const Sequence&, const Sequence&) const;

    private:
	const prrn::SCORE ge_;
};

class PairwiseDPLong : public PairwiseDPAlgorithm
{
    public:
	PairwiseDPLong();
	~PairwiseDPLong(){}

	virtual void showAlignmentMode() const;
	virtual Alignment getAlignment(const Sequence&, const Sequence&) const;

    private:
	const size_t L;	// the number of picewise linear functions
	const std::vector<prrn::SCORE>& go_;	// gap open penalties
	const std::vector<prrn::SCORE>& ge_;	// gap extension penalties
};

class PairwiseDPLoclLong : public PairwiseDPAlgorithm
{
    public:
	PairwiseDPLoclLong();
	~PairwiseDPLoclLong(){}
	virtual void showAlignmentMode() const;
	virtual Alignment getAlignment(const Sequence& , const Sequence&) const;

    private:
	const size_t L;	// the number of picewise linear functions
	const std::vector<prrn::SCORE>& go_;	// gap open penalties
	const std::vector<prrn::SCORE>& ge_;	// gap extension penalties
};

#endif
