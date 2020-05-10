// profile_dp_global.h
//
// Last Modified: 1, Jun 2007
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

#ifndef __PROFILE_DP_GLOBAL_H__
#define __PROFILE_DP_GLOBAL_H__

#include "multiple_dp_algorithm.h"
#include "dynamic_gap_state.h"

class Alignment;
class PhylogeneticTree;
class AverageProfile;
class GapProfile;

class ProfileDPGlobal : public MultipleDPAlgorithm
{
    public:
	ProfileDPGlobal();
	~ProfileDPGlobal(){}
	void showAlignmentMode() const;

    private:
	double getNumGaps(const GapProfile&, const DynamicGapStates&,
		const GapProfile&, const DynamicGapStates&) const;

	prrn::SCORE calc_diag(const AverageProfile&, size_t,
		const AverageProfile&, size_t, const cell&) const;
	prrn::SCORE calc_vert(const GapProfile&, const GapProfile&, const cell&) const;
	prrn::SCORE calc_hori(const GapProfile&, const GapProfile&, const cell&) const;

#ifndef MSA_TEST
	prrn::SCORE getWSPScore(size_t, AverageProfile& rap) const;
#endif

    public:
	Alignment getAlignment(const Alignment&, const Alignment&,
		const PTHolder* = 0) const;

	prrn::SCORE getPSAscore(const Alignment&, size_t, size_t) const;
	prrn::SCORE getInterProfileScore(const AverageProfile&, const AverageProfile&) const;

	prrn::SCORE getSPScore(const Alignment&, const PTHolder* = 0) const;
};

#endif
