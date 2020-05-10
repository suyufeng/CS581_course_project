// profile_dp_long.h
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

#ifndef __PROFILE_DP_LONG_H__
#define __PROFILE_DP_LONG_H__

#include "multiple_dp_algorithm.h"
#include "dynamic_gap_state.h"
#include "count_sequence.h"

class Alignment;
class AverageProfile;
class GapProfile;

class ProfileDPLong : public MultipleDPAlgorithm
{
    public:
	ProfileDPLong();
	~ProfileDPLong(){}
	void showAlignmentMode() const;
    private:
	typedef std::pair<int, int> DGI;

	size_t getStartPosition(size_t, size_t, size_t&,
		std::vector<DGI>::const_reverse_iterator&,
		const std::vector<DGI>::const_reverse_iterator&) const;
	prrn::SCORE getGapPenalty(const GapProfile&, const DynamicGapStates&,
		const GapProfile&, const DynamicGapStates&,
		std::vector<DGI>&, const std::vector<std::pair<size_t, double> >&,
		size_t) const;

	prrn::SCORE calc_diag(const AverageProfile&, size_t, const CountSequence&,
		const AverageProfile&, size_t, const CountSequence&, const cell4lng&) const;
	prrn::SCORE calc_vert(const AverageProfile&, size_t,
		const GapProfile&, const cell4lng&, const CountSequence&) const;
	prrn::SCORE calc_hori(const GapProfile&,
		const AverageProfile&, size_t, const cell4lng&, const CountSequence&) const;

    public:
	Alignment getAlignment(const Alignment&, const Alignment&,
		const PTHolder* = 0) const;
	prrn::SCORE getPSAscore(const Alignment&, size_t, size_t) const;
	prrn::SCORE getSPScore(const Alignment&, const PTHolder* = 0) const;
};

#endif
