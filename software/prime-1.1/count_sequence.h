// count_sequence.h
//
// Last Modified: 31, Jan 2007
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

#ifndef __COUNT_SEQUENCE_H__
#define __COUNT_SEQUENCE_H__

#include "prrn.h"
#include <vector>
#include <cstddef>

class Alignment;

class CountSequence
{
    public:
	CountSequence(const Alignment&, const std::vector<double>&);
	~CountSequence(){}

	const std::vector<std::pair<size_t, double> >& getRangeProfile(size_t pos) const
	{
	    return rangeprof_[pos];
	}

	prrn::SCORE getExtensionPenalty(size_t, size_t) const;

    private:
	void makeRangeList(const Alignment&, const std::vector<double>&);
	CountSequence();	// undefined
	CountSequence(const CountSequence&);	// undefined
	CountSequence& operator=(const CountSequence&);	// undefined

	const size_t gtr_;

	const size_t nseq_;
	const size_t alen_;

	std::vector<std::vector<std::pair<size_t, double> > > rangeprof_;
};

#endif
