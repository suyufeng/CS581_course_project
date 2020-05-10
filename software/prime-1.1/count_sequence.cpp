// count_sequence.cpp
//
// Last Modified: 12, Mar 2007
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

#include "count_sequence.h"
#include "alignment.h"
#include "substitution_matrix.h"

using namespace std;
using namespace prrn;

CountSequence::CountSequence(const Alignment& ra, const vector<double>& weight)
    : gtr_(SubstitutionMatrix::getInstance().getGapThreshold(0))
    , nseq_(ra.getNumSequences())
    , alen_(ra.getLength())
    , rangeprof_(alen_)
{
    if(gtr_ != limitSizeMax)
	makeRangeList(ra, weight);
    else
    {
	const SCORE gef = SubstitutionMatrix::getInstance().getExtensionPenalty(0);
	double tfreq = 0.0;
	for(size_t i = 0; i < alen_; ++i)
	{
	    tfreq = 0.0;
	    for(size_t j = 0; j < nseq_; ++j)
		if(PenaltyGap < ra.getResidueCode(j, i))
		    tfreq += weight[j];
	    rangeprof_[i].push_back(make_pair(0, gef * tfreq));
	}
    }
}

void CountSequence::makeRangeList(const Alignment& ra, const vector<double>& weight)
{
    const SCORE gef = SubstitutionMatrix::getInstance().getExtensionPenalty(0);
    const SCORE ges = SubstitutionMatrix::getInstance().getExtensionPenalty(1);

    vector<size_t> count(nseq_, 0);
    for(size_t i = 0; i < alen_; ++i)
    {
	double tfreq = 0.0;
	for(size_t l = 0; l < nseq_; ++l)
	    if(PenaltyGap < ra.getResidueCode(l, i))
	    {
		++count[l];
		tfreq += weight[l];
	    }
	double ns = tfreq;

	for(int j = i-1; j >= 0; --j)
	{
	    double tns = ns;
	    for(size_t k = 0; k < nseq_; ++k)
	    {
		if(!count[k])
		    continue;

		if(PenaltyGap < ra.getResidueCode(k, j)
			&& gtr_ < ++count[k])
		{
		    ns -= weight[k];
		    count[k] = 0;
		}
	    }

	    if(ns < tns)
	    {
		rangeprof_[i].push_back(make_pair(j+1,
			    gef*tns + ges*(tfreq-tns)));
		if(isApproxEqual(ns, 0.0))
		    break;
	    }
	}
	rangeprof_[i].push_back(make_pair(0,
		    gef*ns + ges*(tfreq-ns)));
	fill(count.begin(), count.end(), 0);
    }
}

SCORE CountSequence::getExtensionPenalty(size_t begin, size_t end) const
{
#ifdef DEBUG
    assert(begin <= end);
    assert(end < alen_);
#endif
    vector<pair<size_t, double> >::const_iterator
	i = rangeprof_[end].begin();
    if(end - begin < gtr_ || i->first <= begin)
	return i->second;

    const vector<pair<size_t, double> >::const_iterator
	ie = rangeprof_[end].end();
    for(++i; i != ie; ++i)
	if(i->first <= begin)
	    return i->second;
    return 0.0; // not executed
}
