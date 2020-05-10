// res_conserv.cpp
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

#include "params.h"
#include "alignment.h"
#include "sequence.h"
#include "res_conserv.h"
#ifdef DEBUG
#include "util.h"
#endif
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cmath>
#ifdef DEBUG
#include <cassert>
#endif

using namespace std;
using namespace prrn;

void calcHenikoffWeight(const Alignment& a, vector<double>& weight)
{
    const size_t nseq = a.getNumSequences();
#ifdef DEBUG
    assert(nseq == weight.size());
#endif
    fill(weight.begin(), weight.end(), 0.0);
    const size_t len = a.getLength();
    map<RESCODE, size_t> rt;
    for(size_t i = 0; i < len; ++i)
    {
	for(size_t j = 0; j < nseq; ++j)
	    ++rt[a.getResidueCode(j, i)];

	const double rtconst = 1.0 / (double)rt.size();
	for(size_t j = 0; j < nseq; ++j)
	    weight[j] += rtconst / (double)rt[a.getResidueCode(j, i)];
	rt.clear();
    }
#ifdef DEBUG
    assert(isApproxEqual(len, accumulate(weight.begin(), weight.end(), 0.0)));
#endif
    transform(weight.begin(), weight.end(), weight.begin(),
	    bind2nd(multiplies<double>(), 1.0/accumulate(weight.begin(), weight.end(), 0.0)));
}

namespace
{
    vector<double> background;
}

void calcBackground(const vector<Sequence>& s, size_t nrt)
{
    background.resize(nrt, 0.0);
    for(vector<Sequence>::const_iterator i = s.begin(), ie = s.end(); i != ie; ++i)
	for(Sequence::const_iterator j = i->begin(), je = i->end(); j != je; ++j)
	    ++background[*j];
    transform(background.begin(), background.end(), background.begin(),
	    bind2nd(multiplies<double>(), 1.0/accumulate(background.begin(), background.end(), 0.0)));
}

void calcBackground(const Alignment& a, size_t nrt)
{
    background.resize(nrt, 0.0);
    const size_t len = a.getLength();
    const size_t nseq = a.getNumSequences();
    for(size_t i = 0; i < nseq; ++i)
	for(size_t j = 0; j < len; ++j)
	    ++background[a.getResidueCode(i, j)];
    // ignore null characters
    background[0] = background[1] = 0.0;
    transform(background.begin(), background.end(), background.begin(),
	    bind2nd(multiplies<double>(), 1.0/accumulate(background.begin(), background.end(), 0.0)));
}

void calcREconservation(const Alignment& a, const vector<double>& weight, vector<double>& rescns)
{
    const size_t len = a.getLength();
#ifdef DEBUG
    assert(len == rescns.size());
#endif
    const size_t nseq = a.getNumSequences();
    fill(rescns.begin(), rescns.end(), 0.0);
    map<RESCODE, double> rt;
    for(size_t i = 0; i < len; ++i)
    {
	for(size_t j = 0; j < nseq; ++j)
	    rt[a.getResidueCode(j, i)] += weight[j];
	if(rt.find(PenaltyGap) == rt.end() && rt.find(nonPenaltyGap) == rt.end())
	{
	    for(map<RESCODE, double>::const_iterator k = rt.begin(), ke = rt.end(); k != ke; ++k)
		rescns[i] += k->second * log(k->second / background[k->first]);
	}
	rt.clear();
    }
    transform(rescns.begin(), rescns.end(), rescns.begin(),
	    bind2nd(multiplies<double>(), 1.0/log(2.0)));
}

void showResidueConservation(const Alignment& a, const vector<double>& rescns, ostream& os)
{
    const size_t len = a.getLength();
    const size_t nseq = a.getNumSequences();
#ifdef DEBUG
    assert(len == rescns.size());
#endif
    const RESIDUE* rescode = getResCode(res_type);
    for(size_t i = 0; i < len; ++i)
    {
	//os << setw(3) << i << setw(3) << " ";
	for(size_t j = 0; j < nseq; ++j)
	    os << (char)rescode[a.getResidueCode(j, i)];
	os << "\t" << rescns[i] << "\n";
    }
}
