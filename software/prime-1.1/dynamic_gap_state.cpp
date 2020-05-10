// dynamic_gap_state.h
//
// Last Modified: 10, May 2007
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

#include "dynamic_gap_state.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#ifdef DEBUG
#include <cassert>
#endif

using namespace std;
using namespace prrn;

bool operator==(const dg& a, const dg& b)
{
    return ((a.static_gap_length_ == b.static_gap_length_)
	    && (a.dynamic_gap_length_ == b.dynamic_gap_length_));
}

void DynamicGapStates::set(const dg& rdg)
{
#ifdef DEBUG
    assert(state_.end() == find(state_.begin(), state_.end(), rdg));
#endif
    state_.push_back(rdg);
}

void DynamicGapStates::update(const GapProfile& rgp, const DynamicGapStates& rdgs)
{
    static vector<dg> tdgs;
    tdgs.resize(rgp.gv_.size()+1);
    iterator k = tdgs.begin();
    *(k++) = dg(0, 0);

    const_iterator j = rdgs.state_.begin();
    const gv_const_iter ie = rgp.gv_.end();
    for(gv_const_iter i = rgp.gv_.begin(); i != ie; ++i)
    {
	// gap profiles must be sorted
	// always find dynamic gap state
	while(i->first != j->static_gap_length_)
	    ++j;
#ifdef DEBUG
	assert(j != rdgs.state_.end());
#endif
	*(k++) = dg(j->static_gap_length_+1,
		j->dynamic_gap_length_);
    }
    state_.swap(tdgs);
}

#ifdef DEBUG
void DynamicGapStates::show() const
{
    if(state_.empty())
    {
	cerr << "emtpy\n";
	return;
    }

    const const_iterator ie = state_.end();
    for(const_iterator i = state_.begin(); i != ie; ++i)
    {
	cerr << "(" << setw(2) << i->static_gap_length_;
	cerr << ", " << setw(2) << i->dynamic_gap_length_ << ") ";
    }
    cerr << "\n";
}
#endif

cell::_cell(SCORE initscore)
    : score(initscore)
    , dgs_a()
    , dgs_b()
    , path(0)
{}

cell4lng::_cell4lng(const SCORE& initscore)
    : c(initscore)
    , dgp_a()
    , dgp_b()
{}

cell4lng::_cell4lng(const _cell4lng& rc)
    : c(rc.c)
    , dgp_a(rc.dgp_a)
    , dgp_b(rc.dgp_b)
{}

#ifdef DEBUG
void cell4lng::show() const
{
    cerr << "\ndgp_a:\n";
    if(dgp_a.empty())
	cerr << "empty";
    for(vector<pair<int, int> >::const_iterator i = dgp_a.begin(); i != dgp_a.end(); ++i)
	cerr << "(" << setw(3) << i->first << ", " << setw(3) << i->second << ") ";
    cerr << "\ndgp_b:\n";
    if(dgp_b.empty())
	cerr << "empty";
    for(vector<pair<int, int> >::const_iterator i = dgp_b.begin(); i != dgp_b.end(); ++i)
	cerr << "(" << setw(3) << i->first << ", " << setw(3) << i->second << ") ";
    cerr << "\n\n";
}
#endif
